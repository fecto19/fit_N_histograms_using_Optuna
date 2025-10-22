#include "functions.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "CommonAnaFunctions.h"
#include "DLM_CkModels.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <vector>
#include <cstdio>
#include <complex>
#include <fstream>
#include <random>


#include <unistd.h>
#include "TH1F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"

using namespace std;

Double_t wrapper_Lednicky_Singlet(double *k_star, double *parameters){
    double source_size = parameters[0];
    double f = parameters[1];
    double d = parameters[2];
    double source_parameters[1] = {source_size};
    double scattering_parameters[2] = {f, d};
    return Lednicky_Singlet(k_star[0], source_parameters, scattering_parameters);
}

int fit_Lednicky_model(string path, string histo_name, double *ranges, int par_count, double initial_r0, double initial_f, double initial_d){
    TFile *fit_Lednicky_Singlet_TF1 = new TFile("fit_Lednicky_Singlet_TF1.root", "RECREATE");

    TFile *data_file = TFile::Open(path.c_str());
    if(!data_file){
        cerr<<"No such file! Exiting program..."<<endl;
        return 1;
    }
    else if(data_file->IsZombie()){
        cerr << "File is Zombie. Exiting..." << endl;
        return 1;
    }

    TH1F *data_histo =(TH1F*) data_file->Get(histo_name.c_str());
    if(!data_histo){
        TString warning_message = TString::Format("There is no TH1F object named %s in this file! ", histo_name.c_str());
        cerr<< warning_message << "Exiting..." << endl;
        return 1;
    }
    //data_file->Close();
    fit_Lednicky_Singlet_TF1->cd();


    TF1 *fit_using_Lednicky_Singlet = new TF1("fit_using_Lednicky_Singlet", wrapper_Lednicky_Singlet, ranges[0], ranges[1], par_count);
    fit_using_Lednicky_Singlet->SetParameters(initial_r0, initial_f, initial_d);
    fit_using_Lednicky_Singlet->FixParameter(0, initial_r0);
    fit_using_Lednicky_Singlet->FixParameter(1, initial_f);
    fit_using_Lednicky_Singlet->FixParameter(2, initial_d);
//    fit_using_Lednicky_Singlet->SetParLimits(1, 12.0, 22.0);
//    fit_using_Lednicky_Singlet->SetParLimits(2, 1.0, 3.0);
//    data_histo->Fit(fit_using_Lednicky_Singlet);
    data_histo->Fit(fit_using_Lednicky_Singlet, "S, N, R, M");
    printf("%s%f\n", "20Mev", fit_using_Lednicky_Singlet->Eval(20.0));
    double final_r0 = fit_using_Lednicky_Singlet->GetParameter(0);
    double final_f = fit_using_Lednicky_Singlet->GetParameter(1);
    double final_d = fit_using_Lednicky_Singlet->GetParameter(2);
    printf("%s%f\n", "r0 =", final_r0);
    printf("%s%f\n", "f =", final_f);
    printf("%s%f\n", "d =", final_d);
    printf("%s%f%f\n", "ranges", ranges[0], ranges[1]);

    int N = 200;
    TGraph *g = new TGraph(N);
    for (int i=0;i<N;i++) {
       double x = 0 + (ranges[1]-ranges[0])*i/(N-1);
       g->SetPoint(i, x, fit_using_Lednicky_Singlet->Eval(x));
    }

    int nbins = data_histo->GetNbinsX();
    double x_min = data_histo->GetXaxis()->GetXmin();
    double x_max = data_histo->GetXaxis()->GetXmax();
    double chi_sqr = 0.0;
    TH1F *calc_corr_h = new TH1F("calc_corr_h", "calc_corr_h", nbins, x_min, x_max);
    for(int i=0; i<nbins; i++){
        double x_i = calc_corr_h->GetBinCenter(i+1);
        double value_i = fit_using_Lednicky_Singlet->Eval(x_i);
        calc_corr_h->SetBinContent(i+1, value_i);
        double old_val_i = data_histo->GetBinContent(i+1);
        double new_val_i = calc_corr_h->GetBinContent(i+1);
        double error_i = data_histo->GetBinError(i+1);
        if (error_i <= 0){
            error_i = 1.0;
        }
        double pre_chi_i = (new_val_i - old_val_i) / error_i;
        chi_sqr = chi_sqr + pow(pre_chi_i, 2);
    }
    int ndf = fit_using_Lednicky_Singlet->GetNDF();
    printf("%s%i\n", "NDF = ", ndf);
    double norm_chi_sqr = chi_sqr / ndf;
    printf("%s%f\n", "chi_sqr (not divided by NDoF): ", chi_sqr);
    printf("%s%f\n", "normalized #chi^2 = ", norm_chi_sqr);
    printf("%s%i\n", "number of bins: ", calc_corr_h->GetNbinsX());

    cout << "ROOT raw chi2 = " << fit_using_Lednicky_Singlet->GetChisquare() << endl;
    cout << "ROOT NDF = " << fit_using_Lednicky_Singlet->GetNDF() << endl;
    cout << "ROOT chi2/NDF = " << fit_using_Lednicky_Singlet->GetChisquare()/fit_using_Lednicky_Singlet->GetNDF() << endl;
    


    fit_Lednicky_Singlet_TF1->cd();
    data_histo->Write();
    calc_corr_h->Write();
    g->Write();
    fit_using_Lednicky_Singlet->Write();
//    fit_Lednicky_Singlet_TF1->Close();
    
    return 0;
}


int fit_Lednicky_model_main(int argc, char *argv[]){
    string path = argv[1];
    string histo_name = argv[2];
    double k_min = stod(argv[3]);
    double k_max = stod(argv[4]);
    int par_count = stoi(argv[5]);
    double initial_r0 = stod(argv[6]);
    double initial_f = stod(argv[7]);
    double initial_d = stod(argv[8]);
    double ranges[2] = {k_min, k_max};
    fit_Lednicky_model(path, histo_name, ranges, par_count, initial_r0, initial_f, initial_d);
    return 0;
}
