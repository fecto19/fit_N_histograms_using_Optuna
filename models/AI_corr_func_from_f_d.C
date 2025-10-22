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
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"

using namespace std;


static double gauss_potential(double* ppars){
// ppars[0] - radius in fm
// ppars[1] - k* in MeV
// ppars[2] - V0 
// ppars[3] - \mu0 in fm
	return ppars[2] * exp(-pow(ppars[0] / ppars[3], 2));	
}

static string strip(string to_be_stripped){
    bool skip_str = false;
    vector<int> indexes_to_erase;
    string stripped = to_be_stripped;
    for(int i=0; i<to_be_stripped.size(); i++){
    if(std::isspace(to_be_stripped[i])== false){
        skip_str = true;
        }
    else if ((skip_str == false) && (std::isspace(to_be_stripped[i]))){
        indexes_to_erase.push_back(i);
        }
    }
    for(int i = 0; i<indexes_to_erase.size(); i++){
        int j = indexes_to_erase[i];
        stripped[j] = '@';
    }
    while(stripped.find('@')!= std::string::npos){
        int k = stripped.find('@');
        stripped = stripped.erase(k, 1);
        }
    return stripped;
}


static vector<double> string_to_vector_double(const string& par_string){
    string stripped_parameters = strip(par_string);
    vector<double> parameters;
    for (int i = 0; i < stripped_parameters.size(); i++){
        if ((stripped_parameters[i] == '[') || (stripped_parameters[i] == ']')){
            stripped_parameters.erase(i, 1);
        }
        while(stripped_parameters.find(',') != std::string::npos){
            int index = stripped_parameters.find(',');
            double par_i = stod(strip(stripped_parameters.substr(0, index))); // does not include comma
            stripped_parameters.erase(0, index+1); // deletes all including first comma
            parameters.push_back(par_i);
        }
        if (!stripped_parameters.empty()) {
            parameters.push_back(stod(strip(stripped_parameters)));
        }
    }
    return parameters;
}

static vector<int> string_to_vector_int(const string& par_string){
    string stripped_parameters = strip(par_string);
    vector<int> parameters;
    for (int i = 0; i < stripped_parameters.size(); i++){
        if ((stripped_parameters[i] == '[') || (stripped_parameters[i] == ']')){
            stripped_parameters.erase(i, 1);
        }
        while(stripped_parameters.find(',') != std::string::npos){
            int index = stripped_parameters.find(',');
            int par_i = stoi(strip(stripped_parameters.substr(0, index))); // does not include comma
            stripped_parameters.erase(0, index+1); // deletes all including first comma
            parameters.push_back(par_i);
        }
        if (!stripped_parameters.empty()) {
            parameters.push_back(stoi(strip(stripped_parameters)));
        }
    }
    return parameters;
}

// a, b - parameters of Gauss potential aexp{-r²/b²}

void correlation_function(double m_red, double a, double b, double source_size, double k_min, double k_max, int cpu, int n_bins){
    CATS ai_cat;
    ai_cat.SetMomBins(n_bins, k_min, k_max);
    ai_cat.SetThetaDependentSource(false);

    CATSparameters source_par(CATSparameters::tSource, 1, true);
    ai_cat.SetAnaSource(GaussSource, source_par);
    ai_cat.SetAnaSource(0, source_size);
    ai_cat.SetAutoNormSource(true);
    ai_cat.SetUseAnalyticSource(true);
    ai_cat.SetMomentumDependentSource(false);
    ai_cat.SetExcludeFailedBins(false);
    ai_cat.SetQ1Q2(0);
    ai_cat.SetQuantumStatistics(true);
    ai_cat.SetRedMass(m_red);

    CATSparameters pot_par(CATSparameters::tPotential, 2, true);
    pot_par.SetParameter(0, a);
    pot_par.SetParameter(1, b);
    ai_cat.SetNumChannels(4);
    ai_cat.SetNumPW(0,1);
    ai_cat.SetShortRangePotential(0, 0, gauss_potential, pot_par);
    ai_cat.SetSpin(0, 0);
    ai_cat.SetSpin(1, 1);
    ai_cat.SetSpin(2, 1);
    ai_cat.SetSpin(3, 1);
    ai_cat.SetChannelWeight(0, 0.25);
    ai_cat.SetChannelWeight(1, 0.25);
    ai_cat.SetChannelWeight(2, 0.25);
    ai_cat.SetChannelWeight(3, 0.25);

    ai_cat.KillTheCat();

    TH1F *corr_func_h = new TH1F("th_histo_1", "th_histo_1", n_bins, k_min, k_max);
    TH1F *phase_shifts_h = new TH1F("phase_shifts", "phase_shifts", n_bins, k_min, k_max);
    for(int i; i<n_bins; i++){
        double k_star_i = corr_func_h->GetBinCenter(i+1);
        double corr_i = ai_cat.EvalCorrFun(k_star_i);
        double phase_i = ai_cat.GetPhaseShift(i, 0, 0);
        phase_shifts_h->SetBinContent(i+1, phase_i);
        cout<<"Corr value (cpu "<<cpu<<"): "<<corr_i<<endl;
        corr_func_h->SetBinContent(i+1, corr_i);
    }
    TFile *file = new TFile(TString::Format("th_model_cpu%i_ai_corr_func_th_histo_1.root", cpu), "RECREATE");
    file->cd();
    corr_func_h->Write();
    phase_shifts_h->Write();
    file->Close();
    delete corr_func_h;
    delete phase_shifts_h;

}

// parameters passed on order: param_string, cpu_str, nbins_str, k_ranges_str

int AI_corr_main(int argc, char *argv[]){
    cout<<"Extracting parameters..."<<endl;
    string parameters_str = argv[1];
    vector<double> parameters_vec = string_to_vector_double(strip(parameters_str));
    double m_red = parameters_vec[0];
    double source_size = parameters_vec[1];
    double a = parameters_vec[2];
    double b = parameters_vec[3];
    int cpu = stoi(argv[2]);
    string n_bins_str = argv[3];
    vector<int> n_bins_vec = string_to_vector_int(strip(n_bins_str));
    int n_bins = n_bins_vec[0];
    string k_ranges_str = strip(argv[4]);
    vector<double> k_ranges_vec = string_to_vector_double(k_ranges_str);
    double k_min = k_ranges_vec[0];
    double k_max = k_ranges_vec[1];
    
    cout<<"Claculating Correlation Function..."<<endl;
    correlation_function(m_red, a, b, source_size, k_min, k_max, cpu, n_bins);
    return 0;
}
