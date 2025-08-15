#include "functions.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "CommonAnaFunctions.h"

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


#include <unistd.h>
#include "TH1F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"

using namespace std;

//double hbarc = 197.3269 // MeV*fm


static double gauss_potential(double* ppars){
// ppars[0] - radius in fm
// ppars[1] - k* in MeV
// ppars[2] - V0 
// ppars[3] - \mu0 in fm
	return ppars[2] * exp(-pow(ppars[0] / ppars[3], 2));
	
}

string strip(string to_be_stripped){
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


vector<double> string_to_vector_double(const string& par_string){
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

vector<int> string_to_vector_int(const string& par_string){
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


/*
int model_phase_shifts(const double f, const double d, const int cpu, const int nbins_ps){
    printf("%s\n", "creating TFile...");
    TFile *th_model_cpu_phase_shifts = new TFile(TString::Format("th_model_cpu%i_phase_shifts_th_histo_1.root", cpu), "RECREATE");
    printf("%s\n", "TFile created");
    double k_min = 0.0;
    double k_max = 150.0;
    TH1F *th_histo_1 = new TH1F("th_histo_1", "th_histo_1", nbins_ps, k_min, k_max);
    for(int i=0; i<nbins_ps; i++){
        double k_star = th_histo_1->GetBinCenter(i+1);
        double f_ = f / hbarc;
        double d_ = d / hbarc;
        double mm_val_i = 1/f_ + 0.5*d_*pow(k_star,2);
        double m_val_i = k_star / mm_val_i;
        double ps_i = atan(m_val_i);
        th_histo_1->SetBinContent(i+1, ps_i);
    }
    
    th_model_cpu_phase_shifts->cd();
    th_histo_1->Write();

    delete th_histo_1;
    printf("%s\n", "closing TFile...");
    th_model_cpu_phase_shifts->Close();
    printf("%s\n", "TFile closed :)");
    return 0;
}
*/


/*
// another model for the phase shifts since it seems it ain't working with the above one...
// this one and the one above get different parameters - must change config file
int model_phase_shifts_gauss(int nbins_ps, double a_gauss, double b_gauss, double source_size, int cpu){
    CATS gauss_cat;
    double k_min = 0;
    double k_max = 150;
    double m1 = 938.272;
    double m2 = 938.272; // must add them later from input
    gauss_cat.SetMomBins(nbins_ps, k_min, k_max);
	gauss_cat.SetThetaDependentSource(false);
	// define parameters for the source of the matching cat :)
	CATSparameters spars(CATSparameters::tSource, 1, true);
	gauss_cat.SetAnaSource(GaussSource, spars);
	gauss_cat.SetAnaSource(0, source_size);
	gauss_cat.SetAutoNormSource(true);
	gauss_cat.SetUseAnalyticSource(true);
	gauss_cat.SetMomentumDependentSource(false);
	gauss_cat.SetExcludeFailedBins(false);
	gauss_cat.SetQ1Q2(0);
	gauss_cat.SetQuantumStatistics(true);
	gauss_cat.SetRedMass((m1 * m2)/(m1 + m2));
	
	// set parameters for the potential for the matching cat
	CATSparameters ppars(CATSparameters::tPotential, 2, true);
	ppars.SetParameter(0, a_gauss);
	ppars.SetParameter(1, b_gauss);
	gauss_cat.SetNumChannels(2);
	gauss_cat.SetNumPW(0, 1);
	gauss_cat.SetShortRangePotential(0, 0, gauss_potential, ppars);
	gauss_cat.SetSpin(0, 0);
	gauss_cat.SetSpin(1, 1);
	gauss_cat.SetChannelWeight(0, 0.25);
	gauss_cat.SetChannelWeight(1, 0.75);

	gauss_cat.KillTheCat();

	TH1F *th_histo_1 = new TH1F("th_histo_1", "th_histo_1", nbins_ps, k_min, k_max);
	for(int i = 0; i<nbins_ps; i++){
	    double phase_i = gauss_cat.GetPhaseShift(i, 0, 0);
	    th_histo_1->SetBinContent(i+1, phase_i);
	}
	
	printf("%s\n", "creating TFile...");
    TFile *th_model_cpu_phase_shifts = new TFile(TString::Format("th_model_cpu%i_phase_shifts_th_histo_1.root", cpu), "RECREATE");
    printf("%s\n", "TFile created");
    
    th_model_cpu_phase_shifts->cd();
    th_histo_1->Write();

    delete th_histo_1;
    printf("%s\n", "closing TFile...");
    th_model_cpu_phase_shifts->Close();
    printf("%s\n", "TFile closed :)");
    return 0;
}
*/

int fit_f_d_from_a_b(int nbins_ps, int nbins_fd, double a_gauss, double b_gauss, double source_size, int cpu){
    nbins_fd = 10;
    double m1 = 938.272;
    double m2 = 938.272;
    double m_red = m1*m2/(m1+m2);
    double k_min = 0;
    double k_max = 80;
    
    CATS scatt_cat;
    scatt_cat.SetMomBins(nbins_fd, k_min, k_max);
    scatt_cat.SetThetaDependentSource(false);

    CATSparameters source_p(CATSparameters::tSource, 1, true);
    scatt_cat.SetAnaSource(GaussSource, source_p);
    scatt_cat.SetAnaSource(0, source_size);
    scatt_cat.SetAutoNormSource(true);
    scatt_cat.SetUseAnalyticSource(true);
    scatt_cat.SetMomentumDependentSource(false);
    scatt_cat.SetExcludeFailedBins(false);
    scatt_cat.SetQ1Q2(0);
    scatt_cat.SetQuantumStatistics(true);
    scatt_cat.SetRedMass(m_red);

    CATSparameters potential_p(CATSparameters::tPotential, 2, true);
    potential_p.SetParameter(0, a_gauss);
    potential_p.SetParameter(1, b_gauss);
    scatt_cat.SetNumChannels(2);
    scatt_cat.SetNumPW(0, 1);
    scatt_cat.SetShortRangePotential(0, 0, gauss_potential, potential_p);
    scatt_cat.SetSpin(0, 0);
    scatt_cat.SetSpin(1, 1);
    scatt_cat.SetChannelWeight(0, 0.25);
    scatt_cat.SetChannelWeight(1, 0.75);

    scatt_cat.KillTheCat();
    
    TH1F *cot_phase_shifts_h = new TH1F("cot_phase_shifts_h", "cot_phase_shifts_h", nbins_fd, k_min, k_max);
    for(int i=0; i<nbins_fd; i++){
        double phase_i = scatt_cat.GetPhaseShift(i, 0, 0);
        double c_phase_i = 1 / tan(phase_i);
        cot_phase_shifts_h->SetBinContent(i+1, c_phase_i);
        cot_phase_shifts_h->SetBinError(i+1, 0.01);
    }

    TF1 *cot_fit = new TF1("cot_fit", "[0]/x + 0.5 * [1] * x", k_min, k_max);
    cot_phase_shifts_h->Fit(cot_fit);
    double f_2B = cot_fit->GetParameter(0);
    double d_2B = cot_fit->GetParameter(1);

    double f_calc = hbarc / f_2B;
    double d_calc = d_2B * hbarc;

    TH1F *scatt_par_h = new TH1F("th_histo_2", "th_histo_2", 2, 0, 50);
    scatt_par_h->SetBinContent(1, f_calc);
    scatt_par_h->SetBinContent(2, d_calc);
    scatt_par_h->SetBinError(1, 0.01);
    scatt_par_h->SetBinError(2, 0.01);

    TH1F *phase_shifts_fd_h = new TH1F("phase_shifts_fd_h", "phase_shifts_fd_h", nbins_ps, k_min, k_max);
    for(int i=0; i<nbins_ps; i++){
        double k_star = phase_shifts_fd_h->GetBinCenter(i+1);
        double f_ = f_calc / hbarc;
        double d_ = d_calc / hbarc;
        double mm_val_i = 1/f_ + 0.5*d_*pow(k_star,2);
        double m_val_i = k_star / mm_val_i;
        double ps_i = atan(m_val_i);
        phase_shifts_fd_h->SetBinContent(i+1, ps_i);
    }

    CATS phase_cat;
    phase_cat.SetMomBins(nbins_ps, k_min, k_max);
    phase_cat.SetThetaDependentSource(false);

//    CATSparameters source_p(CATSparameters::tSource, 1, true);
    phase_cat.SetAnaSource(GaussSource, source_p);
    phase_cat.SetAnaSource(0, source_size);
    phase_cat.SetAutoNormSource(true);
    phase_cat.SetUseAnalyticSource(true);
    phase_cat.SetMomentumDependentSource(false);
    phase_cat.SetExcludeFailedBins(false);
    phase_cat.SetQ1Q2(0);
    phase_cat.SetQuantumStatistics(true);
    phase_cat.SetRedMass(m_red);

//    CATSparameters potential_pp(CATSparameters::tPotential, 2, true);
//    potential_pp.SetParameter(0, a_gauss);
//    potential_pp.SetParameter(1, b_gauss);
    phase_cat.SetNumChannels(2);
    phase_cat.SetNumPW(0, 1);
    phase_cat.SetShortRangePotential(0, 0, gauss_potential, potential_p);
    phase_cat.SetSpin(0, 0);
    phase_cat.SetSpin(1, 1);
    phase_cat.SetChannelWeight(0, 0.25);
    phase_cat.SetChannelWeight(1, 0.75);

    phase_cat.KillTheCat();

    TH1F *phase_shifts_h = new TH1F("th_histo_1", "th_histo_1", nbins_ps, k_min, k_max);
    for(int i=0; i<nbins_ps; i++){
        double phase_i = phase_cat.GetPhaseShift(i, 0, 0);
        phase_shifts_h->SetBinContent(i+1, phase_i);
        phase_shifts_h->SetBinError(i+1, 0.01);
    }


    TFile *th_model_cpu_scattering_parameters = new TFile(TString::Format("th_model_cpu%i_scattering_parameters_th_histo_2.root", cpu), "RECREATE");
    th_model_cpu_scattering_parameters->cd();
    scatt_par_h->Write();
    phase_shifts_fd_h->Write();

    delete scatt_par_h;
    delete phase_shifts_fd_h;
    th_model_cpu_scattering_parameters->Close();

    TFile *th_model_cpu_phase_shifts = new TFile(TString::Format("th_model_cpu%i_phase_shifts_th_histo_1.root", cpu), "RECREATE");
    th_model_cpu_phase_shifts->cd();
    phase_shifts_h->Write();

    delete phase_shifts_h;
    th_model_cpu_phase_shifts->Close();

    return 0;
}


int fit_only_f_d_from_a_b(int nbins_ps, double a_gauss, double b_gauss, double source_size, int cpu){
    nbins_ps = 10;
    double m1 = 938.272;
    double m2 = 938.272;
    double m_red = m1*m2/(m1+m2);
    double k_min = 0;
    double k_max = 80;
    
    CATS scatt_cat;
    scatt_cat.SetMomBins(nbins_ps, k_min, k_max);
    scatt_cat.SetThetaDependentSource(false);

    CATSparameters source_p(CATSparameters::tSource, 1, true);
    scatt_cat.SetAnaSource(GaussSource, source_p);
    scatt_cat.SetAnaSource(0, source_size);
    scatt_cat.SetAutoNormSource(true);
    scatt_cat.SetUseAnalyticSource(true);
    scatt_cat.SetMomentumDependentSource(false);
    scatt_cat.SetExcludeFailedBins(false);
    scatt_cat.SetQ1Q2(0);
    scatt_cat.SetQuantumStatistics(true);
    scatt_cat.SetRedMass(m_red);

    CATSparameters potential_p(CATSparameters::tPotential, 2, true);
    potential_p.SetParameter(0, a_gauss);
    potential_p.SetParameter(1, b_gauss);
    scatt_cat.SetNumChannels(2);
    scatt_cat.SetNumPW(0, 1);
    scatt_cat.SetShortRangePotential(0, 0, gauss_potential, potential_p);
    scatt_cat.SetSpin(0, 0);
    scatt_cat.SetSpin(1, 1);
    scatt_cat.SetChannelWeight(0, 0.25);
    scatt_cat.SetChannelWeight(1, 0.75);

    scatt_cat.KillTheCat();
    
    TH1F *cot_phase_shifts_h = new TH1F("cot_phase_shifts_h", "cot_phase_shifts_h", nbins_ps, k_min, k_max);
    for(int i=0; i<nbins_ps; i++){
        double phase_i = scatt_cat.GetPhaseShift(i, 0, 0);
        double c_phase_i = 1 / tan(phase_i);
        cot_phase_shifts_h->SetBinContent(i+1, c_phase_i);
        cot_phase_shifts_h->SetBinError(i+1, 0.01);
    }

    TF1 *cot_fit = new TF1("cot_fit", "[0]/x + 0.5 * [1] * x", k_min, k_max);
    cot_phase_shifts_h->Fit(cot_fit);
    double f_2B = cot_fit->GetParameter(0);
    double d_2B = cot_fit->GetParameter(1);

    double f_calc = hbarc / f_2B;
    double d_calc = d_2B * hbarc;

    TH1F *scatt_par_h = new TH1F("th_histo_1", "th_histo_1", 2, 0, 50);
    scatt_par_h->SetBinContent(1, f_calc);
    scatt_par_h->SetBinContent(2, d_calc);

    TH1F *phase_shifts_h = new TH1F("phase_shifts_h", "phase_shifts_h", 100, k_min, k_max);
    for(int i=0; i<100; i++){
        double k_star = phase_shifts_h->GetBinCenter(i+1);
        double f_ = f_calc / hbarc;
        double d_ = d_calc / hbarc;
        double mm_val_i = 1/f_ + 0.5*d_*pow(k_star,2);
        double m_val_i = k_star / mm_val_i;
        double ps_i = atan(m_val_i);
        phase_shifts_h->SetBinContent(i+1, ps_i);
    }


    TFile *th_model_cpu_scattering_parameters = new TFile(TString::Format("th_model_cpu%i_scattering_parameters_th_histo_1.root", cpu), "RECREATE");
    th_model_cpu_scattering_parameters->cd();
    scatt_par_h->Write();
    phase_shifts_h->Write();

    delete scatt_par_h;
    delete phase_shifts_h;
    th_model_cpu_scattering_parameters->Close();
    return 0;
}


int OPTUNA_TRY001(int argc, char *argv[]){
    printf("%s\n", "Start:");
    string parameters_str = argv[1];
    vector<double> parameters = string_to_vector_double(parameters_str);
//    double f = parameters[1];
    double a_gauss = parameters[0];
//    double d = parameters[2];
    double b_gauss = parameters[1];
    int cpu = stoi(argv[2]);
    string nbins_str = argv[3];
    vector<int> nbins = string_to_vector_int(nbins_str);
    int nbins_ps = nbins[0];
    int nbins_fd = nbins[1];
    double source_size = 1.25; // must take later from config file
//    model_phase_shifts(f, d, cpu, nbins_ps);
//    model_phase_shifts_gauss(nbins_ps, a_gauss, b_gauss, source_size, cpu);
//    fit_f_d_from_a_b(nbins_ps, nbins_fd, a_gauss, b_gauss, source_size, cpu);
    fit_only_f_d_from_a_b(nbins_fd, a_gauss, b_gauss, source_size, cpu);
    printf("%s\n", "End.");
    return 0;
}
