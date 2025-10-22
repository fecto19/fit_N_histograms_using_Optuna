#include "functions.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CATStools.h"
#include "CATSconstants.h"
//#include "CommonAnaFunctions.h"

#include "/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/CATS/DLM/DLM_FemtoTools/CommonAnaFunctions.h"


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

//double hbarc = 197.327;


int av18_no_coloumb_main(){
// set mass for the interacting particles in MeV
	double m1 = 938;
	double m2 = 938;
	double source_size = 1.5;

	double k_star_min = 0;
	double k_star_max = 200;
	int bin_number = 100;
	int bins_for_r = 1000;
	double r_star_min = 0;
	double r_star_max = 10;

	CATS smelly_cat;
// set how many points, k_min, k-max
	smelly_cat.SetMomBins(bin_number, k_star_min, k_star_max);

// does set up of the source and potential parameters, ...
    DLM_CommonAnaFunctions AnalysisObject; 
	AnalysisObject.SetUpCats_pp(smelly_cat,"AV18","Gauss",0,0);
	smelly_cat.SetAnaSource(0, source_size);

//commented because they are set up with the right values by default
//	smelly_cat.SetChannelWeight(0,0.25);
//	smelly_cat.SetChannelWeight(1,0.25);
//	smelly_cat.SetChannelWeight(2,0.25);
//	smelly_cat.SetChannelWeight(3,0.25);
	
// deleting d wave:
	smelly_cat.RemoveShortRangePotential(0,2);
// deleting p wave:
	smelly_cat.RemoveShortRangePotential(1,1);
	smelly_cat.RemoveShortRangePotential(2,1);
	smelly_cat.RemoveShortRangePotential(3,1);

	smelly_cat.SetQ1Q2(1);
	smelly_cat.SetQuantumStatistics(true);

	smelly_cat.SetChannelWeight(0, 0.25);
    smelly_cat.SetChannelWeight(1, 0.25);
    smelly_cat.SetChannelWeight(2, 0.25);
    smelly_cat.SetChannelWeight(3, 0.25);	

	smelly_cat.KillTheCat();

	TH1F *phase_shifts_av18_no_coloumb_h = new TH1F("phase_shifts_av18_no_coloumb_h", "phase_shifts_av18_no_coloumb_h", bin_number, k_star_min, k_star_max);
	TH1F *cot_ps_av18_no_coloumb_h = new TH1F("cot_ps_av18_no_coloumb_h", "cot_ps_av18_no_coloumb_h", bin_number, k_star_min, k_star_max);
	TH1F *f_d_av18_no_coloumb_h = new TH1F("f_d_av18_no_coloumb_h", "f_d_av18_no_coloumb_h", 2, 0, 50);
	TH1F *calc_ps_h = new TH1F("calc_ps_h", "calc_ps_h", bin_number, k_star_min, k_star_max);
	for(int i=0; i<bin_number; i++){
	    double phase_i = smelly_cat.GetPhaseShift(i, 0, 0);
	    phase_shifts_av18_no_coloumb_h->SetBinContent(i+1, phase_i);
	    phase_shifts_av18_no_coloumb_h->SetBinError(i+1, 0.01);
	    double cot_i = 1.0 / tan(phase_i);
	    cot_ps_av18_no_coloumb_h->SetBinContent(i+1, cot_i);
	    cot_ps_av18_no_coloumb_h->SetBinError(i+1, 0.01);
	}
	TF1 *fit = new TF1("fit", "[0] / x + 0.5 * [1] * x", k_star_min, k_star_max);
	cot_ps_av18_no_coloumb_h->Fit(fit, "R");
    double f_2b = fit->GetParameter(0);
    double d_2B = fit->GetParameter(1);
    double f = (1 / f_2b) * hbarc;
    double d = d_2B * hbarc;

    f_d_av18_no_coloumb_h->SetBinContent(1, f);
    f_d_av18_no_coloumb_h->SetBinContent(2, d);

    for(int i =0; i<bin_number; i++){
        double k_st = calc_ps_h->GetBinCenter(i+1);
        double value = f_2b + k_st * k_st * 0.5 * d_2B;
        double phase_i = atan(k_st / value);
        calc_ps_h->SetBinContent(i+1, phase_i);
        calc_ps_h->SetBinError(i+1, 0.01);
    }

    TH1F *corr_func_h = new TH1F("corr_func_h", "corr_func_h", bin_number, k_star_min, k_star_max);
    for(int i=0; i<bin_number; i++){
        double k_star_i = corr_func_h->GetBinCenter(i+1);
        double corr_i = smelly_cat.EvalCorrFun(k_star_i);
        corr_func_h->SetBinError(i+1, 0.01);
        corr_func_h->SetBinContent(i+1, corr_i);
    }

	TFile *av18_no_coloumb = new TFile(TString::Format("av18_no_coloumb_%f.root", source_size), "RECREATE");
    av18_no_coloumb->cd();
    phase_shifts_av18_no_coloumb_h->Write();
    f_d_av18_no_coloumb_h->Write();
    corr_func_h->Write();
    calc_ps_h->Write();

    delete phase_shifts_av18_no_coloumb_h;
    delete f_d_av18_no_coloumb_h;
    delete calc_ps_h;
    delete corr_func_h;
    av18_no_coloumb->Close();

    return 0;
}
