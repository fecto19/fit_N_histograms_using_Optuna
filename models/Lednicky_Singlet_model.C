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


static string strip(string to_be_stripped){
    bool skip_str = false;
    vector<int> indexes_to_erase;
    string stripped = to_be_stripped;
    for(int i=0; i<to_be_stripped.size(); i++){
    if(std::isspace(to_be_stripped[i]) == false){
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


int lednicky_exp(int nbins, double *source_size, double f, double d, double *range, int cpu){
    double scattering_parameters[2] = {f, d};
    TH1D *corr_func_h = new TH1D("th_histo_1", "th_histo_1)", nbins, range[0], range[1]);
    for(int i=0; i<nbins; i++){
        double k_i = corr_func_h->GetBinCenter(i+1);
        double c_k = Lednicky_Singlet(k_i, source_size, scattering_parameters);
        corr_func_h->SetBinContent(i+1, c_k);
        double error_i = 0.01;
        corr_func_h->SetBinError(i+1, error_i);
    }

    TFile *corr_func_lednicky_exp = new TFile(TString::Format("th_model_cpu%i_th_histo_1.root", cpu), "RECREATE");
    corr_func_lednicky_exp->cd();
    corr_func_h->Write();

    corr_func_lednicky_exp->Close();
    delete corr_func_h;
    
    return 0;
}


int Lednicky_Singlet_model_main(int argc, char *argv[]){
    printf("%s", "Starting Lednicky Singlet model...\n");
    
    string param_string = strip(argv[1]);
    vector<double> parameters = string_to_vector_double(param_string);
    double s_size = parameters[0];
    double f = parameters[1];
    double d = parameters[2];
    
    int cpu = stoi(argv[2]);

    string nbins_str = argv[3];
    vector<int> nbins_list = string_to_vector_int(nbins_str);
    int nbins = nbins_list[0];

    string k_ranges_str = argv[4];
    vector<double> k_ranges_list = string_to_vector_double(k_ranges_str);
    double k_min = k_ranges_list[0];
    double k_max = k_ranges_list[1];

    double source_size[1] = {s_size};

    double range[2] = {k_min, k_max};

    printf("%s%f\n", "source size", source_size[0]);
    printf("%s%f%f\n", "scattering_parameters", f, d);
    
    lednicky_exp(nbins, source_size, f, d, range, cpu);
    return 0;
}
