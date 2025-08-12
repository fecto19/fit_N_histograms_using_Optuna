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

int line(double a, double b, int cpu, int nbins){
    double k_min = 0;
    double k_max = 150;
    TH1F *line_histo = new TH1F("th_histo_1", "th_histo_1", nbins, k_min, k_max);
    for(int i = 0; i<nbins; i++){
        double value_i = a * line_histo->GetBinCenter(i+1) + b;
        line_histo->SetBinContent(i+1, value_i);
    }

    TFile *line_file = new TFile(TString::Format("th_model_cpu%i_th_histo_1.root", cpu), "RECREATE");
    
    line_file->cd();
    line_histo->Write();

    line_file->Close();
    delete line_histo;
    return 0;
}

int parabole(double star, double first, double zero, int nbins, int cpu){
    double k_min = 0;
    double k_max = 150;
    TH1F *parabole_histo = new TH1F("th_histo_2", "th_histo_2", nbins, k_min, k_max);
    for(int i = 0; i<nbins; i++){
        double x_val = parabole_histo->GetBinCenter(i+1);
        double value_i = star * pow(x_val, 2) + first * x_val + zero;
        parabole_histo->SetBinContent(i+1, value_i);
        }

    TFile *parabole_file = new TFile(TString::Format("th_model_cpu%i_th_histo_2.root", cpu), "RECREATE");

    parabole_file->cd();
    parabole_histo->Write();

    parabole_file->Close();
    delete parabole_histo;

    return 0;
}

int OPTUNA_SIMPLE_TRY001(int argc, char *argv[]){
    printf("%s", "Starting simple model...");
    
    string param_string = strip(argv[1]);
    vector<double> parameters = string_to_vector_double(param_string);
    double a = parameters[0];
    double b = parameters[1];
    double star = parameters[2];
    double first = parameters[3];
    double zero = parameters[4];

    int cpu = stoi(argv[2]);

    string nbins_str = argv[3];
    vector<int> nbins_list = string_to_vector_int(nbins_str);
    int nbins_line = nbins_list[0];
    int nbins_parabole = nbins_list[1];

    line(a, b, cpu, nbins_line);
    parabole(star, first, zero, nbins_parabole, cpu);  

    return 0;
}
