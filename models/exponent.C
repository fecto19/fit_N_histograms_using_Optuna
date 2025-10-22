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


cout<<"create_exponent(int nbins, double x_min, double x_max)"<<endl;

int create_exponent(int nbins, double k_min, double k_max){
    TH1F *exponent_h = new TH1F("exponent_h", "exponent_h", nbins, k_min, k_max);
    for(int i; i<nbins; i++){
        double x_i = exponent_h->GetBinCenter(i+1);
        double value_i = 4.0*exp(-2.3 * x_i);
    }
    return 0;
}


int fit_exponent(){
    TFile exponent_data = TFile::Open("");
    return 0;
}
