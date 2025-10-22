#include <TFile.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <vector>

void makeResidualHistos()
{
   
    const int nModels = 120;
    TString truthFile = "output_SEED8.txt";
    TString predPrefix = "predictions_model_";
    TString outputRoot = "Histos.root";

   
    std::ifstream fin(truthFile.Data());
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open " << truthFile << std::endl;
        return;
    }

    std::vector<double> truth1, truth2;
    std::string line;
    bool firstLine = true;

    while (std::getline(fin, line)) {
        if (firstLine) { firstLine = false; continue; } // skip header
        std::stringstream ss(line);
        double c1, c2, c3, c4, c5;
        ss >> c1 >> c2 >> c3 >> c4 >> c5;
        truth1.push_back(c2);
        truth2.push_back(c3);
    }
    fin.close();

    int nEvents = truth1.size();
    std::cout << "Loaded " << nEvents << " truth events" << std::endl;

   
    TFile *fout = new TFile(outputRoot, "RECREATE");

    for (int i = 0; i < nModels; ++i) {
        TString predFile = Form("%s%d.txt", predPrefix.Data(), i);

        std::ifstream fpred(predFile.Data());
        std::vector<double> pred1, pred2;
        double p1, p2;
        while (fpred >> p1 >> p2) {
            pred1.push_back(p1);
            pred2.push_back(p2);
        }
        fpred.close();

        TString hname = Form("hResiduals_model_%d", i);
        TString htitle = Form("Residuals (model %d);(Truth-Pred)/Truth (POT_PAR1);(Truth-Pred)/Truth (POT_PAR2)", i);
        TH2F *h = new TH2F(hname, htitle, 200, -1, 1, 200, -1, 1);

        for (int j = 0; j < nEvents; ++j) {
            double dx = (truth1[j] - pred1[j]) / truth1[j];
            double dy = (truth2[j] - pred2[j]) / truth2[j];
            h->Fill(dx, dy);
        }

        fout->cd();
        h->Write();
        std::cout << "Wrote histogram for model " << i << std::endl;
    }

    fout->Close();
    std::cout << "All histograms saved in " << outputRoot << std::endl;
}
