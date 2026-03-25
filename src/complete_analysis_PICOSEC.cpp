/*Goal: obtain time resolution of picosec detectors.
To compile: 
c++ analysis_picosec.cpp analysis_picosec.cc `root-config --cflags --libs` -o analysis.out*/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TProfile.h>
#include <TF1.h>
#include <vector>
#include <string>
#include <iostream>
#include <TString.h>






int main(){
       // Apri il file ROOT
    TFile *file = TFile::Open("/home/riccardo-speziali/after_waveforms_analysis/sampic_mcp_setup_run222_final.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Errore nell'apertura del file ROOT" << std::endl;
        return 1;
    }

    // Ottieni l'albero
    TTree *tree = (TTree*)file->Get("picoTree");

    // Dichiarazione delle variabili per i branch
    std::vector<double> *Cell0TimeStamp_PICOSEC = nullptr;
    std::vector<double> *cfd_PICO = nullptr;
    std::vector<double> *amplitude_PICOSEC = nullptr;
    std::vector<double> *e_peak_PICOSEC = nullptr;
    
    double Cell0TimeStamp_MCP;
    double cfd_MCP;
    int hits;
    int channel;
    int srs;
    double track_chi_2;

    // Collegamento dei branch
    tree->SetBranchAddress("Cell0timeSTamp_PICOSEC", &Cell0TimeStamp_PICOSEC);
    tree->SetBranchAddress("Cell0timestamp_MCP", &Cell0TimeStamp_MCP);
    tree->SetBranchAddress("pulses_time_cfd30", &cfd_PICO);
    tree->SetBranchAddress("mcp_time_cfd30", &cfd_MCP);
    tree->SetBranchAddress("hit_x_event", &hits);
    tree->SetBranchAddress("chanel_PICOSEC", &channel);
    tree->SetBranchAddress("TriggerIDSRS_MCP", &srs);
    tree->SetBranchAddress("pulses_amplitude", &amplitude_PICOSEC);
    tree->SetBranchAddress("pulses_integral", &e_peak_PICOSEC);
    tree->SetBranchAddress("track_chi_2", &track_chi_2);

    // Loop sugli eventi
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Esempio: stampare i primi valori per controllo
        if (Cell0TimeStamp_PICOSEC->size() > 0) {
            std::cout << "Evento " << i 
                      << ": Cell0TimeStamp_PICOSEC[0] = " << (*Cell0TimeStamp_PICOSEC)[0] 
                      << ", Cell0TimeStamp_MCP = " << Cell0TimeStamp_MCP << std::endl;
        }
    }


    // Istogrammi
    TH1D *h_chi2 = new TH1D("h_chi2", "Distribution of Track Chi^2", 50, 0, 10);
    TH1D *h_chi2_cut = new TH1D("h_chi2_cut", "Distribution of Track Chi^2 after cut", 50, 0, 10);

    // Loop eventi
    Long64_t nentries = tree->GetEntries();

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Rimuovi NaN
        if (std::isnan(track_chi_2)) continue;

        // Istogramma prima del cut
        h_chi2->Fill(track_chi_2);

        // Applica cut chi2 < 3
        if (track_chi_2 < 3) {
            h_chi2_cut->Fill(track_chi_2);
            
        }
    }

    // Canvas 1: prima del cut
    TCanvas *c1 = new TCanvas("c1", "Chi2", 800, 600);
    h_chi2->GetXaxis()->SetTitle("Track Chi^2");
    h_chi2->GetYaxis()->SetTitle("Frequency");
    h_chi2->SetFillColor(kBlue);
    h_chi2->Draw();

    // Canvas 2: dopo cut
    TCanvas *c2 = new TCanvas("c2", "Chi2 after cut", 800, 600);
    h_chi2_cut->GetXaxis()->SetTitle("Track Chi^2");
    h_chi2_cut->GetYaxis()->SetTitle("Frequency");
    h_chi2_cut->SetFillColor(kGreen);
    h_chi2_cut->Draw();


    file->Close();
    delete file;

    return 0;


}

