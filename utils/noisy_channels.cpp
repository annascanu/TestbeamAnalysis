#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>

// c++ noisy_channels.cpp `root-config --cflags --libs` -o noise.out

int main() 
{
    TString filename = "../Testbeam2025/sampic_run22_final.root";
    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) 
    {
        std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
        return 1;
    }
    TTree *tree = (TTree*) f->Get("picoTree");
    if (!tree) 
    {
        std::cerr << "Errore: TTree 'picoTree' non trovato" << std::endl;
        return 1;
    }

    const int MAXPULSES = 200;
    Int_t npulses;
    Int_t Board[MAXPULSES];      // numero del detector (0, 1, 2)
    Int_t Channel[MAXPULSES];    // numero del canale (0 – 63)
    Int_t HitFeb[3];             // numero di hit per ciascun detector
    Bool_t pulses_bad_pulse[MAXPULSES];

    tree->SetBranchAddress("npulses", &npulses);
    tree->SetBranchAddress("Board", Board);
    tree->SetBranchAddress("Channel", Channel);
    tree->SetBranchAddress("HitFeb", HitFeb);
    tree->SetBranchAddress("pulses_bad_pulse", pulses_bad_pulse);

    const int TARGET_BOARD = 0;   // detector 0 → 64 canali
    const int NCHANNELS = 64;
    int singleHitCount[NCHANNELS] = {0};

    Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; i++) 
    {
        tree->GetEntry(i);

        int goodPulse = -1; // Trovo il singolo pulse valido nel detector scelto

        for (int j = 0; j < npulses; j++) 
        {
            if (HitFeb[0] == 1 && HitFeb[1] == 0 && HitFeb[2] == 0) // Just single hit events aka potential noise?
            {
                if (pulses_bad_pulse[j]) 
                    continue;      // scarta pulse cattivi
                if (Board[j] != TARGET_BOARD) 
                    continue; // scarta altri detectors

                if (goodPulse == -1) 
                    goodPulse = j;
                else 
                { 
                    goodPulse = -1; 
                    break; 
                } // più di un pulse valido → skip
                
                if (goodPulse < 0) 
                    continue;

                int ch = Channel[goodPulse];

                if (ch >= 0 && ch < NCHANNELS)
                    singleHitCount[ch]++;
            }
        }
    }
    TH1F *h = new TH1F("hNoisy", Form("Single-hit events per channel (Detector %d);Channel;Counts", TARGET_BOARD), NCHANNELS, 0, NCHANNELS);

    for (int ch = 0; ch < NCHANNELS; ch++)
        h->SetBinContent(ch + 1, singleHitCount[ch]);

    TCanvas *c = new TCanvas("c", "Noisy Channels", 900, 600);
    h->SetFillColor(kOrange-3);
    h->SetLineColor(kRed+1);
    h->SetLineWidth(2);
    h->Draw("HIST");

    c->SaveAs("noisy_channels_detector0.pdf");
    TFile fout("noisy_channels_detector0.root", "RECREATE");
    h->Write();
    c->Write();
    fout.Close();

    std::cout << "Fatto! Salvato: noisy_channels_detector0.pdf" << std::endl;

    return 0;
}