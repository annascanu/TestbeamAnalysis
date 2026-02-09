/*
Analysis for the October 2025 ENUBET testbeam picosec time resolution.
Authors: ...
Goal: obtain time resolution of picosec detectors.
To compile: c++ analysis_runs.cpp `root-config --cflags --libs` -o analysis_runs.out

Roadmap:
- [x] work on amplitude plotting
- [] understand what threshold to use for t_0 determination
- [] understand how to deal with charge sharing
- [] ...
*/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

int main() 
{
    // gROOT->SetBatch(kFALSE); 
    gStyle->SetOptStat(0);

    // -------------------------------------------------
    //     Load run file with processed waveforms
    // -------------------------------------------------
    /*
    Processed files can be found at this path:
    /eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/
    - Muon runs: run 19, 20, 21, 22
    - 15 GeV pion runs: run 23
    */

    TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root";
    cout << "Opening file: " << filename << endl;

    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: cannot open " << filename << endl;
        return 1;
    }
    TTree *tree = (TTree*)file->Get("picoTree");
    if (!tree) {
        cerr << "Error: picoTree not found in " << filename << endl;
        return 1;
    }

    // -----------------------
    //    Set up branches
    // -----------------------
    const int MAXPULSES = 200; // Maybe could be less
    Int_t npulses;
    Double_t pulses_amplitude[MAXPULSES], pulses_integral[MAXPULSES], pulses_time_cfd10[MAXPULSES], pulses_time_cfd20[MAXPULSES], pulses_time_cfd30[MAXPULSES];
    Bool_t  pulses_bad_pulse[MAXPULSES];

    tree->SetBranchAddress("npulses", &npulses);
    tree->SetBranchAddress("pulses_amplitude", &pulses_amplitude);
    tree->SetBranchAddress("pulses_integral", &pulses_integral);
    tree->SetBranchAddress("pulses_time_cfd10", &pulses_time_cfd10);
    tree->SetBranchAddress("pulses_time_cfd20", &pulses_time_cfd20);
    tree->SetBranchAddress("pulses_time_cfd20", &pulses_time_cfd30);
    tree->SetBranchAddress("pulses_bad_pulse", &pulses_bad_pulse);

    // ------------------------------
    //     Initialize histograms
    // ------------------------------
    TH1F *hAmpAll  = new TH1F("hAmpAll", "All pulse amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hAmp1Hit = new TH1F("hAmp1Hit", "Single-Hit event amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hQ       = new TH1F("hQ", "Pulse integral (namely charge);Integral [a.u.];Counts", 1000, 0, 10);
    TH1F *hTime    = new TH1F("hTime", "Constant fraction discriminator set at 20% Time;Time [ns];Counts", 1000, 0, 10);

    // -----------------------
    //        Read tree
    // -----------------------
    Long64_t nentries = tree->GetEntries();
    // cout << "Total entries: " << nentries << endl;

    double ampThreshold = 50.0; // 
    int nGood, idx;

    for (Long64_t i = 0; i < nentries; i++) 
    {
        tree->GetEntry(i);

        nGood = 0;
        idx = -1;

        for (int j = 0; j < npulses; j++) // Loop inside the event
        {
            if (pulses_bad_pulse[j])
                continue;                             // Skip events where the fit didn't work
            //if (pulses_amplitude[j] < ampThreshold) 
            //    continue;                              

            nGood++;
            idx = j;
            hAmpAll->Fill(pulses_amplitude[j]);
        }

        if (nGood == 1 && idx >= 0) // Only one hit
        {
            hAmp1Hit->Fill(pulses_amplitude[idx]);
            hQ->Fill(pulses_integral[idx]);
            hTime->Fill(pulses_time_cfd20[idx]);
        }

        if (i % 100000 == 0) cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }
    cout << endl;

    // -----------------------
    //       Plot results
    // -----------------------
    TCanvas *c1 = new TCanvas("c1", "Amplitude Comparison", 800, 600);
    hAmpAll->SetLineColor(kBlue);
    hAmp1Hit->SetLineColor(kRed);
    hAmpAll->Draw("HIST");
    hAmp1Hit->Draw("HIST SAME");

    auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->AddEntry(hAmpAll, "All events", "l");
    legend->AddEntry(hAmp1Hit, "Single-hit events", "l");
    legend->Draw();

    // -----------------------
    //     Save everything
    // -----------------------
    TFile *fout = new TFile("muon_run22_analysis.root", "RECREATE");
    hAmpAll->Write();
    hAmp1Hit->Write();
    hQ->Write();
    hTime->Write();
    c1->Write();
    fout->Close();

    cout << "Analysis complete. Results saved to muon_run22_analysis.root" << endl;

    file->Close();
    return 0;
}