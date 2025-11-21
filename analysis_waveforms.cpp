/*
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
                Analysis for the October 2025 ENUBET testbeam picosec time resolution.
                               Authors: A. Scanu, R. Speziali
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

Goal: obtain time resolution of picosec detectors.
To compile: c++ analysis_waveforms.cpp `root-config --cflags --libs` -o analysis.out
Roadmap:
- [x] work on amplitude plotting
- [] understand what threshold to use for t_0 determination (ongoing)
- [] understand how to deal with charge sharing (not started)
- [] ...

Things we are trying to understand or need to ask:
- Q: Order of pulse_time and pulse_amplitude inside the array (probably understood: should use Board to discriminate?)
  A: Tentative answer is to try and use the Board array to identify which pulse corresponds to which FEB. What I (Anna) don't understand is if 
- Q: how to take into account the RMS baseline value for the CDF threshold?
  A: (still in progress)
*/

/*
IMPORTANT NOTE ABOUT FILE PATHS:
Processed files can be found at this path:
/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/
- Muon runs: run 19, 20, 21, 22
- 15 GeV pion runs: run 23, ...
*/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <TStyle.h>
#include <TLatex.h>
#include <TF1.h>
#include "TMath.h"

using namespace std;

int main()
{
    // gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(0);

    // -------------------------------------------------
    //     Load run file with processed waveforms
    // -------------------------------------------------

    // TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root"; // Filename while running on lxplus
    TString filename = "/Users/anna/Developing/PhD/Testbeam2025/sampic_run22_final.root"; // Filename while running on Anna's machine
    // TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run19_final.root"; // Filename while running on Riccardo's machine
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
    const int MAXPULSES = 200; // Maybe could be less.
    Int_t npulses;

    Double_t pulses_amplitude[MAXPULSES], pulses_integral[MAXPULSES], pulses_time_cfd10[MAXPULSES], pulses_time_cfd20[MAXPULSES], pulses_time_cfd30[MAXPULSES];
    Bool_t  pulses_bad_pulse[MAXPULSES];
    Int_t HitFeb[3], Board[MAXPULSES];
    int tot_hit_feb;
    Double_t time_12, time_23, time_13, time20_12, time20_13, time20_23;
    Double_t time30_12, time30_13, time30_23;

    tree->SetBranchAddress("npulses", &npulses);
    tree->SetBranchAddress("pulses_amplitude", &pulses_amplitude);
    tree->SetBranchAddress("pulses_integral", &pulses_integral);
    tree->SetBranchAddress("pulses_time_cfd10", &pulses_time_cfd10);
    tree->SetBranchAddress("pulses_time_cfd20", &pulses_time_cfd20);
    tree->SetBranchAddress("pulses_time_cfd30", &pulses_time_cfd30);
    tree->SetBranchAddress("pulses_bad_pulse", &pulses_bad_pulse);
    tree->SetBranchAddress("HitFeb", &HitFeb);
    tree->SetBranchAddress("Board", &Board);

    // ------------------------------
    //     Initialize histograms
    // ------------------------------
    TH1F *hAmpAll  = new TH1F("hAmpAll", "All pulse amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hAmp1Hit = new TH1F("hAmp1Hit", "Single-Hit event amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hQ       = new TH1F("hQ", "Pulse integral (namely charge);Integral [a.u.];Counts", 1000, 0, 10);
    TH1F *hTime    = new TH1F("hTime", "Constant fraction discriminator set at 20% Time;Time [ps];Counts", 1000, 0, 10);

    // Histograms for triple hits
    TH1F *htriple1  = new TH1F("htriple1", "amplitude of a triple", 1000, 0, 500);
    TH1F *htriple2  = new TH1F("htriple2", "amplitude of a triple", 1000, 0, 500);
    TH1F *htriple3  = new TH1F("htriple3", "amplitude of a triple", 1000, 0, 500);

    // Histograms for time differences for the different values of CDF
    TH1F *htime_12    = new TH1F("htime_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime_13    = new TH1F("htime_13", "time diff ;Time [ps];Counts", 100, -2000, 1000);
    TH1F *htime_23    = new TH1F("htime_23", "time diff ;Time [ps];Counts", 100, -2000, 2000);

    TH1F *htime20_12    = new TH1F("htime20_12", "time diff ;Time [ps];Counts", 100, -100, 100);
    TH1F *htime20_13    = new TH1F("htime20_13", "time diff ;Time [ps];Counts", 100, -100, 100);
    TH1F *htime20_23    = new TH1F("htime20_23", "time diff ;Time [ps];Counts", 100, -100, 100);

    TH1F *htime30_12    = new TH1F("htime30_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime30_13    = new TH1F("htime30_13", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime30_23    = new TH1F("htime30_23", "time diff ;Time [ps];Counts", 100, -2000, 2000);

    // -----------------------
    //        Read tree
    // -----------------------
    Long64_t nentries = tree->GetEntries();
    int nGood, idx;

    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        nGood = 0;
        idx = -1;

        for (int j = 0; j < npulses; j++) // Loop inside the event
        {
            if (pulses_bad_pulse[j])
                continue; // Skip events where the fit didn't work

            nGood++;
            idx = j;
            hAmpAll->Fill(pulses_amplitude[j]);
        }

        if (nGood == 1 && idx >= 0)
        {
            hAmp1Hit->Fill(pulses_amplitude[idx]);
            hQ->Fill(pulses_integral[idx]);
            hTime->Fill(pulses_time_cfd20[idx]);
        }

        // ---- Triple hits for all three FEBs ----
        // AS: Replaced assumption about pulse index ordering with explicit Board[] mapping
        // also see the hitfeb = [1, 1, anything] (otherwise we lose statistics because of smaller geometrical acceptance of 3rd detector)
        if (HitFeb[0] == 1 && HitFeb[1] == 1 && nGood == 2)
        {
            // Arrays indexed by FEB number (0,1,2)
            double ampFEB[3]  = {-1.0, -1.0, -1.0};
            double cfd10[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd20[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd30[3]   = {-9999.0, -9999.0, -9999.0};

            // Fill according to Board[]
            for (int j = 0; j < npulses; j++)
            {
                if (pulses_bad_pulse[j])
                    continue;
                int det = Board[j];   // Which detector this pulse belongs to (either 0, 1, 2)

                if (det < 0 || det > 2)
                {
                    cout << "WARNING: unexpected Board index = " << det << " (event " << i << ", pulse " << j << ")" << endl;
                    continue;
                }

                ampFEB[det] = pulses_amplitude[j];
                cfd10[det]  = pulses_time_cfd10[j];
                cfd20[det]  = pulses_time_cfd20[j];
                cfd30[det]  = pulses_time_cfd30[j];
            }

            /* // This is probably not needed anymore
            // Sanity check to make sure we found valid entries for all three FEBs
            bool validTriple = true;
            for (int d = 0; d < 3; ++d)
                if (ampFEB[d] < 0) 
                    validTriple = false;
    
            if (!validTriple) 
                continue; // skip this event if something weird happened
            */

            // Fill amplitude histograms           // AS: why this x 1000 scaling?
            htriple1->Fill(ampFEB[0] * 1000);
            htriple2->Fill(ampFEB[1] * 1000);
            htriple3->Fill(ampFEB[2] * 1000);

            // CFD = 10%
            time_12 = cfd10[0] - cfd10[1];
            time_23 = cfd10[1] - cfd10[2];
            time_13 = cfd10[0] - cfd10[2];

            htime_12->Fill(time_12);
            htime_23->Fill(time_23);
            htime_13->Fill(time_13);

            // CFD = 20%
            time20_12 = cfd20[0] - cfd20[1];
            time20_23 = cfd20[1] - cfd20[2];
            time20_13 = cfd20[0] - cfd20[2];

            htime20_12->Fill(time20_12);
            htime20_23->Fill(time20_23);
            htime20_13->Fill(time20_13);

            // CFD = 30%
            time30_12 = cfd30[0] - cfd30[1];
            time30_23 = cfd30[1] - cfd30[2];
            time30_13 = cfd30[0] - cfd30[2];

            htime30_12->Fill(time30_12);
            htime30_23->Fill(time30_23);
            htime30_13->Fill(time30_13);
        }

        // Progress indicator for sanity :)
        if (i % 100000 == 0) cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }
    
    /// --- TF1 Polya function ---
    TF1 *fpolyaAmpl1 = new TF1("fpolyaAmpl1", "([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])", 0, 800);
    fpolyaAmpl1->SetTitle("Polya fit (ampl.)");
    fpolyaAmpl1->SetParName(0,"c");
    fpolyaAmpl1->SetParName(1,"#bar{Q}");
    fpolyaAmpl1->SetParName(2,"#theta");
    fpolyaAmpl1->SetLineColor(kGreen+2);
    //fpolyaAmpl1->SetRange(minax, amplMax[ci]);
    fpolyaAmpl1->SetNpx(10000);
    //fpolyaAmpl1->SetParameter(0,800);
    fpolyaAmpl1->SetParameter(1,180);
    //fpolyaAmpl1->SetParameter(2,30);
    //fpolyaAmpl1->SetParLimits(2,0.001,5000);
    fpolyaAmpl1->SetRange(0, 800);

    TF1 *fpolyaAmpl2 = (TF1*)fpolyaAmpl1->Clone("fpolyaAmpl2");
    TF1 *fpolyaAmpl3 = (TF1*)fpolyaAmpl1->Clone("fpolyaAmpl3");

    // --- Tre canvas separati ---
    // Canvas 1
    // AS: removed grids, I think the plots look nicer, but can be re-added if needed
    TCanvas *c1 = new TCanvas("c1", "htriple1", 900, 700);
    htriple1->SetLineColor(kAzure+1);
    htriple1->SetLineWidth(3);
    htriple1->SetFillColorAlpha(kAzure-4, 0.35);
    htriple1->SetTitle("Amplitudes of waveforms on the first FEB");
    htriple1->Draw("HIST");
    htriple1->Fit(fpolyaAmpl1,"R"); // Fit su htriple1
    fpolyaAmpl1->Draw("SAME");
    c1->Update();
    c1->SaveAs("hist1.pdf");

    // Canvas 2
    TCanvas *c2 = new TCanvas("c2", "htriple2", 900, 700);
    htriple2->SetLineColor(kGreen+1);
    htriple2->SetLineWidth(3);
    htriple2->SetFillColorAlpha(kGreen-4, 0.35);
    htriple2->SetTitle("Amplitudes of waveforms on the second FEB");
    htriple2->Draw("HIST");
    htriple2->Fit(fpolyaAmpl2,"R"); // Fit su htriple2
    fpolyaAmpl2->Draw("SAME");
    c2->Update();
    c2->SaveAs("hist2.pdf");

    // Canvas 3
    TCanvas *c3 = new TCanvas("c3", "htriple3", 800, 600);
    htriple3->SetLineColor(kRed+1);
    htriple3->SetLineWidth(3);
    htriple3->SetFillColorAlpha(kRed-4, 0.35);
    htriple3->SetTitle("Amplitudes of waveforms on the third FEB");
    htriple3->Draw("HIST");
    htriple3->Fit(fpolyaAmpl3,"R"); // Fit su htriple3
    fpolyaAmpl3->Draw("SAME");
    c3->Update();
    c3->SaveAs("hist3.pdf");

    // --- Canvas unico con tutti e tre ---
    TCanvas *cAll = new TCanvas("cAll", "All Histograms", 900, 700);

    htriple1->SetLineColor(kAzure+1);
    htriple1->SetLineWidth(3);
    htriple1->SetFillColorAlpha(kAzure-4, 0.35);
    htriple1->Draw("HIST");
    htriple1->Fit(fpolyaAmpl1,"R SAME");

    htriple2->SetLineColor(kGreen+1);
    htriple2->SetLineWidth(3);
    htriple2->SetFillColorAlpha(kGreen-4, 0.35);
    htriple2->Draw("HIST SAME");
    htriple2->Fit(fpolyaAmpl2,"R SAME");

    htriple3->SetLineColor(kRed+1);
    htriple3->SetLineWidth(3);
    htriple3->SetFillColorAlpha(kRed-4, 0.35);
    htriple3->Draw("HIST SAME");
    htriple3->Fit(fpolyaAmpl3,"R SAME");

    cAll->Update();
    cAll->SaveAs("histAll.pdf");

    // ------------ Gaussian plot of the time difference ---------------

    TCanvas *ctime12 = new TCanvas("ctime12", "htriple1", 900, 700);
    htime20_12->SetLineColor(kAzure+1);
    htime20_12->SetLineWidth(3);
    htime20_12->SetFillColorAlpha(kAzure-4, 0.35);
    htime20_12->SetTitle("Amplitude");
    htime20_12->Draw("HIST");
    htime20_12->Fit("gaus"); // Fit su htime20_12
    ctime12->Update();

    TF1 *f = htime20_12->GetFunction("gaus");
    double amp_12   = f->GetParameter(0);  ///taking gaus parameters
    double mean_12  = f->GetParameter(1);
    double sigma_12 = f->GetParameter(2);

    double e_amp_12   = f->GetParError(0);
    double e_mean_12  = f->GetParError(1);
    double e_sigma_12 = f->GetParError(2);

    TCanvas *ctime23 = new TCanvas("ctime23", "htriple1", 900, 700);
    htime20_23->SetLineColor(kAzure+1);
    htime20_23->SetLineWidth(3);
    htime20_23->SetFillColorAlpha(kAzure-4, 0.35);
    htime20_23->SetTitle("Amplitude");
    htime20_23->Draw("HIST");
    htime20_23->Fit("gaus"); // Fit su htime20_23
    ctime23->Update();

    TF1 *f2 = htime20_23->GetFunction("gaus");
    double amp_23   = f2->GetParameter(0);  ///taking gaus parameters
    double mean_23  = f2->GetParameter(1);
    double sigma_23 = f2->GetParameter(2);

    double e_amp_23   = f2->GetParError(0);
    double e_mean_23  = f2->GetParError(1);
    double e_sigma_23 = f2->GetParError(2);

    TCanvas *ctime13 = new TCanvas("ctime13", "htriple1", 900, 700);
    htime20_13->SetLineColor(kAzure+1);
    htime20_13->SetLineWidth(3);
    htime20_13->SetFillColorAlpha(kAzure-4, 0.35);
    htime20_13->SetTitle("Amplitude");
    htime20_13->Draw("HIST");
    htime20_13->Fit("gaus"); // Fit su htime20_13
    ctime13->Update();

    // Take gaussian from the correct histogram htime20_13
    TF1 *f3 = htime20_13->GetFunction("gaus");
    double amp_13   = f3->GetParameter(0);  // taking gaus parameters
    double mean_13  = f3->GetParameter(1);
    double sigma_13 = f3->GetParameter(2);

    double e_amp_13   = f3->GetParError(0);
    double e_mean_13  = f3->GetParError(1);
    double e_sigma_13 = f3->GetParError(2);

    // --------------- Time resolution (with error propagation) ---------------

    // Use: sigma_12^2 = sigma1^2 + sigma2^2, etc.
    // Inversion:
    // sigma1^2 = 1/2 * ( sigma_12^2 + sigma_13^2 - sigma_23^2 )
    // sigma2^2 = 1/2 * ( sigma_12^2 + sigma_23^2 - sigma_13^2 )
    // sigma3^2 = 1/2 * ( sigma_13^2 + sigma_23^2 - sigma_12^2 )

    double s12 = sigma_12;
    double s13 = sigma_13;
    double s23 = sigma_23;

    double es12 = e_sigma_12;
    double es13 = e_sigma_13;
    double es23 = e_sigma_23;

    double s1sq = 0.5 * ( s12*s12 + s13*s13 - s23*s23 );
    double s2sq = 0.5 * ( s12*s12 + s23*s23 - s13*s13 );
    double s3sq = 0.5 * ( s13*s13 + s23*s23 - s12*s12 );

    if (s1sq < 0) {
        cout << "WARNING: sigma1^2 < 0 (" << s1sq << "). Setting to 0 to avoid NaN." << endl;
        s1sq = 0;
    }
    if (s2sq < 0) {
        cout << "WARNING: sigma2^2 < 0 (" << s2sq << "). Setting to 0 to avoid NaN." << endl;
        s2sq = 0;
    }
    if (s3sq < 0) {
        cout << "WARNING: sigma3^2 < 0 (" << s3sq << "). Setting to 0 to avoid NaN." << endl;
        s3sq = 0;
    }

    double sigma_1 = TMath::Sqrt(s1sq);
    double sigma_2 = TMath::Sqrt(s2sq);
    double sigma_3 = TMath::Sqrt(s3sq);

    // Error propagation:
    // S1 = sigma1^2 = 1/2*(a + b - c), with a = s12^2, b = s13^2, c = s23^2
    // var(a) = (2*s12*es12)^2, etc.
    // var(S1) = (1/2)^2 * ( var(a) + var(b) + var(c) )
    // => var(S1) = s12^2*es12^2 + s13^2*es13^2 + s23^2*es23^2
    // then error on sigma1: e_sigma1 = sqrt(var(S1)) / (2*sigma1)

    double varS1 = s12*s12*es12*es12 + s13*s13*es13*es13 + s23*s23*es23*es23;
    double varS2 = s12*s12*es12*es12 + s23*s23*es23*es23 + s13*s13*es13*es13; // same form, ordering irrelevant
    double varS3 = s13*s13*es13*es13 + s23*s23*es23*es23 + s12*s12*es12*es12;

    double e_sigma_1 = 0.0;
    double e_sigma_2 = 0.0;
    double e_sigma_3 = 0.0;

    if (sigma_1 > 0) 
        e_sigma_1 = TMath::Sqrt(varS1) / (2.0 * sigma_1);
    else 
        e_sigma_1 = 0.0;

    if (sigma_2 > 0) 
        e_sigma_2 = TMath::Sqrt(varS2) / (2.0 * sigma_2);
    else 
        e_sigma_2 = 0.0;

    if (sigma_3 > 0) 
        e_sigma_3 = TMath::Sqrt(varS3) / (2.0 * sigma_3);
    else 
        e_sigma_3 = 0.0;

    cout.setf(std::ios::fixed); cout.precision(4);
    cout << "Results of time resolution with CFD = 20% (obtained from gaussian fits):\n"
         << "sigma12 (fit) = " << s12 << " +/- " << es12 << " ps\n"
         << "sigma23 (fit) = " << s23 << " +/- " << es23 << " ps\n"
         << "sigma13 (fit) = " << s13 << " +/- " << es13 << " ps\n\n"
         << "Derived per-detector resolutions (sigma):\n"
         << "sigma1 = " << sigma_1 << " +/- " << e_sigma_1 << " ps\n"
         << "sigma2 = " << sigma_2 << " +/- " << e_sigma_2 << " ps\n"
         << "sigma3 = " << sigma_3 << " +/- " << e_sigma_3 << " ps\n";

    // -----------------------
    //       Plot results
    // -----------------------
    TCanvas *c = new TCanvas("c", "Amplitude Comparison", 800, 600);
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
    c->Write();
    c1->Write();
    htriple1->Write();
    htriple2->Write();
    htriple3->Write();
    htime_12->Write();
    htime_23->Write();
    htime_13->Write();
    htime20_12->Write();
    htime20_23->Write();
    htime20_13->Write();
    htime30_12->Write();
    htime30_23->Write();
    htime30_13->Write();

    fout->Close(); // Close output file
    file->Close(); // Close input file
    cout << "Analysis complete. Results saved to muon_run22_analysis.root" << endl;

    return 0;
}