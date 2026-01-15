/*
    This code is to apply corrections to time resolution.
    To compile: c++ time_resolution.cpp `root-config --cflags --libs` -o timeRes.out
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
#include <TH2F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include "TProfile.h"
#include "TGraph2D.h"
#include <vector>
using namespace std;

double media(const vector<double>& v) 
{
    if (v.empty()) return 0; // evitare divisione per zero
    double somma = 0;
    for (double x : v) {
        somma += x;
    }
    return somma / v.size();
}

int main()
{
    gStyle->SetOptStat(0);

    // -------------------------------------------------
    //     Load run file with processed waveforms
    // -------------------------------------------------

    //TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root"; // Filename while running on lxplus
    TString filename = "/Users/anna/Developing/PhD/Testbeam2025/sampic_run22_final.root"; // Filename while running on Anna's machine
    //TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run19_final.root"; // Filename while running on Riccardo's machine
    cout << "Opening file: " << filename << endl;

    // -----------------------
    //    Set up branches
    // -----------------------
    const int MAXPULSES = 200;
    Int_t npulses;

    Double_t pulses_amplitude[MAXPULSES], pulses_integral[MAXPULSES], pulses_time_cfd10[MAXPULSES], pulses_time_cfd20[MAXPULSES], pulses_time_cfd30[MAXPULSES], integral[MAXPULSES], pulses_time_cfd50[MAXPULSES], pulses_time_cfd60[MAXPULSES], Baseline[MAXPULSES], pulses_peak_time[MAXPULSES],  pulses_rise_time[MAXPULSES], mean_1[MAXPULSES], mean_2[MAXPULSES], conteggi_canali1[64], conteggi_canali2[64] ;
    Bool_t  pulses_bad_pulse[MAXPULSES];
    Int_t HitFeb[3], Board[MAXPULSES], Channel[MAXPULSES];
    int tot_hit_feb;
    Double_t time_12, time_23, time_13, time20_12, time20_13, time20_23;
    Double_t time30_12, time30_13, time30_23, time50_12, peso, peso2;
    vector<vector<double>> tabella1(64);
    vector<vector<double>> tabella2(64);

    tree->SetBranchAddress("npulses", &npulses);
    tree->SetBranchAddress("pulses_amplitude", &pulses_amplitude);
    tree->SetBranchAddress("pulses_integral", &pulses_integral);
    tree->SetBranchAddress("pulses_time_cfd10", &pulses_time_cfd10);
    tree->SetBranchAddress("pulses_time_cfd20", &pulses_time_cfd20);
    tree->SetBranchAddress("pulses_time_cfd30", &pulses_time_cfd30);
    tree->SetBranchAddress("pulses_bad_pulse", &pulses_bad_pulse);
    tree->SetBranchAddress("HitFeb", &HitFeb);
    tree->SetBranchAddress("Board", &Board);
    tree->SetBranchAddress("pulses_time_cfd50", &pulses_time_cfd50);
    tree->SetBranchAddress("pulses_time_cfd60", &pulses_time_cfd60);
    tree->SetBranchAddress("pulses_baseline", &Baseline);
    tree->SetBranchAddress("Channel", &Channel);
    tree->SetBranchAddress("pulses_peak_time", &pulses_peak_time);
    tree->SetBranchAddress("pulses_rise_time", &pulses_rise_time);

    // ------------------------------
    //     Initialize histograms
    // ------------------------------
    TH1F *hAmpAll    = new TH1F("hAmpAll", "All pulse amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hAmp1Hit   = new TH1F("hAmp1Hit", "Single-hit event amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hQ         = new TH1F("hQ", "Pulse integral (namely charge);Integral [a.u.];Counts", 1000, 0, 10);
    TH1F *hTime20    = new TH1F("hTime20", "Constant fraction discriminator set at 20% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime30    = new TH1F("hTime30", "Constant fraction discriminator set at 30% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime50    = new TH1F("hTime50", "Constant fraction discriminator set at 50% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime60    = new TH1F("hTime60", "Constant fraction discriminator set at 60% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime60_2  = new TH1F("hTime60_2", "Constant fraction discriminator set at 60% Time;Time [ps];Counts", 1000, 0, 10000);

    // Histograms for triple hits
    TH1F *htriple1   = new TH1F("htriple1", "amplitude of a triple", 100, 0, 0.5);
    TH1F *htriple2   = new TH1F("htriple2", "amplitude of a triple", 100, 0, 0.5);
    TH1F *htriple3   = new TH1F("htriple3", "amplitude of a triple", 100, 0, 0.5);

    // Histograms for time differences for the different values of CDF
    TH1F *htime_12   = new TH1F("htime_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime_13   = new TH1F("htime_13", "time diff ;Time [ps];Counts", 100, -0.2, 0.2);
    TH1F *htime_23   = new TH1F("htime_23", "time diff ;Time [ps];Counts", 100, -0.2, 0.2);

    TH1F *htime20_12 = new TH1F("htime20_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);//unico gragfico che sto considereando
    TH1F *htime20_13 = new TH1F("htime20_13", "time diff ;Time [ps];Counts", 100, -100, 100);
    TH1F *htime20_23 = new TH1F("htime20_23", "time diff ;Time [ps];Counts", 100, -0.2, 0.2);

    TH1F *htime30_12 = new TH1F("htime30_12", "time diff ;Time [ps];Counts", 500, -2000, 2000);
    TH1F *htime30_13 = new TH1F("htime30_13", "time diff ;Time [ps];Counts", 500, -2000, 2000);
    TH1F *htime30_23 = new TH1F("htime30_23", "time diff ;Time [ps];Counts", 500, -2000, 2000);

    // Baseline, peak time, rise time, charge histograms
    TH1F *hcharge1   = new TH1F("hcharge1", "charge1", 1000, 0, 10);
    TH1F *hcharge2   = new TH1F("hcharge2", "charge2", 1000, 0, 10);
    TH1F *hrisetime1 = new TH1F("hrisetime1", "risetime1", 1000, -2000, 2000);
    TH1F *hrisetime2 = new TH1F("hrisetime2", "risetime2", 1000, -2000, 2000);
    TH1F *hbaseline1 = new TH1F("hbaseline", "baseline", 1000, 0.9, 1.1);
    TH1F *hbaseline2 = new TH1F("hbaseline2", "baseline2", 1000, 0.9, 1.1);
    TH1F *hpeaktime1 = new TH1F("hpeaktime1", "peaktime1", 1000, -20000, 20000);
    TH1F *hpeaktime2 = new TH1F("hpeaktime2", "peaktime2", 1000, -20000, 20000);

    TH1F *conteggixcanale1 = new TH1F("conteggixcanale1", "conteggixcanale1", 64, 0, 64);
    TH1F *conteggixcanale2 = new TH1F("conteggixcanale2", "conteggixcanale2", 64, 0, 64);
    TH1F *conteggixcanale3 = new TH1F("conteggixcanale3", "conteggixcanale3", 8, 0, 7);

    ///Amplitude vs integral histogram
    TGraph *hintegral1 = new TGraph();
    TGraph *hintegral2 = new TGraph();
    TGraph *hintegral3 = new TGraph();
    TGraph2D *timevsamplitude =new TGraph2D();
    TGraph *timevsintegral =new TGraph();
    TGraph2D *timevsamplitude2 =new TGraph2D();
    TGraph *timevsintegral2 =new TGraph();
    TGraph *timevsamplitude3 = new TGraph();

    TProfile *prof = new TProfile("prof", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);
    TProfile *prof2 = new TProfile("prof2", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);
    TProfile *prof3 = new TProfile("prof3", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);

    TGraph *cfdvschannel1 = new TGraph();

    // -----------------------
    //        Read tree
    // -----------------------
    Long64_t nentries = tree->GetEntries();
    bool cut[nentries]={false};
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
            //hTime->Fill(pulses_time_cfd20[idx]);
        }

        // ---- Triple hits for all three FEBs ----
        if (HitFeb[0] == 1 && HitFeb[1] == 1 && nGood == 3 && HitFeb[2]==1)
        {
            // Arrays indexed by FEB number (0,1,2)
            double ampFEB[3]  = {-1.0, -1.0, -1.0};
            double cfd10[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd20[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd30[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd50[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd60[3]   = {-9999.0, -9999.0, -9999.0};
            double base[3]    = {-9999.0, -9999.0, -9999.0};
            double peak_time[3]  = {-9999.0, -9999.0, -9999.0};
            double risetime[3] = {-9999.0, -9999.0, -9999.0};

            // Fill according to Board[]
            for (int j = 0; j < npulses; j++)
            {
                if (pulses_bad_pulse[j])
                    continue;
                int det = Board[j]; // Which detector this pulse belongs to (either 0, 1, 2)

                if (det < 0 || det > 2)
                {
                    cout << "WARNING: unexpected Board index = " << det << " (event " << i << ", pulse " << j << ")" << endl;
                    continue;
                }

                ampFEB[det]    = pulses_amplitude[j];
                cfd10[det]     = pulses_time_cfd10[j];
                cfd20[det]     = pulses_time_cfd20[j];
                cfd30[det]     = pulses_time_cfd30[j];
                integral[det]  = pulses_integral[j];
                cfd50[det]     = pulses_time_cfd50[j];
                cfd60[det]     = pulses_time_cfd60[j];
                base[det]      = Baseline[j];
                peak_time[det] = pulses_peak_time[j];
                risetime[det]  = pulses_rise_time[j];
                hTime20    -> Fill(cfd20[0]);
                hTime30    -> Fill(cfd30[0]);
                hTime50    -> Fill(cfd50[0]);
                hTime60    -> Fill(cfd60[0]);
                hTime60_2  -> Fill(cfd60[1]);
                hcharge1   -> Fill(integral[0]);
                hcharge2   -> Fill(integral[1]);
                hrisetime1 -> Fill(risetime[0]);
                hrisetime2 -> Fill(risetime[1]);
            }

            // Skip se un detector non ha dato segnale buono
            if (ampFEB[0] <= 0 || ampFEB[1] <= 0)
                continue;
            
            if(Channel[0]==12 && Channel[1]==12 && Channel[2]==5)
            {
                if(peak_time[1]<0)
                    continue;

                // Fill amplitude histograms; x 1000 scaling to get values in mV
                htriple1->Fill(ampFEB[0]);
                htriple2->Fill(ampFEB[1]);
                htriple3->Fill(ampFEB[2]);
                
                peso = Channel[0] / 8;
                if (peso > 7) 
                    peso = 7;

                peso2 = Channel[1] / 8;
                if (peso2 > 7)
                    peso2 = 7;

                // Average for every channel
                for(int i=0; i<64; i++)
                {
                    if(Channel[0]==i)
                        tabella1[i].push_back(peak_time[0]);
                }
                for(int i=0; i<64; i++)
                {
                    if(Channel[1]==i)
                        tabella2[i].push_back(peak_time[1]);
                }

                double cfd30_co=cfd30[0]-(3438-27.53/ampFEB[0]);
                double cfd30_co1=cfd30[1]-(3819-30.49/ampFEB[1]);
                double cfd30_co2=cfd30[2]-(3902-24.2/ampFEB[2]);

                // ---- CFD = 10% ----
                //time_12 = cfd60[0] - cfd60[1];
                //time_23 = cfd10[1] - cfd10[2];
                //time_13 = cfd10[0] - cfd10[2];

                //htime_12->Fill(time_12);
                //htime_23->Fill(time_23);
                //htime_13->Fill(time_13);

                // ---- CFD = 20% ----
                //time20_12 = cfd20[0] - cfd20[1];
                //if(fabs(time20_12)>)continue;
                //time20_23 = cfd20[1] - cfd20[2];
                //time20_13 = cfd20[0] - cfd20[2];

                //htime20_12->Fill(time20_12);
                //htime20_23->Fill(time20_23);
                //htime20_13->Fill(time20_13);
                hpeaktime1->Fill(peak_time[0]);
                hpeaktime2->Fill(peak_time[1]);

                // ---- CFD = 60% ----
                //double timedif1= peak_time[0]-cfd60[0];
                //double timedif2= peak_time[1]-cfd60[1];
                time30_12 = cfd30_co-cfd30_co1;
                time30_23 = cfd30_co1 - cfd30_co2;
                time30_13 = cfd30_co - cfd30_co2;
                htime30_12->Fill(time30_12);
                htime30_23->Fill(time30_23);
                htime30_13->Fill(time30_13);

                hbaseline1->Fill(base[0]);
                hbaseline2->Fill(base[1]);

                ///Integral vs amplitude distribution
                timevsamplitude->SetPoint(timevsamplitude->GetN(), ampFEB[0], cfd30[0], peso);
                timevsamplitude2->SetPoint(timevsamplitude2->GetN(), ampFEB[1], cfd30[1], peso2);
                timevsamplitude3->SetPoint(timevsamplitude3->GetN(), ampFEB[2], cfd30[2]);

                prof->Fill(ampFEB[0], cfd60[0]);
                prof2->Fill(ampFEB[1], cfd60[1]);
                prof3->Fill(ampFEB[2], cfd60[2]);

                timevsintegral->SetPoint(timevsintegral->GetN(), integral[0], cfd60[0]);
                hintegral1->SetPoint(hintegral1->GetN(),ampFEB[0],integral[0]);
                hintegral2->SetPoint(hintegral2->GetN(),ampFEB[1],integral[1]);
                hintegral3->SetPoint(hintegral3->GetN(),ampFEB[2],integral[2]);

                conteggixcanale1->Fill(Channel[0]);
                conteggixcanale2->Fill(Channel[1]);
                conteggixcanale3->Fill(Channel[2]);

                if(ampFEB[0]<0)
                    cout<< "Amplitude < 0; event number: " << i << endl;
                cfdvschannel1->SetPoint(cfdvschannel1->GetN(), Channel[0], cfd60[0]);
            }
        }
        // Progress indicator for sanity :)
        if (i % 100000 == 0) cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }
    cout << endl << "Finished processing " << nentries << " events." << endl;

    for(int i=0; i<64; i++)
    {
        mean_1[i]=media(tabella1[i]);
        //cout << "media canali primo riv: " << mean_1[i] << endl;
    }



    for(int i=0; i<64; i++){
        mean_2[i]=media(tabella2[i]);
        //cout << "media canali secondo riv: " << mean_2[i] << endl;
    }

    // New tree (is it really needed here though?)
    //TTree *copy = tree->CloneTree();
    TTree *tout = new TTree("channelMeans", "Medie per canale");
    tout->Branch("mean_1", mean_1, "mean_1[64]/D");
    tout->Branch("mean_2", mean_2, "mean_2[64]/D");
    tout->Fill();

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
    TF1 *fpolyacharge1 = (TF1*)fpolyaAmpl1->Clone("fpolyacharge1");
    TF1 *fpolyacharge2 = (TF1*)fpolyaAmpl1->Clone("fpolyacharge2");

    // Canvas 1
    TCanvas *c1 = new TCanvas("c1", "htriple1", 900, 700);
    htriple1->SetLineColor(kAzure+1);
    htriple1->SetLineWidth(3);
    htriple1->SetFillColorAlpha(kAzure-4, 0.35);
    htriple1->SetTitle("Amplitudes of waveforms on the first FEB");
    htriple1->Draw("HIST");
    htriple1->Fit(fpolyaAmpl1,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
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
    //fpolyaAmpl3->Draw("SAME");
    c3->Update();
    c3->SaveAs("hist3.pdf");

    // ------------ Canvas charge1, charge2, risetime1, risetime2, peaktime1, peaktime2 ------------
    TCanvas *ccharge1 = new TCanvas("ccharge1", "charge1", 900, 700);
    hcharge1->SetLineColor(kAzure+1);
    hcharge1->SetLineWidth(3);
    hcharge1->SetFillColorAlpha(kAzure-4, 0.35);
    hcharge1->SetTitle("Charge of waveforms on the first FEB");
    hcharge1->Draw("HIST");
    hcharge1->Fit(fpolyacharge1,"R");
    ccharge1->Update();

    TCanvas *ccharge2 = new TCanvas("ccharge2", "charge2", 900, 700);
    hcharge2->SetLineColor(kAzure+1);
    hcharge2->SetLineWidth(3);
    hcharge2->SetFillColorAlpha(kAzure-4, 0.35);
    hcharge2->SetTitle("Charge of waveforms on the second FEB");
    hcharge2->Draw("HIST");
    hcharge2->Fit(fpolyacharge2,"R");
    ccharge2->Update();

    TCanvas *crisetime1 = new TCanvas("crisetime1", "risetime1", 900, 700);
    hrisetime1->SetLineColor(kAzure+1);
    hrisetime1->SetLineWidth(3);
    hrisetime1->SetFillColorAlpha(kAzure-4, 0.35);
    hrisetime1->SetTitle("risetime of waveforms on the first FEB");
    hrisetime1->Draw("HIST");
    crisetime1->Update();

    TCanvas *crisetime2 = new TCanvas("crisetime2", "risetime2", 900, 700);
    hrisetime2->SetLineColor(kAzure+1);
    hrisetime2->SetLineWidth(3);
    hrisetime2->SetFillColorAlpha(kAzure-4, 0.35);
    hrisetime2->SetTitle("risetime of waveforms on the second FEB");
    hrisetime2->Draw("HIST");
    crisetime2->Update();


    TCanvas *cpeaketime1 = new TCanvas("cpeaktime1", "peaktime1", 900, 700);
    hpeaktime1->SetLineColor(kAzure+1);
    hpeaktime1->SetLineWidth(3);
    hpeaktime1->SetFillColorAlpha(kAzure-4, 0.35);
    hpeaktime1->SetTitle("peaktime of waveforms on the first FEB");
    hpeaktime1->Draw("HIST");
    cpeaketime1->Update();

    TCanvas *cpeaketime2 = new TCanvas("cpeaktime2", "peaktime2", 900, 700);
    hpeaktime2->SetLineColor(kAzure+1);
    hpeaktime2->SetLineWidth(3);
    hpeaktime2->SetFillColorAlpha(kAzure-4, 0.35);
    hpeaktime2->SetTitle("peaktime of waveforms on the second FEB");
    hpeaktime2->Draw("HIST");
    cpeaketime2->Update();

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

    // Amplitude vs integral plot
    TCanvas *cintegral = new TCanvas("cintegral", "amplitudevsintegral", 900, 700);

    hintegral1->SetTitle("Amplitude vs Integral;amplitude(V);integral");
    hintegral1->SetName("hintegral1");
    hintegral1->SetMarkerStyle(20);
    hintegral1->SetMarkerSize(1.0);
    hintegral1->SetMarkerColor(kAzure+2);
    hintegral1->SetLineColor(0);
    hintegral1->SetLineWidth(0);

    double xmin = hintegral1->GetXaxis()->GetXmin();
    double xmax = hintegral1->GetXaxis()->GetXmax();

    TF1* avsi = new TF1("avsi", "[0]+[1]*x", 0, 2);
    //avsi->SetParameter(0, 0.0);    // intercetta iniziale
    avsi->SetParLimits(1, 12, 26);   // pendenza iniziale stimata
    avsi->SetNpx(5000);
    hintegral1->Fit(avsi);
    hintegral1->Draw("AP");
    avsi->SetLineColor(kRed);
    avsi->SetLineWidth(2);
    avsi->Draw("same");
    cintegral->Update();

    TCanvas *cintegral2 = new TCanvas("cintegral2", "amplitudevsintegral", 900, 700);
    hintegral2->SetName("hintegral2");
    hintegral2->SetTitle("Amplitude vs Integral;amplitude;integral");
    hintegral2->SetMarkerStyle(20);
    hintegral2->SetMarkerSize(1.0);
    hintegral2->SetMarkerColor(kAzure+2);
    hintegral2->SetLineColor(0); 
    hintegral2->SetLineWidth(0);

    TF1* avsi2 = new TF1("avsi2", "[0] + [1]*x");
    avsi2->SetNpx(500);
    hintegral2->Fit(avsi2, "R");
    avsi2->SetLineColor(kRed);
    avsi2->SetLineWidth(2);
    hintegral2->Draw("AP");
    hintegral2->Fit(avsi, "R");  // "R" = fit solo nell'intervallo di x specificato
    hintegral2->Draw("AP");
    avsi2->Draw("same");          // adesso la linea del fit viene disegnata*/
    cintegral2->Update();

    // Counts per channel plots
    TCanvas *cconteggi = new TCanvas("cconteggi", "cconteggi", 900, 700);
    conteggixcanale1->SetLineColor(kAzure+1);
    conteggixcanale1->SetLineWidth(3);
    conteggixcanale1->SetFillColorAlpha(kAzure-4, 0.35);
    conteggixcanale1->SetTitle("eventi per canale");
    cconteggi->Update();

    TCanvas *cconteggi2 = new TCanvas("cconteggi2", "cconteggi2", 900, 700);
    conteggixcanale2->SetLineColor(kAzure+1);
    conteggixcanale2->SetLineWidth(3);
    conteggixcanale2->SetFillColorAlpha(kAzure-4, 0.35);
    conteggixcanale2->SetTitle("eventi per canale");
    cconteggi2->Update();

    TCanvas *cconteggi3 = new TCanvas("cconteggi3", "cconteggi3", 900, 700);
    conteggixcanale3->SetLineColor(kAzure+1);
    conteggixcanale3->SetLineWidth(3);
    conteggixcanale3->SetFillColorAlpha(kAzure-4, 0.35);
    conteggixcanale3->SetTitle("eventi per canale");
    cconteggi3->Update();

    // Amplitude vs time plot
    TCanvas *cAmpVsTime = new TCanvas("cAmpVsTime", "Amplitude vs Time", 900, 700);
    cAmpVsTime->SetGrid();
    cAmpVsTime->SetTickx();
    cAmpVsTime->SetTicky();

    timevsamplitude3->SetTitle("Amplitude vs Time CFD 30%;Amplitude [a.u.];Time [ps]");
    timevsamplitude3->SetMarkerStyle(20);
    timevsamplitude3->SetMarkerSize(1.3);
    timevsamplitude3->SetMarkerColor(kAzure+2);
    timevsamplitude3->SetLineColor(kAzure+1);
    timevsamplitude3->SetLineWidth(2);

    timevsamplitude3->GetXaxis()->SetTitleSize(0.045);
    timevsamplitude3->GetYaxis()->SetTitleSize(0.045);
    timevsamplitude3->GetXaxis()->SetLabelSize(0.04);
    timevsamplitude3->GetYaxis()->SetLabelSize(0.04);

    timevsamplitude3->Draw("AP");
    cAmpVsTime->Update();

    // ----------------------------------------------------------------------------------------------------

    cout << "Time walk and Slewing Correction Plots" << endl;
    TCanvas *ctimeamlitude = new TCanvas("ctimeaplitude", "amplitudevsintegral", 900, 700);
    //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);
    timevsamplitude->SetName("timeaplitude");
    timevsamplitude->SetTitle("SAT vs amplitude; Amplitude; SAT[us]");

    gStyle->SetPalette(kViridis);
    //timevsamplitude->Fit("timewalk", "R");
    timevsamplitude->Draw("COLZ");
    ctimeamlitude->Update();

    TCanvas *ctimeamlitude2 = new TCanvas("ctimeaplitude2", "amplitudevsintegral", 900, 700);
    //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);
    timevsamplitude2->SetName("timeaplitude2");
    timevsamplitude2->SetTitle("SAT vs amplitude; Amplitude; SAT[us]");

    gStyle->SetPalette(kViridis);
    //timevsamplitude->Fit("timewalk", "R");
    timevsamplitude2->Draw("COLZ");
    ctimeamlitude2->Update();

    TCanvas *ctvsint = new TCanvas("ctvsint", "Time vs Integral", 900, 700);
    ctvsint ->SetGrid();

    timevsintegral->SetName("timevsintegral");
    timevsintegral->SetLineColor(0);
    timevsintegral->SetLineWidth(0);
    timevsintegral->SetMarkerStyle(20);
    timevsintegral->SetMarkerSize(1.1);
    timevsintegral->SetMarkerColor(kBlue+1);

    timevsintegral->SetTitle("SAT vs integral;Integral;SAT [ps]");
    timevsintegral->GetXaxis()->SetTitleSize(0.045);
    timevsintegral->GetYaxis()->SetTitleSize(0.045);
    timevsintegral->GetXaxis()->SetLabelSize(0.04);
    timevsintegral->GetYaxis()->SetLabelSize(0.04);
    timevsintegral->Draw("AP");
    ctvsint->Update();

    TCanvas *cmediabin = new TCanvas("cmediabin", "Slewing Correction", 800, 600);
    //prof->Fit("pol3");
    TF1 *empirico = new TF1("empirico", "[0] + [1]/x", 0, 1);
    prof->Fit(empirico, "R");
    prof->Draw();
    cmediabin->Update();

    TCanvas *cmediabin2 = new TCanvas("cmediabin2", "Slewing Correction", 800, 600);
    //prof->Fit("pol3");
    TF1 *empirico2 = new TF1("empirico2", "[0] + [1]/x", 0, 1);
    prof2->Fit(empirico2, "R");
    prof2->Draw();
    cmediabin2->Update();

    TCanvas *cmediabin3 = new TCanvas("cmediabin2", "Slewing Correction", 800, 600);
    //prof->Fit("pol3");
    TF1 *empirico3 = new TF1("empirico3", "[0] + [1]/x", 0, 1);
    prof3->Fit(empirico3, "R");
    prof3->Draw();
    cmediabin3->Update();

    /*TCanvas *cintegral3 = new TCanvas("cintegral3", "amplitudevsintegral", 900, 700);
    hintegral3->SetName("hintegral1");
    hintegral3->SetTitle("Amplitude vs Integral;amplitude;integral");

    hintegral3->SetMarkerStyle(20);
    hintegral3->SetMarkerSize(1.0);
    hintegral3->SetMarkerColor(kAzure+2);
    hintegral3->SetLineColor(0);    // nessuna linea
    hintegral3->SetLineWidth(0);    // nessuna linea
    hintegral3->Draw("AP");
    cintegral3->Update();*/

    // Baseline
    TCanvas *cbaseline = new TCanvas("cbaseline", "baseline", 900, 700);
    TF1 *f_gaus_b = new TF1("f_gaus", "gaus", 0.003, 0.009);
    hbaseline1->Fit(f_gaus_b, "R"); // "R" usa il range definito in TF1

    hbaseline1->SetLineColor(kAzure+1);
    hbaseline1->SetLineWidth(3);
    hbaseline1->SetFillColorAlpha(kAzure-4, 0.35);
    hbaseline1->SetTitle("baseline fist det");
    hbaseline1->Draw("HIST");
    cbaseline->Update();

    double amp_baseline   = f_gaus_b->GetParameter(0);
    double mean_baseline  = f_gaus_b->GetParameter(1);
    double sigma_baseline = f_gaus_b->GetParameter(2);

    //baseline
    TCanvas *cbaseline2 = new TCanvas("cbaseline2", "baseline2", 900, 700);
    TF1 *f_gaus_b2 = new TF1("f_gaus", "gaus2", 0.003, 0.009);
    hbaseline1->Fit(f_gaus_b2, "R"); // "R" usa il range definito in TF1

    hbaseline2->SetLineColor(kAzure+1);
    hbaseline2->SetLineWidth(3);
    hbaseline2->SetFillColorAlpha(kAzure-4, 0.35);
    hbaseline2->SetTitle("baseline second det");
    hbaseline2->Draw("HIST");
    cbaseline2->Update();

    // 5. Estrai i risultati dal TF1 dopo il fit (usa l'oggetto f_gaus)
    double amp_baseline2   = f_gaus_b2->GetParameter(0);
    double mean_baseline2 = f_gaus_b2->GetParameter(1); // Questo sarà il risultato del fit
    double sigma_baseline2 = f_gaus_b2->GetParameter(2);


    // ... estrazione degli errori
    /* TCanvas *ctime23 = new TCanvas("ctime23", "htriple1", 900, 700);
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
    double e_sigma_13 = f3->GetParError(2);*/

    TCanvas *ctimeperchannel = new TCanvas("ctimeperchannel", "Time vs Integral", 900, 700);
    ctvsint ->SetGrid();

    cfdvschannel1->SetName("cfdvschannel1");
    cfdvschannel1->SetLineColor(0);
    cfdvschannel1->SetLineWidth(0);
    cfdvschannel1->SetMarkerStyle(20);
    cfdvschannel1->SetMarkerSize(1.1);
    cfdvschannel1->SetMarkerColor(kBlue+1);
    cfdvschannel1->SetTitle("channel vs SAT;Channel;SAT [ps]");
    cfdvschannel1->GetXaxis()->SetTitleSize(0.045);
    cfdvschannel1->GetYaxis()->SetTitleSize(0.045);
    cfdvschannel1->GetXaxis()->SetLabelSize(0.04);
    cfdvschannel1->GetYaxis()->SetLabelSize(0.04);
    cfdvschannel1->Draw("AP");

    ctimeperchannel->Update();

    // ------------ Gaussian plot of the time difference ---------------

    TCanvas *ctime12 = new TCanvas("ctime12", "htriple1", 900, 700);

    TF1 *f_gaus = new TF1("f_gaus", "gaus", -1000, 1000);
    htime30_12->Fit(f_gaus, "R"); // "R" usa il range definito in TF1
    htime30_12->SetLineColor(kAzure+1);
    htime30_12->SetLineWidth(3);
    htime30_12->SetFillColorAlpha(kAzure-4, 0.35);
    htime30_12->SetTitle("time difference 1&2 detector");
    htime30_12->Draw("HIST");
    ctime12->Update();

    double amp_12   = f_gaus->GetParameter(0);
    double mean_12  = f_gaus->GetParameter(1);
    double sigma_12 = f_gaus->GetParameter(2);
    double e_amp_12   = f_gaus->GetParError(0);
    double e_mean_12  = f_gaus->GetParError(1);
    double e_sigma_12 = f_gaus->GetParError(2);

    TCanvas *ctime23 = new TCanvas("ctime23", "htriple1", 900, 700);
    htime30_23->SetLineColor(kAzure+1);
    htime30_23->SetLineWidth(3);
    htime30_23->SetFillColorAlpha(kAzure-4, 0.35);
    htime30_23->SetTitle("time difference 2&3 detector");
    htime30_23->Draw("HIST");
    htime30_23->Fit("gaus"); // Fit su htime20_23
    ctime23->Update();

    TF1 *f2 = htime30_23->GetFunction("gaus");
    double amp_23   = f2->GetParameter(0);  ///taking gaus parameters
    double mean_23  = f2->GetParameter(1);
    double sigma_23 = f2->GetParameter(2);

    double e_amp_23   = f2->GetParError(0);
    double e_mean_23  = f2->GetParError(1);
    double e_sigma_23 = f2->GetParError(2);

    TCanvas *ctime13 = new TCanvas("ctime13", "htriple1", 900, 700);
    htime30_13->SetLineColor(kAzure+1);
    htime30_13->SetLineWidth(3);
    htime30_13->SetFillColorAlpha(kAzure-4, 0.35);
    htime30_13->SetTitle("time difference 1&3 detector");
    htime30_13->Draw("HIST");
    htime30_13->Fit("gaus"); // Fit su htime20_13
    ctime13->Update();

    // Take gaussian from the correct histogram htime20_13
    TF1 *f3 = htime30_13->GetFunction("gaus");
    double amp_13   = f3->GetParameter(0);  // taking gaus parameters
    double mean_13  = f3->GetParameter(1);
    double sigma_13 = f3->GetParameter(2);

    double e_amp_13   = f3->GetParError(0);
    double e_mean_13  = f3->GetParError(1);
    double e_sigma_13 = f3->GetParError(2);

    // --------------- Time resolution (with error propagation) ---------------

    // Use: sigma_12^2 = sigma1^2 + sigma2^2, etc.
    // Inversion:
    //              sigma1^2 = 1/2 * ( sigma_12^2 + sigma_13^2 - sigma_23^2 )
    //              sigma2^2 = 1/2 * ( sigma_12^2 + sigma_23^2 - sigma_13^2 )
    //              sigma3^2 = 1/2 * ( sigma_13^2 + sigma_23^2 - sigma_12^2 )

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
    //     Save everything
    // -----------------------
    TFile *fout = new TFile("muon_run_only2_TProfile_run19_timeresolution_channel12.root", "RECREATE");
    hAmpAll->Write();
    hAmp1Hit->Write();
    hQ->Write();
    hTime20->Write();
    hTime30->Write();
    hTime50->Write();
    hTime60->Write();
    hTime60_2->Write();
    //c->Write();
    c1->Write();
    htriple1->Write();
    htriple2->Write();
    htriple3->Write();
    htime_12->Write();
    //htime_23->Write();
    //htime_13->Write();
    htime20_12->Write();
    //htime20_23->Write();
    //htime20_13->Write();
    htime30_12->Write();
    htime30_23->Write();
    htime30_13->Write();
    hintegral1->Write();
    hintegral2->Write();
    hintegral3->Write();
    hcharge1->Write();
    hcharge2->Write();
    hrisetime1->Write();
    hrisetime2->Write();
    hpeaktime1->Write();
    hpeaktime2->Write();
    timevsamplitude->Write();
    timevsamplitude2->Write();
    timevsamplitude3->Write();
    hbaseline1->Write();
    hbaseline2->Write();
    timevsintegral->Write();
    prof->Write();
    prof2->Write();
    prof3->Write();
    cmediabin->Write();
    cfdvschannel1->Write();
    cintegral->Write();
    cintegral2->Write();
    tout->Write();
    conteggixcanale1->Write();
    conteggixcanale2->Write();
    conteggixcanale3->Write();
    cAmpVsTime->Write();

    fout->Close(); // Close output file
    file->Close(); // Close input file
    cout << "Analysis complete." << endl;

    return 0;
}