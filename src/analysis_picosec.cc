#include "analysis_picosec.h"
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>

using namespace std;

double CalculateAverage(const vector<double>& v) {
    if (v.empty()) return 0;
    double sum = 0;
    for (double x : v) {
        sum += x;
    }
    return sum / v.size();
}

TFile* OpenInputFile(const string &filename) 
{
    cout << "Opening file: " << filename << endl;
    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        cerr << "Error: cannot open " << filename << endl;
        return nullptr;
    }
    return file;
}

void InitializeHistograms(Histograms &h) 
{
    // Basic histograms
    h.hAmpAll = new TH1F("hAmpAll", "All pulse amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    h.hAmp1Hit = new TH1F("hAmp1Hit", "Single-Hit event amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    h.hQ = new TH1F("hQ", "Pulse integral (namely charge);Integral [a.u.];Counts", 1000, 0, 10);
    
    // CFD time histograms
    h.hTime20 = new TH1F("hTime20", "CFD 20% Time;Time [ps];Counts", 1000, 0, 10000);
    h.hTime30 = new TH1F("hTime30", "CFD 30% Time;Time [ps];Counts", 1000, 0, 10000);
    h.hTime50 = new TH1F("hTime50", "CFD 50% Time;Time [ps];Counts", 1000, 0, 10000);
    h.hTime60 = new TH1F("hTime60", "CFD 60% Time;Time [ps];Counts", 1000, 0, 10000);
    h.hTime60_2 = new TH1F("hTime60_2", "CFD 60% Time;Time [ps];Counts", 1000, 0, 10000);
    
    // Triple hit amplitudes
    h.htriple1 = new TH1F("htriple1", "Pulse amplitude for first detector", 1000, 0, 0.5);
    h.htriple2 = new TH1F("htriple2", "Pulse amplitude for second detector", 1000, 0, 0.5);
    h.htriple3 = new TH1F("htriple3", "Pulse amplitude for third detector", 1000, 0, 0.5);
    
    // Time differences
    h.htime_12 = new TH1F("htime_12", "Time difference;Time [ps];Counts", 100, -2000, 2000);
    h.htime_13 = new TH1F("htime_13", "Time difference;Time [ps];Counts", 100, -0.2, 0.2);
    h.htime_23 = new TH1F("htime_23", "Time difference;Time [ps];Counts", 100, -0.2, 0.2);
    
    h.htime20_12 = new TH1F("htime20_12", "Time difference;Time [ps];Counts", 100, -2000, 2000);
    h.htime20_13 = new TH1F("htime20_13", "Time difference;Time [ps];Counts", 100, -100, 100);
    h.htime20_23 = new TH1F("htime20_23", "Time difference;Time [ps];Counts", 100, -0.2, 0.2);
    
    h.htime30_12 = new TH1F("htime30_12", "Time difference;Time [ps];Counts", 100, -2000, 2000);
    h.htime30_13 = new TH1F("htime30_13", "Time difference;Time [ps];Counts", 100, -0.2, 0.2);
    h.htime30_23 = new TH1F("htime30_23", "Time difference;Time [ps];Counts", 100, -0.2, 0.2);
    
    // Charge, risetime, peaktime, baseline
    h.hcharge1 = new TH1F("hcharge1", "charge1", 1000, 0, 10);
    h.hcharge2 = new TH1F("hcharge2", "charge2", 1000, 0, 10);
    h.hrisetime1 = new TH1F("hrisetime1", "risetime1", 1000, -2000, 2000);
    h.hrisetime2 = new TH1F("hrisetime2", "risetime2", 1000, -2000, 2000);
    h.hbaseline1 = new TH1F("hbaseline", "baseline", 1000, 0.9, 1.1);
    h.hbaseline2 = new TH1F("hbaseline2", "baseline2", 1000, 0.9, 1.1);
    h.hpeaktime1 = new TH1F("hpeaktime1", "peaktime1", 1000, -20000, 20000);
    h.hpeaktime2 = new TH1F("hpeaktime2", "peaktime2", 1000, -20000, 20000);
    
    // Channel counts
    h.conteggixcanale1 = new TH1F("conteggixcanale1", "Counts per channel, first detector", 1000, 0, 64);
    h.conteggixcanale2 = new TH1F("conteggixcanale2", "Counts per channel, second detector", 1000, 0, 64);
    h.conteggixcanale3 = new TH1F("conteggixcanale3", "Counts per channel, third detector", 100, 0, 7);
    
    // Hit maps
    h.mapdet1 = new TH2F("hitmapdet1", "Hitmap detector 1;x;y", 8, 0, 8, 8, 0, 8);
    h.mapdet2 = new TH2F("hitmapdet2", "Hitmap detector 2;x;y", 8, 0, 8, 8, 0, 8);
    
    // Graphs
    h.hintegral1 = new TGraph();
    h.hintegral2 = new TGraph();
    h.hintegral3 = new TGraph();
    h.timevsamplitude = new TGraph2D();
    h.timevsamplitude2 = new TGraph2D();
    h.timevsintegral = new TGraph();
    h.timevsintegral2 = new TGraph();
    h.cfdvschannel1 = new TGraph();
    
    // Profiles
    h.prof = new TProfile("prof", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);
    h.prof2 = new TProfile("prof2", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);
}

void SetupTreeBranches(TTree *tree, TreeBranches &b) 
{
    tree->SetBranchAddress("npulses", &b.npulses);
    tree->SetBranchAddress("pulses_amplitude", &b.pulses_amplitude);
    tree->SetBranchAddress("pulses_integral", &b.pulses_integral);
    tree->SetBranchAddress("pulses_time_cfd10", &b.pulses_time_cfd10);
    tree->SetBranchAddress("pulses_time_cfd20", &b.pulses_time_cfd20);
    tree->SetBranchAddress("pulses_time_cfd30", &b.pulses_time_cfd30);
    tree->SetBranchAddress("pulses_time_cfd50", &b.pulses_time_cfd50);
    tree->SetBranchAddress("pulses_time_cfd60", &b.pulses_time_cfd60);
    tree->SetBranchAddress("pulses_bad_pulse", &b.pulses_bad_pulse);
    tree->SetBranchAddress("HitFeb", &b.HitFeb);
    tree->SetBranchAddress("Board", &b.Board);
    tree->SetBranchAddress("Channel", &b.Channel);
    tree->SetBranchAddress("pulses_baseline", &b.Baseline);
    tree->SetBranchAddress("pulses_peak_time", &b.pulses_peak_time);
    tree->SetBranchAddress("pulses_rise_time", &b.pulses_rise_time);
    tree->SetBranchAddress("pulses_channel_x", &b.pulses_channel_x);
    tree->SetBranchAddress("pulses_channel_y", &b.pulses_channel_y);
}

void ProcessEvents(TTree *tree, TreeBranches &b, Histograms &h,
                   vector<vector<double>> &tabella1, 
                   vector<vector<double>> &tabella2) 
{
    Long64_t nentries = tree->GetEntries();
    
    for (Long64_t i = 0; i < nentries; i++) 
    {
        tree->GetEntry(i);
        
        int nGood = 0;
        int idx = -1;
        
        // Fill basic histograms
        for (int j = 0; j < b.npulses; j++) 
        {
            nGood++;
            idx = j;
            h.hAmpAll->Fill(b.pulses_amplitude[j]);
        }
        
        if (nGood == 1 && idx >= 0) 
        {
            h.hAmp1Hit->Fill(b.pulses_amplitude[idx]);
            h.hQ->Fill(b.pulses_integral[idx]);
        }
        
        // Process double/triple hits
        if (b.HitFeb[0] >= 1 && b.HitFeb[1] >= 1) 
        {
            double ampFEB[3] = {-1.0, -1.0, -1.0};
            double cfd10[3] = {-9999.0, -9999.0, -9999.0};
            double cfd20[3] = {-9999.0, -9999.0, -9999.0};
            double cfd30[3] = {-9999.0, -9999.0, -9999.0};
            double cfd50[3] = {-9999.0, -9999.0, -9999.0};
            double cfd60[3] = {-9999.0, -9999.0, -9999.0};
            double base[3] = {-9999.0, -9999.0, -9999.0};
            double peak_time[3] = {-9999.0, -9999.0, -9999.0};
            double risetime[3] = {-9999.0, -9999.0, -9999.0};
            double integral[3] = {-9999.0, -9999.0, -9999.0};
            
            for (int j = 0; j < b.npulses; j++) 
            {
                int det = b.Board[j];
                if (det < 0 || det > 2) continue;
                
                ampFEB[det] = b.pulses_amplitude[j];
                cfd10[det] = b.pulses_time_cfd10[j];
                cfd20[det] = b.pulses_time_cfd20[j];
                cfd30[det] = b.pulses_time_cfd30[j];
                cfd50[det] = b.pulses_time_cfd50[j];
                cfd60[det] = b.pulses_time_cfd60[j];
                base[det] = b.Baseline[j];
                peak_time[det] = b.pulses_peak_time[j];
                risetime[det] = b.pulses_rise_time[j];
                integral[det] = b.pulses_integral[j];
                
                h.hTime20->Fill(cfd20[0]);
                h.hTime30->Fill(cfd30[0]);
                h.hTime50->Fill(cfd50[0]);
                h.hTime60->Fill(cfd60[0]);
                h.hTime60_2->Fill(cfd60[1]);
                h.hcharge1->Fill(integral[0]);
                h.hcharge2->Fill(integral[1]);
                h.hrisetime1->Fill(risetime[0]);
                h.hrisetime2->Fill(risetime[1]);
                
                if (det == 0)
                    h.mapdet1->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);
                if (det == 1)
                    h.mapdet2->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);
            }
            
            h.htriple1->Fill(ampFEB[0]);
            h.htriple2->Fill(ampFEB[1]);
            h.htriple3->Fill(ampFEB[2]);
            
            // Store peak times for channel averaging
            for (int k = 0; k < NCHANNELS; k++) 
            {
                if (b.Channel[0] == k)
                    tabella1[k].push_back(peak_time[0]);
                if (b.Channel[1] == k)
                    tabella2[k].push_back(peak_time[1]);
            }
            
            // Time differences
            double time30_12 = cfd60[0] - cfd60[1];
            h.htime30_12->Fill(time30_12);
            
            h.hpeaktime1->Fill(peak_time[0]);
            h.hpeaktime2->Fill(peak_time[1]);
            h.hbaseline1->Fill(base[0]);
            h.hbaseline2->Fill(base[1]);
            
            // Fill graphs
            double peso = 1.0;
            h.timevsamplitude->SetPoint(h.timevsamplitude->GetN(), ampFEB[0], cfd60[0], peso);
            h.timevsamplitude2->SetPoint(h.timevsamplitude2->GetN(), ampFEB[1], cfd60[1], peso);
            h.prof->Fill(ampFEB[0], cfd60[0]);
            h.prof2->Fill(ampFEB[1], cfd60[1]);
            h.timevsintegral->SetPoint(h.timevsintegral->GetN(), integral[0], cfd60[0]);
            h.hintegral1->SetPoint(h.hintegral1->GetN(), ampFEB[0], integral[0]);
            h.hintegral2->SetPoint(h.hintegral2->GetN(), ampFEB[1], integral[1]);
            h.conteggixcanale1->Fill(b.Channel[0]);
            h.conteggixcanale2->Fill(b.Channel[1]);
            h.conteggixcanale3->Fill(b.Channel[2]);
            h.cfdvschannel1->SetPoint(h.cfdvschannel1->GetN(), b.Channel[0], cfd60[0]);
        }
        
        if (i % 100000 == 0)
            cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }
    cout << endl;
}

void FitHistograms(Histograms &h) 
{
    // Polya function
    TF1 *fpolyaAmpl1 = new TF1("fpolyaAmpl1", 
        "([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])", 
        0, 800);
    fpolyaAmpl1->SetParameter(1, 180);
    fpolyaAmpl1->SetNpx(10000);
    
    h.htriple1->Fit(fpolyaAmpl1, "RQ");
    
    TF1 *fpolyaAmpl2 = (TF1*)fpolyaAmpl1->Clone("fpolyaAmpl2");
    h.htriple2->Fit(fpolyaAmpl2, "RQ");
    
    // Gaussian fit for time differences
    TF1 *f_gaus = new TF1("f_gaus", "gaus", -1000, 1000);
    h.htime30_12->Fit(f_gaus, "RQ");
    
    // Baseline fits
    TF1 *f_gaus_b = new TF1("f_gaus_b", "gaus", 0.9, 1.1);
    h.hbaseline1->Fit(f_gaus_b, "RQ");
    
    TF1 *f_gaus_b2 = new TF1("f_gaus_b2", "gaus", 0.9, 1.1);
    h.hbaseline2->Fit(f_gaus_b2, "RQ");
    
    // Slewing correction
    TF1 *empirico = new TF1("empirico", "[0] + [1]/x", 0, 1);
    h.prof->Fit(empirico, "RQ");
    
    TF1 *empirico2 = new TF1("empirico2", "[0] + [1]/x", 0, 1);
    h.prof2->Fit(empirico2, "RQ");
}

void CreateCanvases(Histograms &h) 
{
    gStyle->SetOptStat(0);
    
    // Amplitude canvases
    TCanvas *c1 = new TCanvas("c1", "htriple1", 900, 700);
    h.htriple1->SetLineColor(kAzure+1);
    h.htriple1->SetLineWidth(3);
    h.htriple1->SetFillColorAlpha(kAzure-4, 0.35);
    h.htriple1->SetTitle("Amplitudes of waveforms on the first FEB");
    h.htriple1->Draw("HIST");
    c1->SaveAs("hist1.pdf");
    delete c1;
    
    TCanvas *c2 = new TCanvas("c2", "htriple2", 900, 700);
    h.htriple2->SetLineColor(kGreen+1);
    h.htriple2->SetLineWidth(3);
    h.htriple2->SetFillColorAlpha(kGreen-4, 0.35);
    h.htriple2->SetTitle("Amplitudes of waveforms on the second FEB");
    h.htriple2->Draw("HIST");
    c2->SaveAs("hist2.pdf");
    delete c2;
    
    // Hit maps
    TCanvas *cmap1 = new TCanvas("cmap1", "Detector map", 800, 750);
    cmap1->SetRightMargin(0.15);
    gStyle->SetPalette(kViridis);
    h.mapdet1->Draw("COLZ TEXT");
    cmap1->SaveAs("hitmap1.pdf");
    delete cmap1;
}

void SaveResults(const string &outputFileName, Histograms &h,
                 const vector<vector<double>> &tabella1,
                 const vector<vector<double>> &tabella2) 
{
    // Calculate averages
    Double_t average_1[NCHANNELS], average_2[NCHANNELS];
    for (int i = 0; i < NCHANNELS; i++) 
    {
        average_1[i] = CalculateAverage(tabella1[i]);
        average_2[i] = CalculateAverage(tabella2[i]);
    }
    
    TFile *fout = new TFile(outputFileName.c_str(), "RECREATE");
    
    // Save all histograms
    h.hTime20->Write();
    h.hTime30->Write();
    h.hTime50->Write();
    h.hTime60->Write();
    h.hTime60_2->Write();
    h.htriple1->Write();
    h.htriple2->Write();
    h.htriple3->Write();
    h.htime_12->Write();
    h.htime20_12->Write();
    h.htime30_12->Write();
    h.hintegral1->Write();
    h.hintegral2->Write();
    h.hcharge1->Write();
    h.hcharge2->Write();
    h.hrisetime1->Write();
    h.hrisetime2->Write();
    h.hpeaktime1->Write();
    h.hpeaktime2->Write();
    h.timevsamplitude->Write();
    h.timevsamplitude2->Write();
    h.hbaseline1->Write();
    h.hbaseline2->Write();
    h.timevsintegral->Write();
    h.prof->Write();
    h.prof2->Write();
    h.cfdvschannel1->Write();
    h.conteggixcanale1->Write();
    h.conteggixcanale2->Write();
    h.conteggixcanale3->Write();
    h.mapdet1->Write();
    h.mapdet2->Write();
    
    // Save channel means tree
    TTree *tout = new TTree("channelMeans", "Medie per canale");
    tout->Branch("average_1", average_1, "average_1[64]/D");
    tout->Branch("average_2", average_2, "average_2[64]/D");
    tout->Fill();
    tout->Write();
    
    fout->Close();
    cout << "Analysis complete. Results saved to " << outputFileName << endl;
}