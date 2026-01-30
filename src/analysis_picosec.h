#ifndef ANALYSIS_H
#define ANALYSIS_H

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

// Constants
const int MAXPULSES = 200;
const int NCHANNELS = 64;
const int NDETECTORS = 3;

// Struct to hold all histograms
struct Histograms {
    // Basic histograms
    TH1F *hAmpAll;
    TH1F *hAmp1Hit;
    TH1F *hQ;
    TH1F *hTime20, *hTime30, *hTime50, *hTime60, *hTime30_2;
    
    // Triple hit amplitudes
    TH1F *htriple1, *htriple2, *htriple3;
    
    // Time differences
    TH1F *htime_12, *htime_13, *htime_23;
    TH1F *htime20_12, *htime20_13, *htime20_23;
    TH1F *htime30_12, *htime30_13, *htime30_23;
    
    // Charge, risetime, peaktime, baseline
    TH1F *hcharge1, *hcharge2;
    TH1F *hrisetime1, *hrisetime2;
    TH1F *hbaseline1, *hbaseline2;
    TH1F *hpeaktime1, *hpeaktime2;
    
    // Channel counts
    TH1F *conteggixcanale1, *conteggixcanale2, *conteggixcanale3;
    
    // Hit maps
    TH2F *mapdet1, *mapdet2, *mapdet1_afterselection, *mapdet2_afterselection;
    
    // Graphs
    TGraph *hintegral1, *hintegral2, *hintegral3;
    TGraph2D *timevsamplitude, *timevsamplitude2;
    TGraph *timevsintegral, *timevsintegral2;
    TGraph *cfdvschannel1;
    
    // Profiles
    TProfile *prof, *prof2;
};

// Struct to hold tree branches
struct TreeBranches {
    Double_t pulses_amplitude[MAXPULSES];
    Double_t pulses_integral[MAXPULSES];
    Double_t pulses_time_cfd10[MAXPULSES];
    Double_t pulses_time_cfd20[MAXPULSES];
    Double_t pulses_time_cfd30[MAXPULSES];
    Double_t pulses_time_cfd50[MAXPULSES];
    Double_t pulses_time_cfd60[MAXPULSES];
    Double_t Baseline[MAXPULSES];
    Double_t pulses_peak_time[MAXPULSES];
    Double_t pulses_rise_time[MAXPULSES];
    Double_t pulses_channel_x[MAXPULSES];
    Long64_t pulses_channel_y[MAXPULSES];
    Bool_t pulses_bad_pulse[MAXPULSES];
    Int_t npulses;
    Int_t HitFeb[3];
    Int_t Board[MAXPULSES];
    Int_t Channel[MAXPULSES];
};


struct FitResults {
    double prof_p0, prof_p0_err;
    double prof_p1, prof_p1_err;
    double prof_chi2;
    int    prof_ndf;

    double prof2_p0, prof2_p0_err;
    double prof2_p1, prof2_p1_err;
    double prof2_chi2;
    int    prof2_ndf;
};





// Function declarations
void InitializeHistograms(Histograms &hists);
void SetupTreeBranches(TTree *tree, TreeBranches &branches);
double CalculateAverage(const std::vector<double>& v);
void ProcessEvents(TTree *tree, TreeBranches &branches, Histograms &hists, 
                   std::vector<std::vector<double>> &tabella1, 
                   std::vector<std::vector<double>> &tabella2);
void FitHistograms(Histograms &hists, FitResults &fit);
void CreateCanvases(Histograms &hists);
void SaveResults(const std::string &outputFileName, Histograms &hists, 
                 const std::vector<std::vector<double>> &tabella1,
                 const std::vector<std::vector<double>> &tabella2, FitResults &fit);
void time_res(TTree *tree, TreeBranches &branches, TTree *tree_correction);
void ProcessEvents3rd_picosec(TTree *tree, TreeBranches &branches, Histograms &hists, 
                   std::vector<std::vector<double>> &tabella1, 
                   std::vector<std::vector<double>> &tabella2); //farlo diventare un int invece che un void per avere il nome del file di output
void ProcessEvents_november(TTree *tree, TreeBranches &branches, Histograms &hists, 
                   std::vector<std::vector<double>> &tabella1, 
                   std::vector<std::vector<double>> &tabella2);
TFile* OpenInputFile(const std::string &filename);

#endif // ANALYSIS_H
