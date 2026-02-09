#ifndef RAWANALYSIS_H
#define RAWANALYSIS_H

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

using namespace std;

// Constants
const int MAXPULSES = 200;
const int NCHANNELS = 64;
const int NDETECTORS = 3;

// Struct to hold all histograms
struct Histograms 
{
    // All amplitudes
    TH1F *hAmpAll1, *hAmpAll2, *hAmpAll3;
    
    // Baseline
    TH1F *hBaseline1, *hBaseline2, *hBaseline3;
    
    // Channel counts
    TH1F *hChannelHits1, *hChannelHits2, *hChannelHits3;
    
    // Amplitudes per channel
    // TH1F *hAmpChannel1[NCHANNELS], *hAmpChannel2[NCHANNELS], *hAmpChannel3[NCHANNELS];

    // Hit maps
    TH2F *mapDet1, *mapDet2, *mapDet3;
};

// Struct to hold tree branches
struct TreeBranches 
{
    Int_t ArraySize;
    Int_t Board[MAXPULSES];
    Int_t HitFeb[3];
    Int_t Channel[MAXPULSES];

    Double_t Cell0TimeStamp[MAXPULSES];
    // Double_t UnixTime[MAXPULSES];
    Double_t TimeInstant[MAXPULSES];

    Float_t Amplitude[MAXPULSES];
    Float_t Baseline[MAXPULSES];
    Float_t PeakValue[MAXPULSES];
    Float_t TOTValue[MAXPULSES];
    Float_t Waveform[MAXPULSES][NCHANNELS];
};

// Function declarations
// (LATER) Function to convert from binary to root
// We assume event building has already been done (https://gitlab.cern.ch/enubet-neutrino/testbeam2025/-/tree/master/Code/SAMPICRead?ref_type=heads)

TFile* OpenInputFile(const string &filename);
void InitializeHistograms(Histograms &h);
void SetupTreeBranches(TTree *tree, TreeBranches &b);
void ProcessEvents(TTree *tree, TreeBranches &b, Histograms &h, vector<vector<int>> &coordinates);
void CreateCanvasesAndSaveResults(const string &outputFileName, Histograms &hists);
vector<vector<int>> ReadFile();

// vogliamo vedere ampiezze (all)
// vedere le hit per canale, per capire se abbiamo canali rumorosi
// una per ampiezze per canale
// cell0timestamp, tempo assoluto del primo sample del sampic, cioè primo punto che acquisisce. da aggiungere poi al tempo
// del segnale per fare l'analisi temporale

// unix time: c'è un tempo di root e un tempo di unix, quindi devi fare una conversione temporale, altrimenti
// c'è qualche problema (chiediamo a Thomas? o vedi codice Alexandra?)

#endif // RAWANALYSIS_H
