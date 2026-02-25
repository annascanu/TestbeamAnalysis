#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include "TFile.h"
#include "TTree.h"
#include <cmath> 
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <TROOT.h>
#include <TApplication.h>
#include <array>

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
#include "TTreeIndex.h"

#include "TApplication.h"
#include "TH1I.h"
#include "TCanvas.h"

#include "TApplication.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TSystem.h"




using namespace std;



//aprire i file .root con feb0(:contiene l'MCP) feb1 e feb3 che corrispondono ai dati del PICOSEC

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

struct matchedEvent {
   
    double Cell0TimeStamp_corr; // <-- aggiunta per il timestamp corretto
    int channel;
    
    float TOTValue;
    Int_t TriggerIDSRS;
    
    float Waveform[64];
};


struct WaveformRecord {
    double Cell0TimeStamp;
    double Cell0TimeStamp_corr; // <-- aggiunta per il timestamp corretto
    int channel;
    double UnixTime;
    float TOTValue;
    double TimeInstant;
    float Baseline;
    float PeakValue;
    float Amplitude;
    int DataSize;
    float Waveform[64];
};


struct building{
    double Cell0timestamp_MCP;
    double Cell0timeSTamp_PICOSEC[10];
    int chanel_PICOSEC[10];
    float TOTValue;
    float Waveform_MCP;
    float Waveform_PICOSEC[10][64];
    int hit_x_event;


};




int main() 
{
    double t_mcp=0;
    int channel_mcp=0;
    float tot_mcp=0;
    vector<float> waveform_temp;
    Int_t SRS;

    cout<<"FINE!"<<endl;
     // Dichiarazioni vettori per output
vector<double> cell0;
vector<int> channel_picosec;
vector<float> tot_picosec;
vector<float> waveform_picosec;
int hitxevent;  // nuovo branch

TGraph *cell0timehisto = new TGraph();
TGraph *dtmcp = new TGraph();

// Riserva memoria per ridurre reallocazioni
channel_picosec.reserve(5);
tot_picosec.reserve(5);
waveform_picosec.reserve(5);
cell0.reserve(5);


    // Apri i file ROOT
    TString filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/ordered_feb1.root";
    TString filename_feb3 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/ordered_feb3.root";
    TString matching_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/MCPtoSRS_run222.root";

    TString output_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/eventbuilding.root";
   
    //check
    TFile *matching_file = OpenInputFile(matching_filename.Data());
    if (!matching_file) {
        cerr << "Error: cannot open " << matching_filename << endl;
        return 1;
    }
    TTree *matching_tree = (TTree*)matching_file->Get("eventTree");
    if (!matching_tree) {
        cerr << "Error: eventTree not found in " << matching_filename << endl;
        return 1;
    }

    TFile *file_feb1 = OpenInputFile(filename_feb1.Data());
    TFile *file_feb3 = OpenInputFile(filename_feb3.Data());
    if (!file_feb1 || !file_feb3)     {
        cerr << "Error: cannot open one or more input files." << endl;
        return 1;
    }   

    TTree *tree_feb1 = (TTree*)file_feb1->Get("eventbuilding");
    TTree *tree_feb3 = (TTree*)file_feb3->Get("eventbuilding");
    if (!tree_feb1 || !tree_feb3) {
        cerr << "Error: picoTreewithCorr not found in one or more files." << endl;
        return 1;
    }

    // Qui puoi aggiungere il codice per leggere i dati dai file e analizzarli

    WaveformRecord rec1, rec3;
    matchedEvent mcp_event;
    building res;
        cout<<"FINE!"<<endl;


    //tree_feb1->SetBranchAddress("Cell0TimeStamp", &rec1.Cell0TimeStamp  );
    tree_feb1->SetBranchAddress("Cell0TimeStamp_corr_FEB1", &rec1.Cell0TimeStamp_corr);
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    tree_feb1->SetBranchAddress("Channel_FEB1", &rec1.channel);
    tree_feb1->SetBranchAddress("TOTValue_FEB1", &rec1.TOTValue);
    tree_feb1->SetBranchAddress("TimeInstant_FEB1", &rec1.TimeInstant);
    tree_feb1->SetBranchAddress("Baseline_FEB1", &rec1.Baseline);
    tree_feb1->SetBranchAddress("PeakValue_FEB1", &rec1.PeakValue);
    tree_feb1->SetBranchAddress("Amplitude_FEB1", &rec1.Amplitude);
    tree_feb1->SetBranchAddress("Waveform_FEB1", rec1.Waveform);

    //tree_feb3->SetBranchAddress("Cell0TimeStamp_FEB3", &rec3.Cell0TimeStamp);
    tree_feb3->SetBranchAddress("Cell0TimeStamp_corr_FEB3", &rec3.Cell0TimeStamp_corr);
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    tree_feb3->SetBranchAddress("Channel_FEB3", &rec3.channel);
    tree_feb3->SetBranchAddress("TOTValue_FEB3", &rec3.TOTValue);
    tree_feb3->SetBranchAddress("TimeInstant_FEB3", &rec3.TimeInstant);
    tree_feb3->SetBranchAddress("Baseline_FEB3", &rec3.Baseline);
    tree_feb3->SetBranchAddress("PeakValue_FEB3", &rec3.PeakValue);
    tree_feb3->SetBranchAddress("Amplitude_FEB3", &rec3.Amplitude);
    tree_feb3->SetBranchAddress("Waveform_FEB3", rec3.Waveform);

    matching_tree->SetBranchAddress("Cell0TimeStamp_corr", &mcp_event.Cell0TimeStamp_corr);
    matching_tree->SetBranchAddress("Channel", &mcp_event.channel);
    matching_tree->SetBranchAddress("TOTValue", &mcp_event.TOTValue);
    matching_tree->SetBranchAddress("TriggerIDSRS", &mcp_event.TriggerIDSRS);
    matching_tree->SetBranchAddress("Waveform", mcp_event.Waveform);    



    //output file
    TFile *output_file = new TFile(output_filename.Data(), "RECREATE");
    TTree *output_tree = new TTree("eventbuilding", "eventbuilding");
    //output_tree->Branch("Cell0TimeStamp", &rec0.Cell0TimeStamp,"Cell0TimeStamp/D");
    output_tree->Branch("Cell0TimeStamp_corr_MCP", &t_mcp,"Cell0TimeStamp_corr/D");
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
   output_tree->Branch("Channel_MCP", &channel_mcp,"Channel_MCP/I");
   output_tree->Branch("TOTValue_MCP", &tot_mcp,"TOTValue_MCP/F");
   output_tree->Branch("TriggerIDSRS_MCP", &SRS,"TriggerIDSRS_MCP/I");
   output_tree->Branch("Waveform_MCP", &waveform_temp);
    output_tree->Branch("Cell0TimeStamp_PICOSEC", &cell0);
    output_tree->Branch("Channel_PICOSEC", &channel_picosec);
    output_tree->Branch("TOTValue_PICOSEC", &tot_picosec);
    output_tree->Branch("Waveform_PICOSEC", &waveform_picosec);
    output_tree->Branch("hitxevent", &hitxevent);

        cout<<"FINE!"<<endl;

    //

    Long64_t nentries_feb1 = tree_feb1->GetEntries();
    Long64_t nentries_feb3 = tree_feb3->GetEntries();
    Long64_t nentries_matching = matching_tree->GetEntries();

    // Qui puoi aggiungere il codice per leggere i dati dai file e analizzarli
    
 

    cout<<"FINE!"<<endl;
    double tmcp_prevous = 0;
    double timefeb1=0, timefeb3=0;
  

double time_window = 50.0; // finestra temporale in ns
Long64_t j1 = 0, j3 = 0;   // indici globali per FEB1 e FEB3

int fakeargc = 0;
char** fakeargv = nullptr;
TApplication app("app", &fakeargc, fakeargv);

TH1I *hOcc = new TH1I("hOcc",
                      "Channel Occupancy (Live);Channel;Counts",
                      128, 0, 128);

TCanvas *cOcc = new TCanvas("cOcc","Live Occupancy",900,600);
hOcc->Draw();
cOcc->Update();

for (Long64_t i = 0; i < nentries_matching; i++) {

    matching_tree->GetEntry(i);

    // Pulizia vettori ad ogni evento MCP
    cell0.clear();
    channel_picosec.clear();
    tot_picosec.clear();
    waveform_picosec.clear();

    t_mcp = mcp_event.Cell0TimeStamp_corr;
    if(t_mcp < tmcp_prevous) {
        cout << "Warning: MCP timestamps not in order at entry " << i << endl;
    }
    //if(i%10000==0){cout<<"Entry: " << i<< "tmcp_prev: " << tmcp_prevous << endl;}
    tmcp_prevous = t_mcp;
    hitxevent = 0;

    //filling MCP Data
    channel_mcp = mcp_event.channel;
    //tot_mcp = mcp_event.TOTValue;
    waveform_temp.clear();
    for (int k = 0; k < 64; k++)
        waveform_temp.push_back(mcp_event.Waveform[k]);
    SRS = mcp_event.TriggerIDSRS;
    if(i%100000==0){
        cout<<"Entry: " << i<< "srs: " << SRS << endl;
    }

    // =========================
    // FEB1
    // =========================
    while (j1 < nentries_feb1) {
        tree_feb1->GetEntry(j1);
        double dt = rec1.Cell0TimeStamp_corr - t_mcp;
        if(rec1.Cell0TimeStamp_corr < timefeb1) {
            cout << "Warning: FEB1 timestamps not in order at entry " << j1 << endl;
        }
        timefeb1 = rec1.Cell0TimeStamp_corr;

        if (dt < -time_window) { j1++; continue; } // troppo vecchio, avanti
        if (dt > time_window) break;               // oltre finestra, stop loop

        // Hit valido
        cell0.push_back(rec1.Cell0TimeStamp_corr);
        channel_picosec.push_back(rec1.channel);
        tot_picosec.push_back(rec1.TOTValue);
        for (int k = 0; k < 64; k++)
            waveform_picosec.push_back(rec1.Waveform[k]);

        hitxevent++;
        j1++;
    }

    // =========================
    // FEB3
    // =========================
    while (j3 < nentries_feb3) {
        tree_feb3->GetEntry(j3);
        double dt = rec3.Cell0TimeStamp_corr - t_mcp;
        

        if (dt < -time_window) { j3++; continue; } // troppo vecchio
        if (dt > time_window) break;               // oltre finestra

        // Hit valido
        cell0.push_back(rec3.Cell0TimeStamp_corr);
        channel_picosec.push_back(rec3.channel + 63);
        tot_picosec.push_back(rec3.TOTValue);
        for (int k = 0; k < 64; k++)
            waveform_picosec.push_back(rec3.Waveform[k]);

        hitxevent++;
        j3++;
    }

    output_tree->Fill();


    // Riempimento live
if(hitxevent < 5 ){
for (auto ch : channel_picosec) {
    hOcc->Fill(ch);
}
}
// Aggiorna ogni 500 eventi (evita rallentamenti)
if (i % 500 == 0) {
    cOcc->Modified();
    cOcc->Update();
    gSystem->ProcessEvents();
}

    if (i % 1000 == 0) {
        cout << "Events left: " << nentries_matching - i 
             << " | Hits this event: " << hitxevent << "\r" << flush;
    }



    /*if(t_mcp> 8.75e12 && hitxevent>0) 
   // if(t_mcp<0.1e12 && hitxevent<4)
    { // esempio di condizione per debug
        cout << "Debug: MCP timestamp " << t_mcp << " at entry " << i << endl;
        cout << "  FEB1 hits in window: " << cell0.size() << endl;
        for (size_t idx = 0; idx < cell0.size(); idx++) {
            cout << "    Hit " << idx 
                 << ": Cell0TimeStamp_corr=" << cell0[idx] 
                 << ", Channel=" << channel_picosec[idx] 
                 << ", TOTValue=" << tot_picosec[idx] 
                 << endl;
        }
    }*/
}

cout << "Event building finito." << endl;
app.Run();



}