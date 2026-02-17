#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include "TFile.h"
#include "TTree.h"
#include <cmath> 

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
    uint64_t TriggerIDSRS;
    
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




int main() 
{
    // Apri i file ROOT
    TString filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb1_Corr.root";
    TString filename_feb3 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb3_Corr.root";
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

    TTree *tree_feb1 = (TTree*)file_feb1->Get("picoTreewithCorr");
    TTree *tree_feb3 = (TTree*)file_feb3->Get("picoTreewithCorr");
    if (!tree_feb1 || !tree_feb3) {
        cerr << "Error: picoTreewithCorr not found in one or more files." << endl;
        return 1;
    }

    // Qui puoi aggiungere il codice per leggere i dati dai file e analizzarli

    WaveformRecord rec1, rec3;
    matchedEvent mcp_event;


    tree_feb1->SetBranchAddress("Cell0TimeStamp", &rec1.Cell0TimeStamp  );
    tree_feb1->SetBranchAddress("Cell0TimeStamp_corr", &rec1.Cell0TimeStamp_corr);
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    tree_feb1->SetBranchAddress("Channel", &rec1.channel);
    tree_feb1->SetBranchAddress("TOTValue", &rec1.TOTValue);
    tree_feb1->SetBranchAddress("TimeInstant", &rec1.TimeInstant);
    tree_feb1->SetBranchAddress("Baseline", &rec1.Baseline);
    tree_feb1->SetBranchAddress("PeakValue", &rec1.PeakValue);
    tree_feb1->SetBranchAddress("Amplitude", &rec1.Amplitude);
    tree_feb1->SetBranchAddress("Waveform", rec1.Waveform);

    tree_feb3->SetBranchAddress("Cell0TimeStamp", &rec3.Cell0TimeStamp);
    tree_feb3->SetBranchAddress("Cell0TimeStamp_corr", &rec3.Cell0TimeStamp_corr);
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    tree_feb3->SetBranchAddress("Channel", &rec3.channel);
    tree_feb3->SetBranchAddress("TOTValue", &rec3.TOTValue);
    tree_feb3->SetBranchAddress("TimeInstant", &rec3.TimeInstant);
    tree_feb3->SetBranchAddress("Baseline", &rec3.Baseline);
    tree_feb3->SetBranchAddress("PeakValue", &rec3.PeakValue);
    tree_feb3->SetBranchAddress("Amplitude", &rec3.Amplitude);
    tree_feb3->SetBranchAddress("Waveform", rec3.Waveform);

    matching_tree->SetBranchAddress("Cell0TimeStamp_corr", &mcp_event.Cell0TimeStamp_corr);
    matching_tree->SetBranchAddress("Channel", &mcp_event.channel);
    matching_tree->SetBranchAddress("TOTValue", &mcp_event.TOTValue);
    matching_tree->SetBranchAddress("TriggerIDSRS", &mcp_event.TriggerIDSRS);
    matching_tree->SetBranchAddress("Waveform", mcp_event.Waveform);    



    //output file
    TFile *output_file = new TFile(output_filename.Data(), "RECREATE");
    TTree *output_tree = new TTree("eventbuilding", "eventbuilding");
    //output_tree->Branch("Cell0TimeStamp", &rec0.Cell0TimeStamp,"Cell0TimeStamp/D");
    output_tree->Branch("Cell0TimeStamp_corr_MCP", &mcp_event.Cell0TimeStamp_corr,"Cell0TimeStamp_corr/D");
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    output_tree->Branch("Channel_MCP", &mcp_event.channel,"Channel_MCP/I");
    output_tree->Branch("TOTValue_MCP", &mcp_event.TOTValue,"TOTValue_MCP/F");
    output_tree->Branch("TriggerIDSRS", &mcp_event.TriggerIDSRS,"TriggerIDSRS_MCP/I");
    output_tree->Branch("Waveform_MCP", mcp_event.Waveform, "Waveform[64]/F");
    output_tree->Branch("Cell0TimeStamp_corr_FEB1", &rec1.Cell0TimeStamp_corr,"Cell0TimeStamp_corr_FEB1/D");
    output_tree->Branch("Channel_FEB1", &rec1.channel,"Channel_FEB1/I");
    output_tree->Branch("TOTValue_FEB1", &rec1.TOTValue,"TOTValue_FEB1/F");
    
    output_tree->Branch("Waveform_FEB1", rec1.Waveform, "Waveform[64]/F");
    output_tree->Branch("Cell0TimeStamp_corr_FEB3", &rec3.Cell0TimeStamp_corr,"Cell0TimeStamp_corr_FEB3/D");
    output_tree->Branch("Channel_FEB3", &rec3.channel,"Channel_FEB3/I");
    output_tree->Branch("TOTValue_FEB3", &rec3.TOTValue,"TOTValue_FEB3/F");
    
    output_tree->Branch("Waveform_FEB3", rec3.Waveform, "Waveform[64]/F");


    Long64_t nentries_feb1 = tree_feb1->GetEntries();
    Long64_t nentries_feb3 = tree_feb3->GetEntries();
    Long64_t nentries_matching = matching_tree->GetEntries();

    // Qui puoi aggiungere il codice per leggere i dati dai file e analizzarli
    
    Long64_t j1 = 0, j3 = 0;

    for (Long64_t i = 0; i < nentries_matching; i++) {

    matching_tree->GetEntry(i);

    




    





}
    
