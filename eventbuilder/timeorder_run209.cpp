#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <array>
#include <cmath>
#include <filesystem>
#include <string>

#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TProfile.h>
#include <TF1.h>
#include "TTreeIndex.h"

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

struct WaveformRecord 
{
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

struct EntryRef 
{
    double time;
    Long64_t entry;
};

int main(int argc, char* argv[]) 
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <run_number>\n";
        return 1;
    }

    int run_number = std::stoi(argv[1]);
    std::filesystem::create_directories("/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/run" + std::to_string(run_number));
    TString filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run" + std::to_string(run_number) + "/sampic_run1_feb1_Corr.root";
    TFile* file_feb1 = OpenInputFile(filename_feb1.Data());
   
    if (!file_feb1) 
    {
        cerr << "Error: could not open file for FEB1." << endl;
        return 1;
    }
    
    TTree *tree_feb1 = (TTree*)file_feb1->Get("picoTreewithCorr");
    
    WaveformRecord rec1, rec3;

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

    //output file: ordered by Cell0TimeStamp_corr
    TString output_filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/run" + std::to_string(run_number) + "/ordered_feb1.root";
    TFile *output_file = new TFile(output_filename_feb1.Data(), "RECREATE");
    TTree *output_tree = new TTree("eventbuilding", "eventbuilding");
    output_tree->Branch("Cell0TimeStamp_corr_FEB1", &rec1.Cell0TimeStamp_corr,"Cell0TimeStamp_corr_FEB1/D");
    output_tree->Branch("Channel_FEB1", &rec1.channel,"Channel_FEB1/I");
    output_tree->Branch("TOTValue_FEB1", &rec1.TOTValue,"TOTValue_FEB1/F");
    output_tree->Branch("TimeInstant_FEB1", &rec1.TimeInstant,"TimeInstant_FEB1/D");
    output_tree->Branch("Baseline_FEB1", &rec1.Baseline,"Baseline_FEB1/F");
    output_tree->Branch("PeakValue_FEB1", &rec1.PeakValue,"PeakValue_FEB1/F");
    output_tree->Branch("Amplitude_FEB1", &rec1.Amplitude,"Amplitude_FEB1/F");
    output_tree->Branch("Waveform_FEB1", rec1.Waveform, "Waveform_FEB1[64]/F");
    
    Long64_t nentries_feb1 = tree_feb1->GetEntries();

    std::vector<EntryRef> feb1_entries;
    feb1_entries.reserve(nentries_feb1);

    // Costruisci vettore (timestamp, entry_number)
    for (Long64_t i = 0; i < nentries_feb1; i++) 
    {
        tree_feb1->GetEntry(i);

        EntryRef e;
        e.time  = rec1.Cell0TimeStamp_corr;
        e.entry = i;

        feb1_entries.push_back(e);
    }

    // Ordina per timestamp
    std::sort(feb1_entries.begin(), feb1_entries.end(),
            [](const EntryRef &a, const EntryRef &b) {
                return a.time < b.time;
            });

    // Scrivi tree ordinato
    for (const auto &e : feb1_entries) 
    {
        tree_feb1->GetEntry(e.entry);
        output_tree->Fill();
    }

    double prev = -1e20;
    for (int i=0; i<100; i++) 
    {
        output_tree->GetEntry(i);
        if (rec1.Cell0TimeStamp_corr < prev)
            std::cout << "Not ordered!" << std::endl;
        prev = rec1.Cell0TimeStamp_corr;
    }

    output_file->cd();
    output_tree->Write();
    output_file->Close();

}