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

struct WaveformRecord 
{
    double Cell0TimeStamp;
    double Cell0TimeStamp_corr; // added for corrected timestamp
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

struct TriggerEntry 
{
    uint64_t TriggerIDSRS;
    double timestamp; // timestamp in ns
};

int main(int argc, char* argv[]) 
{
    int run_number = std::stoi(argv[1]);
    //int subrun_number = std::stoi(argv[2]);

    int conta=0;
    // Apri i file ROOT
    TString filename_feb0 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run" + std::to_string(run_number) + "/sampic_run1_feb0_Corr.root";
    TString filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run" + std::to_string(run_number) + "/sampic_run1_feb1_Corr.root";
    TString filename_feb3 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run" + std::to_string(run_number) + "/sampic_run1_feb3_Corr.root";
    TString filename_trigger = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/trigger/sampic_trigger_run" + std::to_string(run_number) + ".root";
    TString output_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/MCPtoSRS_run" + std::to_string(run_number) + ".root";
    TFile *file_feb0 = OpenInputFile(filename_feb0.Data());
    //TFile *file_feb1 = OpenInputFile(filename_feb1.Data());
    //TFile *file_feb3 = OpenInputFile(filename_feb3.Data());
    TFile *trigger_file = OpenInputFile(filename_trigger.Data());

    TTree *tree_feb0 = (TTree*)file_feb0->Get("picoTreewithCorr");
    //TTree *tree_feb1 = (TTree*)file_feb1->Get("picoTreewithCorr");
    //TTree *tree_feb3 = (TTree*)file_feb3->Get("picoTreewithCorr");
    TTree *trigger_tree = (TTree*)trigger_file->Get("triggerTree");

    WaveformRecord rec0, rec1, rec3;
    TriggerEntry trig;

    tree_feb0->SetBranchAddress("Cell0TimeStamp", &rec0.Cell0TimeStamp);
    tree_feb0->SetBranchAddress("Cell0TimeStamp_corr", &rec0.Cell0TimeStamp_corr);
    // tree_feb1->SetBranchAddress("UnixTime", &rec.UnixTime);
    tree_feb0->SetBranchAddress("Channel", &rec0.channel);
    tree_feb0->SetBranchAddress("TOTValue", &rec0.TOTValue);
    tree_feb0->SetBranchAddress("TimeInstant", &rec0.TimeInstant);
    tree_feb0->SetBranchAddress("Baseline", &rec0.Baseline);
    tree_feb0->SetBranchAddress("PeakValue", &rec0.PeakValue);
    tree_feb0->SetBranchAddress("Amplitude", &rec0.Amplitude);
    tree_feb0->SetBranchAddress("Waveform", rec0.Waveform);

    trigger_tree->SetBranchAddress("TriggerIDSRS", &trig.TriggerIDSRS);
    trigger_tree->SetBranchAddress("timestamp_ns", &trig.timestamp);

    TFile *output_file = new TFile(output_filename.Data(), "RECREATE");
    TTree *output_tree = new TTree("eventTree", "Combined_Event_Data");
    //output_tree->Branch("Cell0TimeStamp", &rec0.Cell0TimeStamp,"Cell0TimeStamp/D");
    output_tree->Branch("Cell0TimeStamp_corr", &rec0.Cell0TimeStamp_corr,"Cell0TimeStamp_corr/D");
    output_tree->Branch("TimeInstant", &rec0.TimeInstant,"TimeInstant/D");
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    output_tree->Branch("Channel", &rec0.channel,"Channel/I");
    output_tree->Branch("TOTValue", &rec0.TOTValue,"TOTValue/F");
    output_tree->Branch("TriggerIDSRS", &trig.TriggerIDSRS,"TriggerIDSRS/I");
    output_tree->Branch("Waveform", rec0.Waveform, "Waveform[64]/F");

    Long64_t nentries_feb0 = tree_feb0->GetEntries();
    //Long64_t nentries_feb1 = tree_feb1->GetEntries();
    //Long64_t nentries_feb3 = tree_feb3->GetEntries();
    Long64_t nentries_trigger = trigger_tree->GetEntries();

    const double EPS = 1e10;   // matching window in ns (10 microseconds)
    Long64_t j = 0;  // index for feb0

    for (Long64_t i = 0; i < nentries_trigger; i++) 
    {
        trigger_tree->GetEntry(i);

        // Avanza feb0 finché è indietro
        while (j < nentries_feb0) 
        {
            tree_feb0->GetEntry(j);
            double diff = rec0.Cell0TimeStamp_corr - trig.timestamp;

            if (fabs(diff) < EPS) 
            {
                conta++;
                j++; 

                // open output file and insert data of rec0, rec1 e rec3 in one TTree with same entries
                // of the trigger; thanks to if statement, bind srsID to events in feb0.
                output_tree->Fill();
                break;
            }

            if (rec0.Cell0TimeStamp_corr < trig.timestamp - EPS) 
            {
                j++; // feb0 is behind, advance 
            }
            else 
            {    
                break; // feb0 is ahead, wait for next trigger
            }
        }

        cout << "events left: " << nentries_trigger - i << "\r" << flush;
    }

    // cout<<"Numero di match trovati: " << conta << endl;

    file_feb0->Close();
    trigger_file->Close();
    output_file->cd();
    output_tree->Write();
    output_file->Close();

    return 0;
}