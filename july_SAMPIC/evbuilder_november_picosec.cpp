#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>
#include <cmath> 
#include <array>

#include "TFile.h"
#include "TTree.h"

using namespace std;

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

struct matchedEvent 
{   
    double Cell0TimeStamp;
    int channel;
    float TOTValue;
    Int_t TriggerIDSRS;
    float Waveform[64];
};

struct WaveformRecord 
{
    double Cell0TimeStamp;
    double Cell0TimeStamp_corr;
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

struct building
{
    double Cell0timestamp_MCP;
    double Cell0timeSTamp_PICOSEC[50];
    int chanel_PICOSEC[50];
    float TOTValue_MCP;
    float TOTValue_PICOSEC[50];
    float Waveform_MCP[64];
    float Waveform_PICOSEC[50][64];
    int hit_x_event;
    int SRS;
};

int main(int argc, char* argv[]) 
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <run_number>\n";
        return 1;
    }

    int run_number = std::stoi(argv[1]);

    TString filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/july2026/eventbuilder/run" + std::to_string(run_number) + "/ordered_feb1.root";
    TString filename_feb3 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/july2026/eventbuilder/run" + std::to_string(run_number) + "/ordered_feb3.root";
    TString matching_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/july_SAMPIC/MCPtoSRS_run" + std::to_string(run_number) + ".root";
    TString output_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/run" + std::to_string(run_number) + "/eventbuilding.root";

    TFile *matching_file = OpenInputFile(matching_filename.Data());
    TFile *file_feb1 = OpenInputFile(filename_feb1.Data());
    TFile *file_feb3 = OpenInputFile(filename_feb3.Data());

    if (!matching_file || !file_feb1 || !file_feb3) {
        cerr << "Error: Impossible opening one or more input files." << endl;
        return 1;
    }

    TTree *matching_tree = (TTree*)matching_file->Get("eventTree");
    TTree *tree_feb1 = (TTree*)file_feb1->Get("eventbuilding");
    TTree *tree_feb3 = (TTree*)file_feb3->Get("eventbuilding");

    if (!matching_tree || !tree_feb1 || !tree_feb3) {
        cerr << "Error: Trees not found in input files." << endl;
        return 1;
    }

    WaveformRecord rec1, rec3;
    matchedEvent mcp_event;
    building built;

    tree_feb1->SetBranchAddress("Cell0TimeStamp_FEB1", &rec1.Cell0TimeStamp);
    tree_feb1->SetBranchAddress("Channel_FEB1", &rec1.channel);
    tree_feb1->SetBranchAddress("TOTValue_FEB1", &rec1.TOTValue);
    tree_feb1->SetBranchAddress("TimeInstant_FEB1", &rec1.TimeInstant);
    tree_feb1->SetBranchAddress("Baseline_FEB1", &rec1.Baseline);
    tree_feb1->SetBranchAddress("PeakValue_FEB1", &rec1.PeakValue);
    tree_feb1->SetBranchAddress("Amplitude_FEB1", &rec1.Amplitude);
    tree_feb1->SetBranchAddress("Waveform_FEB1", rec1.Waveform);

    tree_feb3->SetBranchAddress("Cell0TimeStamp_FEB3", &rec3.Cell0TimeStamp);
    tree_feb3->SetBranchAddress("Channel_FEB3", &rec3.channel);
    tree_feb3->SetBranchAddress("TOTValue_FEB3", &rec3.TOTValue);
    tree_feb3->SetBranchAddress("TimeInstant_FEB3", &rec3.TimeInstant);
    tree_feb3->SetBranchAddress("Baseline_FEB3", &rec3.Baseline);
    tree_feb3->SetBranchAddress("PeakValue_FEB3", &rec3.PeakValue);
    tree_feb3->SetBranchAddress("Amplitude_FEB3", &rec3.Amplitude);
    tree_feb3->SetBranchAddress("Waveform_FEB3", rec3.Waveform);

    matching_tree->SetBranchAddress("Cell0TimeStamp", &mcp_event.Cell0TimeStamp);
    matching_tree->SetBranchAddress("Channel", &mcp_event.channel);
    matching_tree->SetBranchAddress("TOTValue", &mcp_event.TOTValue);
    matching_tree->SetBranchAddress("TriggerIDSRS", &mcp_event.TriggerIDSRS);
    matching_tree->SetBranchAddress("Waveform", mcp_event.Waveform);    

    // Apri il file di output PRIMA di creare il TTree di output
    TFile *output_file = TFile::Open(output_filename.Data(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        cerr << "Error: cannot create output file " << output_filename << endl;
        return 1;
    }

    TTree *output_tree = new TTree("eventbuilding", "eventbuilding");

    output_tree->Branch("Cell0TimeStamp_MCP", &built.Cell0timestamp_MCP, "Cell0TimeStamp_MCP/D");
    output_tree->Branch("TOTValue_MCP", &built.TOTValue_MCP, "TOTValue_MCP/F");
    output_tree->Branch("TriggerIDSRS_MCP", &built.SRS, "TriggerIDSRS_MCP/I");
    output_tree->Branch("Waveform_MCP", built.Waveform_MCP, "Waveform_MCP[64]/F");
    output_tree->Branch("Cell0TimeStamp_PICOSEC", built.Cell0timeSTamp_PICOSEC, "Cell0TimeStamp_PICOSEC[hitxevent]/D");
    output_tree->Branch("Channel_PICOSEC", built.chanel_PICOSEC, "Channel_PICOSEC[hitxevent]/I");
    output_tree->Branch("TOTValue_PICOSEC", built.TOTValue_PICOSEC, "TOTValue_PICOSEC[hitxevent]/F");
    output_tree->Branch("Waveform_PICOSEC", built.Waveform_PICOSEC, "Waveform_PICOSEC[hitxevent][64]/F");
    output_tree->Branch("hitxevent", &built.hit_x_event, "hitxevent/I");

    Long64_t nentries_feb1 = tree_feb1->GetEntries();
    Long64_t nentries_feb3 = tree_feb3->GetEntries();
    Long64_t nentries_matching = matching_tree->GetEntries();

    double time_window = 50.0; // ns
    Long64_t j1_start = 0, j3_start = 0; 

    for (Long64_t i = 0; i < nentries_matching; i++) 
    {
        matching_tree->GetEntry(i);

        double t_mcp = mcp_event.Cell0TimeStamp;
        built.Cell0timestamp_MCP = t_mcp;
        built.TOTValue_MCP = mcp_event.TOTValue;
        built.SRS = mcp_event.TriggerIDSRS;

        for (int k = 0; k < 64; k++)         
            built.Waveform_MCP[k] = mcp_event.Waveform[k];
        
        int hitxevent = 0;

        // --- FEB1 Matching ---
        Long64_t j1 = j1_start;
        while (j1 < nentries_feb1) 
        {
            tree_feb1->GetEntry(j1);
            double dt = rec1.Cell0TimeStamp - t_mcp;

            if (dt < -time_window) 
            { 
                j1_start = j1 + 1; // Avanza l'indice di partenza per i prossimi eventi
                j1++; 
                continue; 
            } 
            if (dt > time_window) 
            { 
                break; // Troppo avanti nel tempo per questo t_mcp
            }

            if (hitxevent < 50) {
                built.Cell0timeSTamp_PICOSEC[hitxevent] = rec1.Cell0TimeStamp;
                built.chanel_PICOSEC[hitxevent] = rec1.channel;
                built.TOTValue_PICOSEC[hitxevent] = rec1.TOTValue;
                for (int k = 0; k < 64; k++)            
                    built.Waveform_PICOSEC[hitxevent][k] = rec1.Waveform[k];
                hitxevent++;
            }
            j1++;
        }

        // --- FEB3 Matching ---
        Long64_t j3 = j3_start;
        while (j3 < nentries_feb3) 
        {
            tree_feb3->GetEntry(j3);
            double dt = rec3.Cell0TimeStamp - t_mcp;
            
            if (dt < -time_window) 
            { 
                j3_start = j3 + 1; 
                j3++; 
                continue; 
            } 
            if (dt > time_window) 
                break;        

            if (hitxevent < 50) {
                built.Cell0timeSTamp_PICOSEC[hitxevent] = rec3.Cell0TimeStamp;
                built.chanel_PICOSEC[hitxevent] = rec3.channel + 64;
                built.TOTValue_PICOSEC[hitxevent] = rec3.TOTValue;
                for (int k = 0; k < 64; k++)
                    built.Waveform_PICOSEC[hitxevent][k] = rec3.Waveform[k];
                hitxevent++;
            }
            j3++;
        }

        built.hit_x_event = hitxevent;
        output_tree->Fill();

        if (i % 1000 == 0) 
        {
            cout << "Eventi analizzati: " << i << " / " << nentries_matching << "\r" << flush;
        }   
    }

    cout << "\nMatching completato. Salvataggio file..." << endl;

    // Scrittura fondamentale del file ROOT
    output_file->cd();
    output_tree->Write();
    output_file->Close();

    file_feb1->Close();
    file_feb3->Close();
    matching_file->Close();

    return 0;
}