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

struct building{
    double Cell0timestamp_MCP;
    double Cell0timeSTamp_PICOSEC[50];
    int chanel_PICOSEC[50];
    float TOTValue;
    float Waveform_MCP[64];
    float Waveform_PICOSEC[50][64];
    int hit_x_event;
    int SRS;


};

struct GEM_event {
    int SRS_trigger_ctr;
    int SRS_Timestamp;
    int SRS_Timestamp_ns;
    int n_tracks;
    int track_num;
    double track_chi_2;
    int ndetsintrack;
    vector<vector<double>> hits; // [n_tracks][ndetsintrack]
    vector<double> distnextcluster;
    vector<double> totchanexcluster;
};


int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run_number>\n";
        return 1;
    }
    int run_number = std::stoi(argv[1]);
    //int subrun_number = std::stoi(argv[2]);

    TString filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/run" + std::to_string(run_number) + "/eventbuilding.root";
    TString trackdata_filename = "/home/riccardo-speziali/Scrivania/bin_file/Run"+std::to_string(run_number)+"_true/Run"+std::to_string(run_number)+"/anaRun0"+std::to_string(run_number)+".root";
    TString output_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/run" + std::to_string(run_number) + "/eventbuilding_withtracks_run"+std::to_string(run_number)+".root";

    TFile *file = OpenInputFile(filename.Data());
    if (!file) {
        cerr << "Error: cannot open " << filename << endl;
        return 1;
    }   
    TTree *tree = (TTree*)file->Get("eventbuilding");
    if (!tree) {
        cerr << "Error: eventbuilding tree not found in " << filename << endl;
        return 1;
    }

        TFile *trackdata_file = OpenInputFile(trackdata_filename.Data());
    if (!trackdata_file) {
        cerr << "Error: cannot open " << trackdata_filename << endl;
        return 1;
    }
    TTree *trackdata_tree = (TTree*)trackdata_file->Get("tracks");
    if (!trackdata_tree) {
        cerr << "Error: trackTree not found in " << trackdata_filename << endl;
        return 1;
    }

    int SRS_mcp=0;
    int SRS_track=0;

    // Leggi i dati dal TTree

    building eventwotrack;
    tree->SetBranchAddress("Cell0TimeStamp_corr_MCP", &eventwotrack.Cell0timestamp_MCP);
    tree->SetBranchAddress("TriggerIDSRS_MCP", &eventwotrack.SRS);
    tree->SetBranchAddress("TOTValue_MCP", &eventwotrack.TOTValue);
    tree->SetBranchAddress("Waveform_MCP", &eventwotrack.Waveform_MCP);
    tree->SetBranchAddress("hitxevent", &eventwotrack.hit_x_event);
    tree->SetBranchAddress("Cell0TimeStamp_PICOSEC", eventwotrack.Cell0timeSTamp_PICOSEC);
    tree->SetBranchAddress("Channel_PICOSEC", eventwotrack.chanel_PICOSEC);
    tree->SetBranchAddress("Waveform_PICOSEC", eventwotrack.Waveform_PICOSEC); 


    GEM_event trackevent;
    trackdata_tree->SetBranchAddress("srstriggerctr", &trackevent.SRS_trigger_ctr);
    trackdata_tree->SetBranchAddress("srstimestamp", &trackevent.SRS_Timestamp);
    trackdata_tree->SetBranchAddress("srstimestampnsec", &trackevent.SRS_Timestamp_ns);
    trackdata_tree->SetBranchAddress("ntracks", &trackevent.n_tracks);
    trackdata_tree->SetBranchAddress("tracknumber", &trackevent.track_num);
    trackdata_tree->SetBranchAddress("trackchi2", &trackevent.track_chi_2);
    trackdata_tree->SetBranchAddress("ndetsintrack", &trackevent.ndetsintrack);
    //trackdata_tree->SetBranchAddress("hits", &trackevent.hits);
    //trackdata_tree->SetBranchAddress("distnextcluster", &trackevent.distnextcluster);
    //trackdata_tree->SetBranchAddress("totchanexcluster", &trackevent.totchanexcluster);

    TFile *output = new TFile(output_filename.Data(), "RECREATE");
    TTree *output_tree = new TTree("eventbuilding_withtracks", "eventbuilding_withtracks");
    output_tree->Branch("Cell0timestamp_MCP", &eventwotrack.Cell0timestamp_MCP,"Cell0timestamp_MCP/D");
    output_tree->Branch("TriggerIDSRS_MCP", &eventwotrack.SRS,"TriggerIDSRS_MCP/I");
    output_tree->Branch("SRS_trigger_ctr", &trackevent.SRS_trigger_ctr,"SRS_trigger_ctr/I");
    output_tree->Branch("SRS_Timestamp", &trackevent.SRS_Timestamp,"SRS_Timestamp/I");
    output_tree->Branch("SRS_Timestamp_ns", &trackevent.SRS_Timestamp_ns,"SRS_Timestamp_ns/I");
    output_tree->Branch("n_tracks", &trackevent.n_tracks,"n_tracks/I");
    output_tree->Branch("track_num", &trackevent.track_num,"track_num/I");
    output_tree->Branch("track_chi_2", &trackevent.track_chi_2,"track_chi_2/D");
    output_tree->Branch("ndetsintrack", &trackevent.ndetsintrack,"ndetsintrack/I");
    output_tree->Branch("hits", &trackevent.hits);
    output_tree->Branch("distnextcluster", &trackevent.distnextcluster);
    output_tree->Branch("totchanexcluster", &trackevent.totchanexcluster);
    output_tree->Branch("TOTValue", &eventwotrack.TOTValue,"TOTValue/F");
    output_tree->Branch("Waveform_MCP", &eventwotrack.Waveform_MCP,"Waveform_MCP[64]/F");
    output_tree->Branch("hit_x_event", &eventwotrack.hit_x_event,"hit_x_event/I");
    output_tree->Branch("Cell0timeSTamp_PICOSEC", eventwotrack.Cell0timeSTamp_PICOSEC,"Cell0timeSTamp_PICOSEC[50]/D");
    output_tree->Branch("chanel_PICOSEC", eventwotrack.chanel_PICOSEC,"chanel_PICOSEC[50]/I");
    output_tree->Branch("Waveform_PICOSEC", eventwotrack.Waveform_PICOSEC,"Waveform_PICOSEC[50][64]/F");      


    int j=0; // Indice per scorrere il TTree dei track
    int srs_mcp_prev=-1; // Variabile per tenere traccia dell'ultimo SRS MCP processato
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        // Ora puoi accedere ai dati di ogni evento tramite eventwotrack.Cell0timestamp
        
        // Trova l'evento di tracking corrispondente
        SRS_mcp = eventwotrack.SRS;
        bool track_found = false;
        if(SRS_mcp<srs_mcp_prev){
            j=0;
            cout << "Resetting track index to 0 because SRS_mcp is smaller than previous value." << endl;
        }
        srs_mcp_prev = SRS_mcp;
        //cout << "srs mcp: " << SRS_mcp << endl;
        
        while (j < trackdata_tree->GetEntries()) {
            trackdata_tree->GetEntry(j);
            SRS_track = trackevent.SRS_trigger_ctr;
            //cout << "srs track: " << SRS_track << endl;
            if(SRS_track < SRS_mcp) {
                j++;
                continue; // Continua a cercare
            }
            if(SRS_track > SRS_mcp) {
                break; // Non c'è corrispondenza, esci dal ciclo
            }   
            if (SRS_track == SRS_mcp) {
                track_found = true;
                output_tree->Fill();
                j++; // Incrementa j per non ripetere lo stesso track
                break;
            }
            if(i%10000==0){
               // cout << "Processing event " << i << "/" << tree->GetEntries() << "\r" << flush;
            }
            
        }
        if (!track_found) {
            //cout << "Warning: no track found for SRS " << SRS_mcp << endl;
        }
        
    }
//closw file
    output->Write();
    output->Close();
    file->Close();
    trackdata_file->Close();

    return 0;   
}

    
