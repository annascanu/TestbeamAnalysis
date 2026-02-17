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

struct TriggerEntry {
    uint64_t TriggerIDSRS;
    double timestamp; // timestamp in ns
};



int main() 
{
    int conta=0;
    // Apri i file ROOT
    TString filename_feb0 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb0_Corr.root";
    TString filename_feb1 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb1_Corr.root";
    TString filename_feb3 = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb3_Corr.root";
    TString filename_trigger = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/trigger/sampic_trigger_run222.root";
    TString output_filename = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/eventbuilder/MCPtoSRS_run222.root";
    TFile *file_feb0 = OpenInputFile(filename_feb0.Data());
    //TFile *file_feb1 = OpenInputFile(filename_feb1.Data());
    //TFile *file_feb3 = OpenInputFile(filename_feb3.Data());
    TFile *trigger_file = OpenInputFile(filename_trigger.Data());

    /*if (!file_feb0 || !file_feb1 || !file_feb3 || !trigger_file) 
    {
        cerr << "Error: one or more files could not be opened." << endl;
        return 1;
    }*/

    // Qui puoi aggiungere il codice per leggere i dati dai file e analizzarli


    TTree *tree_feb0 = (TTree*)file_feb0->Get("picoTreewithCorr");
    //TTree *tree_feb1 = (TTree*)file_feb1->Get("picoTreewithCorr");
    //TTree *tree_feb3 = (TTree*)file_feb3->Get("picoTreewithCorr");
    TTree *trigger_tree = (TTree*)trigger_file->Get("triggerTree");

    /*if (!tree_feb0 || !tree_feb1 || !tree_feb3 || !trigger_tree) {
        cerr << "Error: picoTree or triggerTree not found in one or more files." << endl;
        return 1;
    }*/

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

 /*   tree_feb1->Branch("Cell0TimeStamp", &rec1.Cell0TimeStamp,"Cell0TimeStamp/D");
    tree_feb1->Branch("Cell0TimeStamp_corr", &rec1.Cell0TimeStamp_corr,"Cell0TimeStamp_corr/D");
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    tree_feb1->Branch("Channel", &rec1.channel,"Channel/I");
    tree_feb1->Branch("TOTValue", &rec1.TOTValue,"TOTValue/F");
    tree_feb1->Branch("TimeInstant", &rec1.TimeInstant,"TimeInstant/D");
    tree_feb1->Branch("Baseline", &rec1.Baseline,"Baseline/F");
    tree_feb1->Branch("PeakValue", &rec1.PeakValue,"PeakValue/F");
    tree_feb1->Branch("Amplitude", &rec1.Amplitude,"Amplitude/F");
    tree_feb1->Branch("Waveform", rec1.Waveform, "Waveform[64]/F");

    tree_feb3->Branch("Cell0TimeStamp", &rec3.Cell0TimeStamp,"Cell0TimeStamp/D");
    tree_feb3->Branch("Cell0TimeStamp_corr", &rec3.Cell0TimeStamp_corr,"Cell0TimeStamp_corr/D");
    // tree_feb1->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    tree_feb3->Branch("Channel", &rec3.channel,"Channel/I");
    tree_feb3->Branch("TOTValue", &rec3.TOTValue,"TOTValue/F");
    tree_feb3->Branch("TimeInstant", &rec3.TimeInstant,"TimeInstant/D");
    tree_feb3->Branch("Baseline", &rec3.Baseline,"Baseline/F");
    tree_feb3->Branch("PeakValue", &rec3.PeakValue,"PeakValue/F");
    tree_feb3->Branch("Amplitude", &rec3.Amplitude,"Amplitude/F");
    tree_feb3->Branch("Waveform", rec3.Waveform, "Waveform[64]/F");
    */

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

   const double EPS = 1e10;   // tua finestra di matching

Long64_t j = 0;

for (Long64_t i = 0; i < nentries_trigger; i++) {

    trigger_tree->GetEntry(i);

    // Avanza feb0 finché è indietro
    while (j < nentries_feb0) {

        tree_feb0->GetEntry(j);

        double diff = rec0.Cell0TimeStamp_corr - trig.timestamp;

        if (fabs(diff) < EPS) {
            conta++;
            j++; 

            //aprire file output e inserire i dati di rec0, rec1 e rec3 in un unico albero con le stesse entry del trigger egrazie all'if legare srsID con gli eventi di feb0.
            
            output_tree->Fill();


            
        
            

            break;
        }

        if (rec0.Cell0TimeStamp_corr < trig.timestamp - EPS) {
            
            j++;           // feb0 troppo indietro → avanza feb0
        }
        else {
            
            break;         // feb0 è già avanti → passa al prossimo trigger
        }
    }

    cout << "events left: "
         << nentries_trigger - i
         << "\r" << flush;
}
            /*for(int k = 0; k < nentries_feb1; k++)
                for(int l = 0; l < nentries_feb3; l++)
                {
                    
                }*/
    

    





    cout<<"Numero di match trovati: " << conta << endl;
    // Chiudi i file alla fine
    file_feb0->Close();
    //file_feb1->Close();
    //file_feb3->Close();
    trigger_file->Close();
    output_file->cd();
    output_tree->Write();
    output_file->Close();

    return 0;
}