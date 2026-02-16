#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include "TFile.h"
#include "TTree.h"


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


int main() 
{
    // Apri i file ROOT
    TFile *file_feb0 = OpenInputFile("/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb0_Corr.root");
    TFile *file_feb1 = OpenInputFile("/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb1_Corr.root");
    TFile *file_feb3 = OpenInputFile("/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run222/sampic_run1_feb3_Corr.root");
    TFile *trigger_file = OpenInputFile("/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/trigger/sampic_trigger_run222.root");

    if (!file_feb0 || !file_feb1 || !file_feb3 || !trigger_file) 
    {
        cerr << "Error: one or more files could not be opened." << endl;
        return 1;
    }

    // Qui puoi aggiungere il codice per leggere i dati dai file e analizzarli


    TTree *tree_feb0 = (TTree*)file_feb0->Get("picoTreewithCorr");
    TTree *tree_feb1 = (TTree*)file_feb1->Get("picoTreewithCorr");
    TTree *tree_feb3 = (TTree*)file_feb3->Get("picoTreewithCorr");
    TTree *trigger_tree = (TTree*)trigger_file->Get("triggerTree");

    if (!tree_feb0 || !tree_feb1 || !tree_feb3 || !trigger_tree) {
        cerr << "Error: picoTree or triggerTree not found in one or more files." << endl;
        return 1;
    }


    long nentries_feb0 = tree_feb0->GetEntries();
    long nentries_feb1 = tree_feb1->GetEntries();
    long nentries_feb3 = tree_feb3->GetEntries();
    long nentries_trigger = trigger_tree->GetEntries();


    for (int i = 0; i < nentries_trigger; i++)
        for(int j = 0; j < nentries_feb0; j++){
            if()
        }







            for(int k = 0; k < nentries_feb1; k++)
                for(int l = 0; l < nentries_feb3; l++)
                {
                    
                }
    

    





    
    // Chiudi i file alla fine
    file_feb0->Close();
    file_feb1->Close();
    file_feb3->Close();

    return 0;
}