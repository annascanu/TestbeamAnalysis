#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>   // std::min, std::max
#include <cstdint>     // uint64_t
#include <atomic>      // std::memory_order
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Rtypes.h> 

using namespace std;





int main()
{

TString filename = "/home/riccardo-speziali/Scrivania/git/trigger_for_november/code/output_new.root";

cout << "Opening file: " << filename << endl;

    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: cannot open " << filename << endl;
        return 1;
    }
    TTree *tree = (TTree*)file->Get("picoTree");
    if (!tree) {
        cerr << "Error: picoTree not found in " << filename << endl;
        return 1;
    }




    // Variabili per leggere i branch
    Char_t TriggerIDFPGA;     // /B
    Short_t TriggerIDSRSRaw;   // /S
    Long64_t timestampRaw;     // /L

    // Collega i branch
    tree->SetBranchAddress("TriggerIDFPGA", &TriggerIDFPGA);
    tree->SetBranchAddress("TriggerIDSRSRaw", &TriggerIDSRSRaw);
    tree->SetBranchAddress("timestampRaw", &timestampRaw);

    // Numero di entry
    Long64_t nentries = tree->GetEntries();

    // Loop sugli eventi
   //       
    for (Long64_t i = 0; i < nentries; i++) {
        
        tree->GetEntry(i);
        if(i<10000 && i>9000){

        std::cout << "Entry " << i << std::endl;
        std::cout << "TriggerIDFPGA: " << (int)TriggerIDFPGA << std::endl;
        std::cout << "TriggerIDSRSRaw: " << TriggerIDSRSRaw << std::endl;
        std::cout << "timestampRaw: " << timestampRaw << std::endl;
        std::cout << "------------------------" << std::endl;
    }}

    file->Close();
    return 0;
}