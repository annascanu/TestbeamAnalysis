/*
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
                Analysis for the October 2025 ENUBET testbeam picosec time resolution.
                               Authors: A. Scanu, R. Speziali
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

Goal: obtain time resolution of picosec detectors.
To compile: c++ raw_analysis.cpp raw_analysis.cc `root-config --cflags --libs` -o raw.out

Processed files can be found at:
/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/
- Muon runs: run 19, 20, 21, 22
- 15 GeV pion runs: run 23, 24, ...
*/

#include "raw_analysis.h"
#include <iostream>
#include <TString.h>
#include <vector>

using namespace std;

int main() {
    // ------------------------------------------------
    //                Configuration
    // ------------------------------------------------
    
    // Choose your input file path
    //TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root"; // when running on lxplus
    //TString filename = "/Users/anna/Developing/PhD/Testbeam2025/sampic_run22_final.root"; // when running on Anna's machine
    //TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run22_final.root"; // when running on Riccardo's machine

    //const string filen = "/home/riccardo-speziali/Scrivania/november_2025/root_tree/sampic_run222_final.root";
    TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run22_final.root"; // when running on Riccardo's machine

    string outputFileName = "online_code_debug_1.root";

    // ------------------------------------------------
    //         Open input file and get tree
    // ------------------------------------------------
    
    TFile *file = OpenInputFile(filename.Data());
    if (!file) return 1;
    
    TTree *tree = (TTree*)file->Get("picoTree");
    if (!tree) {
        cerr << "Error: picoTree not found in " << filename << endl;
        return 1;
    }
    
    // ------------------------------------------------
    //         Initialize analysis structures
    // ------------------------------------------------
    
    Histograms hists;
    TreeBranches branches;
        
    InitializeHistograms(hists);
    SetupTreeBranches(tree, branches);
    
    // ------------------------------------------------
    //               Process all events
    // ------------------------------------------------
    
    cout << "\nProcessing events..." << endl;
    // 
    
    // ------------------------------------------------
    //                 Fit histograms
    // ------------------------------------------------
    
    // (LATER) Fit histograms to extract time resolution
    
    // ------------------------------------------------
    //            Create and save canvases
    // ------------------------------------------------
    
    cout << "\nCreating canvases..." << endl;
    //CreateCanvasesAndSaveResults();
    
    // ------------------------------------------------
    //        DEBUG: print some histogram contents to check
    // ------------------------------------------------

    int num_rows = 64;
    int info = 4;
    vector<vector<int>> coords;
    coords = ReadFile(); 

    cout << "\n\n\nEsempio di coordinate: Canale " << coords[2][0] << " x: " << coords[2][1] << " y: " << coords[2][2] << " board: " << coords[2][3] << "\n\n\n" << endl;
    ProcessEvents(tree, branches, hists, coords);
    CreateCanvasesAndSaveResults(outputFileName, hists);
        
    // ------------------------------------------------
    //                     Cleanup
    // ------------------------------------------------
    
    file->Close();
    delete file;
    
    cout << "\n=== Analysis completed successfully ===" << endl;
    
    return 0;
}
