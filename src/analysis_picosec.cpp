/*
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
                Analysis for the October 2025 ENUBET testbeam picosec time resolution.
                               Authors: A. Scanu, R. Speziali
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

Goal: obtain time resolution of picosec detectors.
To compile: c++ analysis_picosec.cpp analysis_picosec.cc `root-config --cflags --libs` -o analysis.out

Processed files can be found at:
/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/
- Muon runs: run 19, 20, 21, 22
- 15 GeV pion runs: run 23, 24, ...
*/

#include "analysis_picosec.h"
#include <iostream>
#include <TString.h>

using namespace std;

int main() {
    // ------------------------------------------------
    //                Configuration
    // ------------------------------------------------
    
    // Choose your input file path
    //TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root"; // when running on lxplus
    //TString filename = "/Users/anna/Developing/PhD/Testbeam2025/sampic_run22_final.root"; // when running on Anna's machine
    //TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run22_final.root"; // when running on Riccardo's machine

    TString filename = "/home/riccardo-speziali/Scrivania/november_2025/root_tree/sampic_run209_final.root"; // when running on Riccardo's machine

    string outputFileName = "muon_run209_final.root";

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
    FitResults fit;
    vector<vector<double>> tabella1(NCHANNELS);
    vector<vector<double>> tabella2(NCHANNELS);
    
    InitializeHistograms(hists);
    SetupTreeBranches(tree, branches);
    
    // ------------------------------------------------
    //               Process all events
    // ------------------------------------------------
    
    cout << "\nProcessing events..." << endl;
    ProcessEvents(tree, branches, hists, tabella1, tabella2);
    
    // ------------------------------------------------
    //                 Fit histograms
    // ------------------------------------------------
    
    cout << "\nFitting histograms..." << endl;
    FitHistograms(hists, fit);
    
    // ------------------------------------------------
    //            Create and save canvases
    // ------------------------------------------------
    
    cout << "\nCreating canvases..." << endl;
    CreateCanvases(hists);
    
    // ------------------------------------------------
    //        Save all results to output file
    // ------------------------------------------------
    
    cout << "\nSaving results..." << endl;
    SaveResults(outputFileName, hists, tabella1, tabella2, fit);
    
    // ------------------------------------------------
    //                     Cleanup
    // ------------------------------------------------
    
    file->Close();
    delete file;
    
    cout << "\n=== Analysis completed successfully ===" << endl;
    
    return 0;
}
