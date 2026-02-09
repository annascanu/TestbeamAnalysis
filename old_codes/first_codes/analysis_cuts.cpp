#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;

int main()
{
    // --------------------------------------------------------------------------
    //  Read parameters from the other code (from fitting of amplitude–integral)
    // --------------------------------------------------------------------------
    double a_cut = 0, b_cut = 0;

    ifstream fin("fit_parameters.txt");
    if (!fin.is_open()) {
        cout << "Cannot find an input file!" << endl;
        return 1;
    }

    fin >> a_cut >> b_cut;
    fin.close();
    cout << "Parameters: a = " << a_cut << " & b = " << b_cut << endl;


    // -----------------------------------------------------------------------
    //                        Opening .root file
    // -----------------------------------------------------------------------
    // TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root"; // Filename while running on lxplus
    TString filename = "/Users/anna/Developing/PhD/Testbeam2025/sampic_run22_final.root"; // Filename while running on Anna's machine
    // TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run19_final.root"; // Filename while running on Riccardo's machine
    cout << "Opening file: " << filename << endl;

    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) 
    {
        cout << "Cannot open ROOT file." << endl;
        return 1;
    }
    TTree *tree = (TTree*)file->Get("picoTree");
    if (!tree) {
        cerr << "Error: picoTree not found in " << filename << endl;
        return 1;
    }

    const int MAXPULSES = 200;
    Int_t npulses;
    Double_t pulses_amplitude[MAXPULSES], pulses_integral[MAXPULSES], pulses_time_20pe[MAXPULSES], pulses_t0[MAXPULSES];
    Bool_t pulses_bad_pulse[MAXPULSES];

    tree->SetBranchAddress("npulses", &npulses);
    tree->SetBranchAddress("pulses_amplitude", pulses_amplitude);
    tree->SetBranchAddress("pulses_integral", pulses_integral);
    tree->SetBranchAddress("pulses_time_20pe", pulses_time_20pe);
    tree->SetBranchAddress("pulses_t0", pulses_t0);
    tree->SetBranchAddress("pulses_bad_pulse", pulses_bad_pulse);

    // -----------------------------------------------------------------------
    //     3. ISTOGRAMMI (IDENTICI ALLA PRIMA MACRO)
    // -----------------------------------------------------------------------
    TH1D *hAmpAll_cut    = new TH1D("hAmpAll_cut", "Amplitude all pulses (after cut);Amplitude [mV]", 400, 0, 400);
    TH2D *hAmpVsIntegral_cut = new TH2D("hAmpVsIntegral_cut", "Amplitude vs Integral (after cut);Amp;Integral", 400, 0, 400, 200, 0, 1000);

    TH1D *hAmpSingle_cut = new TH1D("hAmpSingle_cut", "Amplitude single hits (after cut)",400, 0, 400);

    TH1D *hAmpTriple0_cut = new TH1D("hAmpTriple0_cut", "Triple hit pulse0 (after cut)", 400, 0, 400);
    TH1D *hAmpTriple1_cut = new TH1D("hAmpTriple1_cut", "Triple hit pulse1 (after cut)", 400, 0, 400);
    TH1D *hAmpTriple2_cut = new TH1D("hAmpTriple2_cut", "Triple hit pulse2 (after cut)", 400, 0, 400);

    // -----------------------------------------------------------------------
    //     Analysis with the applied cut
    // -----------------------------------------------------------------------

    int good = 0;
    double A, Q;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        good = 0;

        for (int j = 0; j < npulses; j++)
        {
            if (pulses_bad_pulse[j]) 
                continue;

            A = pulses_amplitude[j];
            Q = pulses_integral[j];

            // -----------------------------------------------
            //   Apply cut from amplitude–integral fitting
            // -----------------------------------------------
            if (Q > a_cut * A + b_cut)
                continue;

            good++;

            hAmpAll_cut->Fill(A);
            hAmpVsIntegral_cut->Fill(A, Q);
        }

        // --------------- single pulse -----------------
        if (good == 1)
        {
            for (int j=0; j<npulses; j++)
            {
                if (pulses_bad_pulse[j]) continue;

                A = pulses_amplitude[j];
                Q = pulses_integral[j];

                if (Q <= a_cut*A + b_cut)
                    hAmpSingle_cut->Fill(A);
            }
        }

        // ----------- triple hit ------------
        if (good == 3)
        {
            int k = 0;
            for (int j=0;j<npulses;j++)
            {
                if (pulses_bad_pulse[j])
                    continue;

                A = pulses_amplitude[j];
                Q = pulses_integral[j];

                if (Q > a_cut*A + b_cut) 
                    continue;

                if (k == 0) 
                    hAmpTriple0_cut->Fill(A);
                if (k == 1) 
                    hAmpTriple1_cut->Fill(A);
                if (k == 2) 
                    hAmpTriple2_cut->Fill(A);
                k++;
            }
        }
    }


    // ------------------------------------------------------------
    //            Save histograms in new root file
    // ------------------------------------------------------------
    TFile *fout = new TFile("analysis_cut.root","RECREATE");

    hAmpAll_cut->Write();
    hAmpVsIntegral_cut->Write();

    hAmpSingle_cut->Write();
    hAmpTriple0_cut->Write();
    hAmpTriple1_cut->Write();
    hAmpTriple2_cut->Write();

    fout->Close();

    cout << "\n----------------------------------------\n";
    cout << "   Output: analysis_cut.root\n";
    cout << "----------------------------------------\n";

    return 0;
}