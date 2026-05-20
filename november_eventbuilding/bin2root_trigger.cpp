#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <cmath>
#include <array> 
#include <sstream>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>
#include <TApplication.h>  
//#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TProfile.h>
#include <TF1.h>
#include "TTreeIndex.h"
#include "TH1I.h" 

struct TriggerEntry 
{
    uint64_t TriggerIDSRS;
    double timestamp; // timestamp in ns
};

double getPeriodFactor(const std::string& settingsPath)
{
    std::ifstream file(settingsPath);
    std::string line;
    double samplingFrequency = 0;

    while (std::getline(file, line))
    {
        if (line.find("SamplingFrequency") != std::string::npos)
        {
            size_t pos = line.find(":");
            if (pos != std::string::npos)
            {
                std::string afterColon = line.substr(pos + 1);

                std::stringstream ss(afterColon);
                ss >> samplingFrequency;   // legge 8512
            }
            break;
        }
    }

    if (samplingFrequency == 0)
    {
        std::cerr << "SamplingFrequency not found!\n";
        return 1.0;
    }

    double periodFactor = (64.0 * 1000.0) / samplingFrequency;

    return periodFactor;
}

std::vector<TriggerEntry> readSAMPICTriggerBinary(const std::string& filename, int& run_number) 
{
    std::vector<TriggerEntry> entries;

    std::ifstream file(filename, std::ios::binary);
    if (!file) 
    {
        std::cerr << "Failed to open file: " << filename << "\n";
        return entries;
    }

    uint8_t TriggerIDFPGA;
    uint16_t TriggerIDSRSRaw;
    uint64_t timestampRaw;

    //double periodFactor = (64.0 * 1000.0) / 8512.0;

    uint64_t timestamp_prev = 0;
    uint64_t timestamp_ov = 0;

    uint16_t event_id_prev = 0;
    uint16_t event_id_ov = 0;

    while (file) 
    {
        // --- Leggi TriggerIDFPGA ---
        file.read(reinterpret_cast<char*>(&TriggerIDFPGA), sizeof(uint8_t));
        if (!file) break;

        // --- Leggi TriggerIDSRSRaw ---
        file.read(reinterpret_cast<char*>(&TriggerIDSRSRaw), sizeof(uint16_t));
        if (!file) break;

        // --- Leggi timestampRaw (5 byte) ---
        uint8_t buffer[5];
        file.read(reinterpret_cast<char*>(buffer), 5);
        if (!file) break;

        double periodFactor = getPeriodFactor("/home/riccardo-speziali/Scrivania/bin_file/Run"+std::to_string(run_number)+"_true/Run"+std::to_string(run_number)+"/sampic_run1/Run_Settings.txt");

        timestampRaw =
              (uint64_t)buffer[0]
            | ((uint64_t)buffer[1] << 8)
            | ((uint64_t)buffer[2] << 16)
            | ((uint64_t)buffer[3] << 24)
            | ((uint64_t)buffer[4] << 32);

        // --- Conta overflow timestamp ---
        if (timestampRaw < timestamp_prev) {
            timestamp_ov++;
        }
        timestamp_prev = timestampRaw;

        uint64_t timestamp_full = timestamp_ov * 1099511627776ULL + timestampRaw;
        double timestamp_ns = timestamp_full * periodFactor;

        // --- Conta overflow TriggerIDSRSRaw ---
        if (TriggerIDSRSRaw < (event_id_prev - 500)) {
            event_id_ov++;
        }
        event_id_prev = TriggerIDSRSRaw;
        uint64_t TriggerIDSRS = event_id_ov * 65536ULL + TriggerIDSRSRaw;

        // --- Salva entry ---
        entries.push_back({TriggerIDSRS, timestamp_ns});
    }

    file.close();
    return entries;
}

int main(int argc, char* argv[]) 
{
    if (argc < 1) 
    {
        std::cerr << "Usage: " << argv[0] << " <run_number> <subrun_number>\n";
        return 1;
    }

    int run_number = std::stoi(argv[1]);

    std::string filename = "/home/riccardo-speziali/Scrivania/bin_file/Run" + std::to_string(run_number) + "_true/Run" + std::to_string(run_number) + "/sampic_run1/sampic_run1_trigger_data.bin";
    auto entries = readSAMPICTriggerBinary(filename, run_number);
    std::cout << "Read " << entries.size() << " entries\n";

    // --- Crea un ROOT file ---
    std::string output_filename = "sampic_trigger_run" + std::to_string(run_number) + ".root";
    TFile *f = new TFile(output_filename.c_str(), "RECREATE");
    TTree *tree = new TTree("triggerTree", "SAMPIC Trigger Data");

    // Variabili per il TTree
    uint64_t TriggerIDSRS;
    double timestamp;

    // Branch
    tree->Branch("TriggerIDSRS", &TriggerIDSRS);
    tree->Branch("timestamp_ns", &timestamp);

    // Riempi il TTree
    for (auto &entry : entries) {
        TriggerIDSRS = entry.TriggerIDSRS;
        timestamp   = entry.timestamp;
        tree->Fill();
    }

    // Salva e chiudi
    tree->Write();
    f->Close();

    std::cout << "Saved data to sampic_trigger.root\n";

    return 0;
}