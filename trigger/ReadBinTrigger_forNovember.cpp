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


struct WaveformRecord {
    uint8_t  TriggerIDFPGA;
    uint16_t TriggerIDSRSRaw;
    uint64_t timestampRaw;
};



std::vector<WaveformRecord> read_waveform_file(const std::string& filename);

void ReadBinTrigger_forNovember() {

    std::string filename = "/home/riccardo-speziali/Scrivania/bin_file/Run222_true/Run222/sampic_run1/sampic_run1_trigger_data.bin";
    auto records = read_waveform_file(filename);

    
    TString outfile_name = "output_new.root";
    TFile outfile(outfile_name,"RECREATE");

    TTree *sampic_tree = new TTree("picoTree", "SAMPIC output in ROOT format");

    WaveformRecord rec;

   sampic_tree->Branch("TriggerIDFPGA", &rec.TriggerIDFPGA, "TriggerIDFPGA/B");
   sampic_tree->Branch("TriggerIDSRSRaw", &rec.TriggerIDSRSRaw, "TriggerIDSRSRaw/S");
   sampic_tree->Branch("timestampRaw", &rec.timestampRaw, "timestampRaw/L");


    for (const auto& r : records) {
        rec = r;
        sampic_tree->Fill();
    }

    sampic_tree->Write();
    outfile.Close();
}




std::vector<WaveformRecord> read_waveform_file(const std::string& filename) {

    std::vector<WaveformRecord> records;
    std::ifstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return records;
    }

    while (true) {
        WaveformRecord record;

        if (!file.read(reinterpret_cast<char*>(&record.TriggerIDFPGA), sizeof(UChar_t)))
            break;

        if (!file.read(reinterpret_cast<char*>(&record.TriggerIDSRSRaw), sizeof(UShort_t)))
            break;

        uint8_t buffer[5];
        if (!file.read(reinterpret_cast<char*>(buffer), 5))
            break;

        record.timestampRaw =
              (ULong64_t)buffer[0]
            | ((ULong64_t)buffer[1] << 8)
            | ((ULong64_t)buffer[2] << 16)
            | ((ULong64_t)buffer[3] << 24)
            | ((ULong64_t)buffer[4] << 32);

        records.push_back(record);
    }

    return records;
}
