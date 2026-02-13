#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

using namespace std;

struct TriggerEntry {
    uint64_t TriggerIDSRS;
    double timestamp; // timestamp in ns
};

std::vector<TriggerEntry> readSAMPICTriggerBinary(const std::string& filename) {
    std::vector<TriggerEntry> entries;

    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return entries;
    }

    uint8_t TriggerIDFPGA;
    uint16_t TriggerIDSRSRaw;
    uint64_t timestampRaw;

    double periodFactor = (64.0 * 1000.0) / 8512.0;

    uint64_t timestamp_prev = 0;
    uint64_t timestamp_ov = 0;

    uint16_t event_id_prev = 0;
    uint16_t event_id_ov = 0;

    while (file) {
        // --- Leggi TriggerIDFPGA ---
        file.read(reinterpret_cast<char*>(&TriggerIDFPGA), sizeof(uint8_t));
        cout << "Read TriggerIDFPGA: " << static_cast<int>(TriggerIDFPGA) << std::endl;
        if (!file) break;

        // --- Leggi TriggerIDSRSRaw ---
        file.read(reinterpret_cast<char*>(&TriggerIDSRSRaw), sizeof(uint16_t));
        cout << "Read TriggerIDSRSRaw: " << TriggerIDSRSRaw << std::endl;
        if (!file) break;

        // --- Leggi timestampRaw (5 byte) ---
        uint8_t buffer[5];
        file.read(reinterpret_cast<char*>(buffer), 5);
        if (!file) break;

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

int main() {
    std::string filename = "/home/riccardo-speziali/Scrivania/bin_file/Run222_true/Run222/sampic_run1/sampic_run1_trigger_data.bin";
    auto entries = readSAMPICTriggerBinary(filename);

    std::cout << "Read " << entries.size() << " entries\n";
    for (size_t i = 0; i < std::min(entries.size(), size_t(10)); i++) {
        std::cout << "Entry " << i 
                  << ": TriggerIDSRS = " << entries[i].TriggerIDSRS 
                  << ", timestamp = " << entries[i].timestamp << " ns\n";
    }
}
