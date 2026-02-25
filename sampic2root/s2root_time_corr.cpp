#include <iostream>
#include <fstream>
#include <string>
#include <vector>



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

std::vector<WaveformRecord> read_waveform_file(const std::string& filename);

void s2root_time_corr(int run_number) {
    std::string filename = "/home/riccardo-speziali/Scrivania/bin_file/Run" + std::to_string(run_number) + "_true/Run" + std::to_string(run_number) + "/sampic_run/feb0/sampic_run_feb0.bin";

    auto records = read_waveform_file(filename);

    std::cout << "Read " << records.size() << " waveform record(s)\n";


    ////////////////////////////////
    // -- Output to root TTree -- //
    ////////////////////////////////


    TString outfile_name = "/home/riccardo-speziali/Scrivania/git/TestbeamAnalysis/sampic2root/root_file/run" + std::to_string(run_number) + "/sampic_run1_feb0_Corr.root";
    TFile outfile(outfile_name,"RECREATE");
    TTree * sampic_tree = new TTree("picoTreewithCorr", "SAMPIC output in ROOT format");

    WaveformRecord rec;

    sampic_tree->Branch("Cell0TimeStamp", &rec.Cell0TimeStamp,"Cell0TimeStamp/D");
    sampic_tree->Branch("Cell0TimeStamp_corr", &rec.Cell0TimeStamp_corr,"Cell0TimeStamp_corr/D");
    // sampic_tree->Branch("UnixTime", &rec.UnixTime,"UnixTime/D");
    sampic_tree->Branch("Channel", &rec.channel,"Channel/I");
    sampic_tree->Branch("TOTValue", &rec.TOTValue,"TOTValue/F");
    sampic_tree->Branch("TimeInstant", &rec.TimeInstant,"TimeInstant/D");
    sampic_tree->Branch("Baseline", &rec.Baseline,"Baseline/F");
    sampic_tree->Branch("PeakValue", &rec.Baseline,"PeakValue/F");
    sampic_tree->Branch("Amplitude", &rec.Amplitude,"Amplitude/F");
    sampic_tree->Branch("Waveform", rec.Waveform, "Waveform[64]/F");

    for (size_t i = 0; i < records.size(); ++i) {
        rec = records[i];
        sampic_tree->Fill();

        // std::cout << "\nRecord #" << i + 1 << ":\n";
        // std::cout << "  Channel: " << rec.channel << "\n";
        // std::cout << "  TimeInstant: " << rec.TimeInstant << "\n";
        // std::cout << "  TOTValue: " << rec.TOTValue << "\n";
        // std::cout << "  Peak: " << rec.PeakValue << "\n";
        // std::cout << "  Amplitude: " << rec.Amplitude << "\n";
        // std::cout << "  DataSize: " << rec.DataSize << "\n";
        // std::cout << "  First 5 Samples: ";
        // for (int j = 0; j < std::min(5, rec.DataSize); ++j)
        // std::cout << rec.Waveform[j] << " ";
        // std::cout << "...\n";
    }


    sampic_tree->Write();





}



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



std::vector<WaveformRecord> read_waveform_file(const std::string& filename) {

    std::vector<WaveformRecord> records;
    std::ifstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return records;
    }
    const double TIMESTAMP_MAX = 1099511627776.0; // 2^40
   // const double periodFactor = (64*1000)/(samplingFrequency/1000000); // <-- metti il tuo valore corretto se necessario
    double periodFactor = getPeriodFactor("/home/riccardo-speziali/Scrivania/bin_file/Run222_true/Run222/sampic_run1/Run_Settings.txt");
    uint64_t timestamp_ov = 0;
    double timestamp_prev = 0;

    // file.seekg(4, std::ios::beg);

    

    while (file) {
        WaveformRecord record;

        // file.read(reinterpret_cast<char*>(&record.UnixTime), sizeof(double));
        file.read(reinterpret_cast<char*>(&record.channel), sizeof(int));
        file.read(reinterpret_cast<char*>(&record.Cell0TimeStamp), sizeof(double));
        file.read(reinterpret_cast<char*>(&record.TOTValue), sizeof(float));
        file.read(reinterpret_cast<char*>(&record.TimeInstant), sizeof(double));
        file.read(reinterpret_cast<char*>(&record.Baseline), sizeof(float));
        file.read(reinterpret_cast<char*>(&record.PeakValue), sizeof(float));
        file.read(reinterpret_cast<char*>(&record.Amplitude), sizeof(float));
        file.read(reinterpret_cast<char*>(&record.DataSize), sizeof(int));

        for(int i=0; i < record.DataSize; i++){
            float temp;
            file.read(reinterpret_cast<char*>(&temp), sizeof(float));
            record.Waveform[i] = temp;
        }

        if (!file) break;
        // ==============================
        //   TIMESTAMP OVERFLOW FIX
        // ==============================

        double timestampRaw = record.Cell0TimeStamp;

        if (timestampRaw < (timestamp_prev - 5000)) {
            timestamp_ov++;
        }
        timestamp_prev = timestampRaw;


       



        record.Cell0TimeStamp_corr = timestamp_ov * TIMESTAMP_MAX * periodFactor + timestampRaw;
       
        //record.Cell0TimeStamp;            
        
        //timestamp_ov * TIMESTAMP_MAX * periodFactor + timestampRaw;

        // ==============================

        //if(record.channel ==19) continue; //canale rumoroso
        if(record.channel !=0) continue; //per feb0 

        
        records.push_back(std::move(record));
       
    }

    return records;
}




