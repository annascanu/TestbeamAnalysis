#include "raw_analysis_april_tb.h"
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <TROOT.h>
#include <TApplication.h>
#include <array>
#include <fstream>
#include <sstream>
#include <vector>
#include <TSystem.h>

using namespace std;

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

void InitializeHistograms(Histograms &h) 
{
    // Basic histograms
    h.hAmpAll1 = new TH1F("hAmpAll1", "All pulse amplitudes for first detector;Amplitude [a.u.];Counts",  1000, 0, 1);
    h.hAmpAll2 = new TH1F("hAmpAll2", "All pulse amplitudes for second detector;Amplitude [a.u.];Counts", 1000, 0, 1);
    h.hAmpAll3 = new TH1F("hAmpAll3", "All pulse amplitudes for third detector;Amplitude [a.u.];Counts",  1000, 0, 1);
    
    // Baseline
    /* h.hBaseline1 = new TH1F("hBaseline1", "Baseline for first detector;Amplitude [a.u.];Counts",  1000, 0.9, 1.1);
    h.hBaseline2 = new TH1F("hBaseline2", "Baseline for second detector;Amplitude [a.u.];Counts", 1000, 0.9, 1.1);
    h.hBaseline3 = new TH1F("hBaseline3", "Baseline for third detector;Amplitude [a.u.];Counts",  1000, 0.9, 1.1);*/

    // Channel counts
    h.hChannelHits1 = new TH1F("hChannelHits1", "Channel hits for first detector;Channel;Counts", 96, 0, 96);
    h.hChannelHits2 = new TH1F("hChannelHits2", "Channel hits for second detector;Channel;Counts", 96, 0, 96);
    h.hChannelHits3 = new TH1F("hChannelHits3", "Channel hits for third detector;Channel;Counts", 96, 0, 96);

    // Amplitudes per channel
    
    for (int i = 0; i < 96; i++) 
    {
        h.hAmpChannel1[i] = new TH1F(Form("hAmpChannel1_%d", i), Form("Amplitudes for channel %d of first detector;Amplitude [a.u.];Counts", i), 100, 0, 0.5);
        h.hAmpChannel2[i] = new TH1F(Form("hAmpChannel2_%d", i), Form("Amplitudes for channel %d of second detector;Amplitude [a.u.];Counts", i), 100, 0, 0.5);
        h.hAmpChannel3[i] = new TH1F(Form("hAmpChannel3_%d", i), Form("Amplitudes for channel %d of third detector;Amplitude [a.u.];Counts", i), 100, 0, 0.5);
    }

    // Hit maps
    h.mapDet1 = new TH2F("hitmapdet1", "Hit map for detector 1;x;y", 10, 0, 10, 10, 0, 10);
    h.mapDet2 = new TH2F("hitmapdet2", "Hit map for detector 2;x;y", 10, 0, 10, 10, 0, 10);
    h.mapDet3 = new TH2F("hitmapdet3", "Hit map for detector 3;x;y", 30, -5, 5, 30, -5, 5);
}

void SetupTreeBranches(TTree *tree, TreeBranches &b) 
{
    cout << "\n\n                  Setting up tree branches...\n\n\n" << endl;
    // tree->Print();

    tree -> SetBranchAddress("ArraySize", &b.ArraySize);
    //tree -> SetBranchAddress("Board", &b.Board);
    tree -> SetBranchAddress("Feb", &b.Feb);
    tree -> SetBranchAddress("Channel", &b.Channel);

    tree -> SetBranchAddress("Cell0TimeStamp", &b.Cell0TimeStamp);
    // tree -> SetBranchAddress("UnixTime", &b.UnixTime);
    tree -> SetBranchAddress("TimeInstant", &b.TimeInstant);

    tree -> SetBranchAddress("Amplitude", &b.Amplitude);
    tree -> SetBranchAddress("Baseline", &b.Baseline);
    tree -> SetBranchAddress("PeakValue", &b.PeakValue);
    tree -> SetBranchAddress("TOTValue", &b.TOTValue);
    tree -> SetBranchAddress("Waveform", &b.Waveform);
    tree -> SetBranchAddress("xCoord", &b.xCoord);
    tree -> SetBranchAddress("yCoord", &b.yCoord);
    tree -> SetBranchAddress("Detector", &b.Detector);
}

void ProcessEvents(TTree *tree, TreeBranches &b, Histograms &h) ///aggiungere vettore di vettori
{   
    
    
    // Reference code for the polya fit on one single channel
    /*
    TF1 *fpolyaAmpl1 = new TF1("fpolyaAmpl1", 
        "([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])", 
        0, 0.5);
    //fpolyaAmpl1->SetParLimits(2, 0.05, 1.5);

    //fpolyaAmpl1->SetParameter(1, 180);
    fpolyaAmpl1->SetNpx(10000);
    fpolyaAmpl1->SetParameter(0,1000); // amplitude
    fpolyaAmpl1->SetParameter(1,30); //Charge
    fpolyaAmpl1->SetParameter(2,2); //theta  im putting theta on 2 from 30 for debugging
    fpolyaAmpl1->SetParLimits(2,0.001,5000); //theta
    
    //h.hcharge1->Fit(fpolyaAmpl1, "RQ");
    */

    Long64_t nentries = tree->GetEntries();
    cout << "Total number of entries: " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++) 
    {
        tree->GetEntry(i);
        
        // Basic histograms for amplitude
        for (int j = 0; j < b.ArraySize; j++) //arraysize : number of hits per event
        {
            //check
            
            if(b.ArraySize!=3)continue; //just to check the events with 3 hits, then we can remove this condition
            if (b.Detector[j] == 1)
            {
            //if(b.xCoord[j]==4 && b.yCoord[j]==3)continue;
            //if(b.xCoord[j]==1 && b.yCoord[j]==4)continue;             


                            
                /*if(b.Channel[j]>=31){
                        
                                if(b.Channel[j]>=47){
                                int diff=b.Channel[j]-63;
                                b.Channel[j]= b.Channel[j]-32+diff;
                                int coord_diff_x=b.xCoord[j]-9;
                                b.xCoord[j]=b.xCoord[j]-9+coord_diff_x;
                                if(b.yCoord[j]>=4){
                                int coord_diff_y=b.yCoord[j]-3;
                                b.yCoord[j]=b.yCoord[j]-3+coord_diff_y;
                                }
                                else{
                                int coord_diff_y=b.yCoord[j]-4;
                                b.yCoord[j]=b.yCoord[j]-4+coord_diff_y;
                                }
                                }
                                else{
                                int diff=b.Channel[j]-63;
                                b.Channel[j]=b.Channel[j]+32-diff;
                                }
                            }*/
                h.hAmpAll1->Fill(abs(b.Amplitude[j]));
                //cout << "Amplitude: " << b.Amplitude[j] << endl;
                //h.hBaseline1->Fill(b.Baseline[j]);
                //if(b.Feb[j]==3){}
                if(b.Feb[j]==1){
                    b.Channel[j]=b.Channel[j]+64;
                }



                h.hChannelHits1->Fill(b.Channel[j]);
                for(int k=0; k<96; k++){
                    if(b.Channel[j]==k){
                        h.hAmpChannel1[k]->Fill(abs(b.Amplitude[j]));
                        //h.hAmpChannel1[k]->Fit(&polya_fits[k], "RQ");
                    }
                }
                
            }
            else if (b.Detector[j] == 2)
            {
                h.hAmpAll2->Fill(abs(b.Amplitude[j]));
                //h.hBaseline2->Fill(b.Baseline[j]);ù
                if(b.Feb[j]==0){
                    b.Channel[j]=b.Channel[j]+32;
                }
                h.hChannelHits2->Fill(b.Channel[j]);
                for(int k=0; k<96; k++){
                    if(b.Channel[j]==k){
                        h.hAmpChannel2[k]->Fill(abs(b.Amplitude[j]));
                        //h.hampChannel2[k]->Fit(fpolyaAmpl1, "RQ");
                    }
                }
            }
            else if (b.Detector[j] == 3)
            {
                h.hAmpAll3->Fill(abs(b.Amplitude[j]));
                //h.hBaseline3->Fill(b.Baseline[j]);
                h.hChannelHits3->Fill(b.Channel[j]);
                for(int k=0; k<64; k++){
                    if(b.Channel[j]==k){
                        h.hAmpChannel3[k]->Fill(abs(b.Amplitude[j]));
                    }
                }
            }
            

        }
          /*  h.hAmpAll1->GetXaxis()->SetRange(
            h.hAmpAll1->FindFirstBinAbove(0),
            h.hAmpAll1->FindLastBinAbove(0)
            );
                    h.hAmpAll2->GetXaxis()->SetRange(
            h.hAmpAll2->FindFirstBinAbove(0),
            h.hAmpAll2->FindLastBinAbove(0)
            );
                    h.hAmpAll3->GetXaxis()->SetRange(
            h.hAmpAll3->FindFirstBinAbove(0),
            h.hAmpAll3->FindLastBinAbove(0)
            );*/
        //cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
        // Channels 
        //if (b.HitFeb[0] >= 1 && b.HitFeb[1] >= 1 && b.HitFeb[2] >= 1) //added the parto of hitfeb[2]
        //if(b.HitFeb[2]==1)//check con ultima analisis con trigger sul terzo det
        if(true)
        {
            double amp_FEB[3] = {-1.0, -1.0, -1.0};
            double base[3] = {-9999.0, -9999.0, -9999.0};
            

            for (int j = 0; j < b.ArraySize; j++) 
            {
                            if(b.ArraySize!=3)continue; //just to check the events with 3 hits, then we can remove this condition

                int det = b.Detector[j];
                //cout << "Board: " << det << endl;

               // if (det < 1 || det > 3) continue;
               
                // Just testing low amplitude cut
                //if (b.pulses_amplitude[j] < 0.05) continue;
                
                // amp_FEB[det] = b.pulses_amplitude[j];
                //base[det] = b.Baseline[j];
                // peak_time[det] = b.pulses_peak_time[j];

               // for(int k=0; k<64; k++){
                    //if(b.Channel[j]==k){
                        if(det==1){
                            /*if(b.xCoord[j]==5 && b.yCoord[j]==4)continue;
                            if(b.xCoord[j]==8 && b.yCoord[j]==5)continue;*/


                                //if(b.xCoord[j]==4 && b.yCoord[j]==3)continue;
                                //if(b.xCoord[j]==1 && b.yCoord[j]==4)continue;
                                h.mapDet1->Fill(b.xCoord[j], b.yCoord[j]);
                           
                            
                    }
                        else if(det==2){
                            h.mapDet2->Fill(b.xCoord[j], b.yCoord[j]);
                        }
                        else if(det==3){
                            h.mapDet3->Fill(b.xCoord[j], b.yCoord[j]);
                        }
                    //}
                //}


                //cout << "Channel: " << b.Channel[j] << endl;

                /*
                if (det == 0)
                    h.mapDet1->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);
                if (det == 1)
                    h.mapDet2->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);
                if (det == 2)
                    h.mapdet3->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);   for third detector*/
            }
            
           /* h.hBaseline1->Fill(base[0]);
            h.hBaseline2->Fill(base[1]);
            h.hBaseline3->Fill(base[2]); *///for the third det
            
        }
        
        if (i % 100000 == 0)
            cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }

    cout << endl;
    cout << "Finished processing all events." << endl;
    cout<<"till here fine"<<endl;

}

void CreateCanvasesAndSaveResults(const std::string &outputFileName, Histograms &hists, int &run)
{
        cout<<"till here fine"<<endl;

    TString outDir = Form("/home/riccardo-speziali/enubet/online_analysis/Results/run%d", run);
    gSystem->mkdir(outDir, kTRUE);  // kTRUE = crea anche path intermedi
    // ---------------------------- Create one PDF file for each detector ----------------------------
    TCanvas *cSavePDF_Detector1 = new TCanvas("cSavePDF", "Saving PDF", 1200, 800);
    cSavePDF_Detector1 -> SaveAs(Form("%s/Amplitude_Distribution_Detector1.pdf[", outDir.Data())); // Just open the PDF file for writing (append mode)
    TCanvas *cSavePDF_Detector2 = new TCanvas("cSavePDF", "Saving PDF", 1200, 800);
    cSavePDF_Detector2 -> SaveAs(Form("%s/Amplitude_Distribution_Detector2.pdf[", outDir.Data())); 

    TCanvas *cSavePDF_Detector3 = new TCanvas("cSavePDF", "Saving PDF", 1200, 800);
    cSavePDF_Detector3 -> SaveAs(Form("%s/Amplitude_Distribution_Detector3.pdf[", outDir.Data())); 
    cout<<"till here fine"<<endl;

    // ----------------------------------------------------------------------------------------------------------------
    
    // Canvas for amplitude distributions
    TCanvas *c1 = new TCanvas("c1", "Amplitude distribution for first detector", 1200, 800);
    hists.hAmpAll1->Draw();
    c1 -> SaveAs(Form("%s/Amplitude_Distribution_Detector1.pdf", outDir.Data())); // Save the first page of the PDF file

    TCanvas *c2 = new TCanvas("c2", "Amplitude distribution for second detector", 1200, 800);
    hists.hAmpAll2->Draw();
    c2 -> SaveAs(Form("%s/Amplitude_Distribution_Detector2.pdf", outDir.Data()));
    
    TCanvas *c3 = new TCanvas("c3", "Amplitude distribution for third detector", 1200, 800);
    hists.hAmpAll3->Draw();
    c3 -> SaveAs(Form("%s/Amplitude_Distribution_Detector3.pdf", outDir.Data()));
    
    c1->Update();
    c2->Update();
    c3->Update();

    cSavePDF_Detector1 -> SaveAs(Form("%s/Amplitude_Distribution_Detector1.pdf]", outDir.Data()));
    cSavePDF_Detector2 -> SaveAs(Form("%s/Amplitude_Distribution_Detector2.pdf]", outDir.Data()));
    cSavePDF_Detector3 -> SaveAs(Form("%s/Amplitude_Distribution_Detector3.pdf]", outDir.Data()));
    cout<<"till here fine"<<endl;

    // Canvas for hit maps
    //creazione unico pdf con tutti i canvas delle hitmap

    TCanvas *cSavePDF_HitMaps = new TCanvas("cSavePDF_HitMaps", "Saving PDF", 1200, 800);
    cSavePDF_HitMaps -> SaveAs(Form("%s/Hit_Maps.pdf[", outDir.Data()));

    TCanvas *cMap1 = new TCanvas("cMap1", "Hit map for DAVIDE detector", 1200, 800);
    hists.mapDet1->Draw("COLZ");
    cMap1 -> SaveAs(Form("%s/Hit_Maps.pdf", outDir.Data()));   

    TCanvas *cMap2 = new TCanvas("cMap2", "Hit map for GOLIA detector", 1200, 800);
    hists.mapDet2->Draw("COLZ");
    cMap2 -> SaveAs(Form("%s/Hit_Maps.pdf", outDir.Data()));

    TCanvas *cMap3 = new TCanvas("cMap3", "Hit map for 7PAD detector", 1200, 800);
    hists.mapDet3->Draw("COLZ");
    cMap3 -> SaveAs(Form("%s/Hit_Maps.pdf", outDir.Data()));

    cMap1->Update();
    cMap2->Update();    
    cMap3->Update();    

    cout<<"till here fine"<<endl;



    cSavePDF_HitMaps -> SaveAs(Form("%s/Hit_Maps.pdf]", outDir.Data()));




    //save channel hits
    TCanvas *cSavePDF_ChannelHits = new TCanvas("cSavePDF_ChannelHits", "Saving PDF", 1200, 800);
    cSavePDF_ChannelHits -> SaveAs(Form("%s/Channel_Hits.pdf[", outDir.Data()));

    TCanvas *cChannelHits1 = new TCanvas("cChannelHits1", "Channel Hits for first detector", 1200, 800);
    hists.hChannelHits1->Draw();
    cChannelHits1 -> SaveAs(Form("%s/Channel_Hits.pdf", outDir.Data()));

    TCanvas *cChannelHits2 = new TCanvas("cChannelHits2", "Channel Hits for second detector", 1200, 800);
    hists.hChannelHits2->Draw();
    cChannelHits2 -> SaveAs(Form("%s/Channel_Hits.pdf", outDir.Data()));

    TCanvas *cChannelHits3 = new TCanvas("cChannelHits3", "Channel Hits for third detector", 1200, 800);
    hists.hChannelHits3->Draw();
    cChannelHits3 -> SaveAs(Form("%s/Channel_Hits.pdf", outDir.Data()));

    cChannelHits1->Update();
    cChannelHits2->Update();
    cChannelHits3->Update();    
    cout<<"till here fine"<<endl;

    cSavePDF_ChannelHits -> SaveAs(Form("%s/Channel_Hits.pdf]", outDir.Data()));


    TCanvas *cSavePDF_AmplitudesPerChannel = new TCanvas("cSavePDF_AmplitudesPerChannel", "Saving PDF", 1200, 800);
    cSavePDF_AmplitudesPerChannel -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector1.pdf[", outDir.Data()));

    TCanvas *cSavePDF_AmplitudesPerChannel2 = new TCanvas("cSavePDF_AmplitudesPerChannel2", "Saving PDF", 1200, 800);
    cSavePDF_AmplitudesPerChannel2 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector2.pdf[", outDir.Data()));

    TCanvas *cSavePDF_AmplitudesPerChannel3 = new TCanvas("cSavePDF_AmplitudesPerChannel3", "Saving PDF", 1200, 800);
    cSavePDF_AmplitudesPerChannel3 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector3.pdf[", outDir.Data()));


       // Polya function
    std::vector<TF1*>* polya_fits = new std::vector<TF1*>();
    for (int i = 0; i < 96; i++) 
    {
        TF1 *fpolya = new TF1(Form("fpolya_%d", i), 
            "([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])", 
            0.02, 0.3);
        fpolya->SetNpx(10000);
        //fpolya->SetParameter(0,1000); // amplitude
        fpolya->SetParameter(1,0.1); //Charge
        fpolya->SetParLimits(1,0.001,0.2); //Charge
        fpolya->SetParameter(2,0.1); //theta
       // fpolya->SetParLimits(2,0.001,5000); //theta

        fpolya->SetParName(0,"c");
        fpolya->SetParName(1,"#bar{Q}");
        fpolya->SetParName(2,"#theta");
        polya_fits->push_back(fpolya);
    }

    for(int k=0; k<96; k++){
        /////polya fits for each channel to be added here, then save the canvases with the fits
     




        TCanvas *cAmpChannel1 = new TCanvas(Form("cAmpChannel1_%d", k), Form("Amplitudes for channel %d of first detector", k), 1200, 800);
        hists.hAmpChannel1[k]->Draw();
        hists.hAmpChannel1[k]->Fit(polya_fits->at(k), "RQ");
        //hists.hAmpChannel1[k]->Legend()->SetTextSize(0.02);
        cAmpChannel1 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector1.pdf", outDir.Data()));

        TCanvas *cAmpChannel2 = new TCanvas(Form("cAmpChannel2_%d", k), Form("Amplitudes for channel %d of second detector", k), 1200, 800);
        hists.hAmpChannel2[k]->Draw();
        cAmpChannel2 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector2.pdf", outDir.Data()));

        TCanvas *cAmpChannel3 = new TCanvas(Form("cAmpChannel3_%d", k), Form("Amplitudes for channel %d of third detector", k), 1200, 800);
        hists.hAmpChannel3[k]->Draw();
        cAmpChannel3 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector3.pdf", outDir.Data()));
    }


    cSavePDF_AmplitudesPerChannel -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector1.pdf]", outDir.Data()));
    cSavePDF_AmplitudesPerChannel2 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector2.pdf]", outDir.Data()));
    cSavePDF_AmplitudesPerChannel3 -> SaveAs(Form("%s/Amplitude_Distribution_Per_Channel_detector3.pdf]", outDir.Data()));


    cout<<"till here fine"<<endl;

    //cell0time stamp vs amplitude per channel DA AGGIUNGERE

}





/*vector<vector<int>> ReadFile() {

    ifstream channel_file("channel_mapping.txt");
    if (!channel_file.is_open()) {
        cerr << "Error: could not open channel_mapping.txt" << endl;
        return {};
    }

    channel_file.ignore(1000, '\n');

    int channel, board;
    double x, y;

    vector<vector<int>> coordinates;

    while (channel_file >> detector >> FEB >> channel >> x >> y ) {
        coordinates.push_back({
            detector,
            FEB
            channel,
            static_cast<int>(x),
            static_cast<int>(y),
            
        });
    }

    cout << "\n\n\n        ciao!" << endl;

    return coordinates;
}*/



/*per processe events november:
1)ricordarsi che su feb0 channel 0 ci sta l'mcp mentre i segnali stanno su feb1 e feb3
2)abbiamo solo un detector 
3) run 209 solo feb 1
4) riscrivere il txt file delle coordinate
5) */


/*void ProcessEvents_November(TTree *tree, TreeBranches &b, Histograms &h, vector<vector<int>> &coordinates) ///aggiungere vettore di vettori
{   
    Long64_t nentries = tree->GetEntries();
    cout << "Total number of entries: " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++) 
    {
        tree->GetEntry(i);
        
        // Basic histograms for amplitude
        for (int j = 0; j < b.ArraySize; j++) 
        {
            if (b.Board[j] == 0)
            {
                if(b.Channel[j]==0){
                h.hAmpAll1->Fill(b.Amplitude[j]);
                h.hBaseline1->Fill(b.Baseline[j]);
                h.hChannelHits1->Fill(b.Channel[j]);
                h.hAmpChannel1[0]->Fill(b.Amplitude[j]);
                }
                }
            }
            else if (b.Board[j] == 1)
            {
                h.hAmpAll2->Fill(b.Amplitude[j]);
                h.hBaseline2->Fill(b.Baseline[j]);
                h.hChannelHits2->Fill(b.Channel[j]);
                for(int k=0; k<64; k++){
                    if(b.Channel[j]==k){
                        h.hAmpChannel2[k]->Fill(b.Amplitude[j]);
                    }
                }
            }
            else if (b.Board[j] == 3)
            {
                h.hAmpAll3->Fill(b.Amplitude[j]);
                h.hBaseline3->Fill(b.Baseline[j]);
                h.hChannelHits3->Fill(b.Channel[j]);
                for(int k=0; k<64; k++){
                    if(b.Channel[j]==k){
                        h.hAmpChannel3[k]->Fill(b.Amplitude[j]);
                    }
                }
            }
            

        }
            h.hAmpAll1->GetXaxis()->SetRange(
            h.hAmpAll1->FindFirstBinAbove(0),
            h.hAmpAll1->FindLastBinAbove(0)
            );
                    h.hAmpAll2->GetXaxis()->SetRange(
            h.hAmpAll2->FindFirstBinAbove(0),
            h.hAmpAll2->FindLastBinAbove(0)
            );
                    h.hAmpAll3->GetXaxis()->SetRange(
            h.hAmpAll3->FindFirstBinAbove(0),
            h.hAmpAll3->FindLastBinAbove(0)
            );
        //cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
        // Channels 
        //if (b.HitFeb[0] >= 1 && b.HitFeb[1] >= 1 && b.HitFeb[2] >= 1) //added the parto of hitfeb[2]
        //if(b.HitFeb[2]==1)//check con ultima analisis con trigger sul terzo det
        if(true)
        {
            double amp_FEB[3] = {-1.0, -1.0, -1.0};
            double base[3] = {-9999.0, -9999.0, -9999.0};
            

            for (int j = 0; j < b.ArraySize; j++) 
            {
                int det = b.Board[j];
                //cout << "Board: " << det << endl;

                if (det < 0 || det > 2) continue;
               
                // Just testing low amplitude cut
                //if (b.pulses_amplitude[j] < 0.05) continue;
                
                // amp_FEB[det] = b.pulses_amplitude[j];
                base[det] = b.Baseline[j];
                // peak_time[det] = b.pulses_peak_time[j];

                for(int k=0; k<64; k++){
                    if(b.Channel[j]==k){
                        if(det==0){
                            h.mapDet1->Fill(coordinates[k][1], coordinates[k][2]);
                        }
                        else if(det==1){
                            h.mapDet2->Fill(coordinates[k][1], coordinates[k][2]);
                        }
                        else if(det==3){
                            h.mapDet3->Fill(coordinates[k][1], coordinates[k][2]);
                        }
                    }
                }

                //cout << "Channel: " << b.Channel[j] << endl;

                
                if (det == 0)
                    h.mapDet1->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);
                if (det == 1)
                    h.mapDet2->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);
                if (det == 3)
                    h.mapDet3->Fill(b.pulses_channel_x[j], b.pulses_channel_y[j]);   for third detectors
            }
            
            h.hBaseline1->Fill(base[0]);
            h.hBaseline2->Fill(base[1]);
            h.hBaseline3->Fill(base[2]); //for the third det
            
        }
        
        if (i % 100000 == 0)
            cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }*/





    