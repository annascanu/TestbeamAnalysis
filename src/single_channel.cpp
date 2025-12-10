/*
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
                Analysis for the October 2025 ENUBET testbeam picosec time resolution.
                               Authors: A. Scanu, R. Speziali
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

Goal: obtain time resolution of picosec detectors.
To compile: c++ analysis_waveforms.cpp `root-config --cflags --libs` -o analysis.out
Roadmap:
- [x] work on amplitude plotting
- [] understand what threshold to use for t_0 determination (ongoing)
- [] understand how to deal with charge sharing (not started)
- [] ...

Things we are trying to understand or need to ask:
- Q: Order of pulse_time and pulse_amplitude inside the array (probably understood: should use Board to discriminate?)
  A: Tentative answer is to try and use the Board array to identify which pulse corresponds to which FEB. What I (Anna) don't understand is if
- Q: how to take into account the RMS baseline value for the CDF threshold?
  A: (still in progress)
*/

/*
IMPORTANT NOTE ABOUT FILE PATHS:
Processed files can be found at this path:
/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/
- Muon runs: run 19, 20, 21, 22
- 15 GeV pion runs: run 23, ...
*/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <TStyle.h>
#include <TLatex.h>
#include <TF1.h>
#include "TMath.h"
#include <TH2F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include "TProfile.h"
#include "TGraph2D.h"
#include <vector>
using namespace std;

double media(const vector<double>& v) {
    if (v.empty()) return 0; // evitare divisione per zero
    double somma = 0;
    for (double x : v) {
        somma += x;
    }
    return somma / v.size();
}

// AS: What is this?
/*
TGraph* gGraph = 0;
Double_t myLine(Double_t *x, Double_t *par)
{
    double xx = x[0];

    // Trovo l’indice del punto più vicino
    int n = gGraph->GetN();
    double gx, gy;

    // ROOT passa i punti in ordine → cerco il punto con la x corrispondente
    for (int i = 0; i < n; i++)
    {
        gGraph->GetPoint(i, gx, gy);
        if (fabs(gx - xx) < 1e-9)
        {   // punto giusto
            // Limiti su y
            double ymin = 0;       // <-- METTI TU
            double ymax = 30;     // <-- METTI TU

            if (gy < ymin || gy > ymax)
            {
                TF1::RejectPoint();   // scarta il punto
                return 0;
            }
            break;
        }
    }

    // Modello lineare
    return par[0] + par[1] * xx;
}*/



int main()
{
    // gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(0);

    // -------------------------------------------------
    //     Load run file with processed waveforms
    // -------------------------------------------------

    // TString filename = "/eos/experiment/neutplatform/enubet/testbeam2025/picosec_data/sampic_runs/rootSampicData/processed_waveforms/sampic_run22_final.root"; // Filename while running on lxplus
    //TString filename = "/Users/anna/Developing/PhD/Testbeam2025/sampic_run22_final.root"; // Filename while running on Anna's machine
     TString filename = "/home/riccardo-speziali/Scrivania/October_2025/root_tree/sampic_run19_final.root"; // Filename while running on Riccardo's machine
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


    // -----------------------
    //    Set up branches
    // -----------------------
    const int MAXPULSES = 200; // Maybe could be less.
    Int_t npulses;

    Double_t pulses_amplitude[MAXPULSES], pulses_integral[MAXPULSES], pulses_time_cfd10[MAXPULSES], pulses_time_cfd20[MAXPULSES], pulses_time_cfd30[MAXPULSES], integral[MAXPULSES], pulses_time_cfd50[MAXPULSES], pulses_time_cfd60[MAXPULSES], Baseline[MAXPULSES], pulses_peak_time[MAXPULSES],  pulses_rise_time[MAXPULSES], mean_1[MAXPULSES], mean_2[MAXPULSES], conteggi_canali1[64], conteggi_canali2[64] ;
    Bool_t  pulses_bad_pulse[MAXPULSES];
    Int_t HitFeb[3], Board[MAXPULSES], Channel[MAXPULSES];
    int tot_hit_feb;
    Double_t time_12, time_23, time_13, time20_12, time20_13, time20_23;
    Double_t time30_12, time30_13, time30_23, time50_12, peso, peso2;
    vector<vector<double>> tabella1(64);
     vector<vector<double>> tabella2(64);

    tree->SetBranchAddress("npulses", &npulses);
    tree->SetBranchAddress("pulses_amplitude", &pulses_amplitude);
    tree->SetBranchAddress("pulses_integral", &pulses_integral);
    tree->SetBranchAddress("pulses_time_cfd10", &pulses_time_cfd10);
    tree->SetBranchAddress("pulses_time_cfd20", &pulses_time_cfd20);
    tree->SetBranchAddress("pulses_time_cfd30", &pulses_time_cfd30);
    tree->SetBranchAddress("pulses_bad_pulse", &pulses_bad_pulse);
    tree->SetBranchAddress("HitFeb", &HitFeb);
    tree->SetBranchAddress("Board", &Board);
    tree->SetBranchAddress("pulses_time_cfd50", &pulses_time_cfd50);
    tree->SetBranchAddress("pulses_time_cfd60", &pulses_time_cfd60);
    tree->SetBranchAddress("pulses_baseline", &Baseline);
    tree->SetBranchAddress("Channel", &Channel);
    tree->SetBranchAddress("pulses_peak_time", &pulses_peak_time);
    tree->SetBranchAddress("pulses_rise_time", &pulses_rise_time);




    // ------------------------------
    //     Initialize histograms
    // ------------------------------
    TH1F *hAmpAll  = new TH1F("hAmpAll", "All pulse amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hAmp1Hit = new TH1F("hAmp1Hit", "Single-Hit event amplitudes;Amplitude [a.u.];Counts", 1000, 0, 5);
    TH1F *hQ       = new TH1F("hQ", "Pulse integral (namely charge);Integral [a.u.];Counts", 1000, 0, 10);
    TH1F *hTime20   = new TH1F("hTime20", "Constant fraction discriminator set at 20% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime30    = new TH1F("hTime30", "Constant fraction discriminator set at 30% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime50    = new TH1F("hTime50", "Constant fraction discriminator set at 50% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime60    = new TH1F("hTime60", "Constant fraction discriminator set at 60% Time;Time [ps];Counts", 1000, 0, 10000);
    TH1F *hTime60_2   = new TH1F("hTime60_2", "Constant fraction discriminator set at 60% Time;Time [ps];Counts", 1000, 0, 10000);

    // Histograms for triple hits
    TH1F *htriple1  = new TH1F("htriple1", "amplitude of a triple", 100, 0, 0.5);
    TH1F *htriple2  = new TH1F("htriple2", "amplitude of a triple", 100, 0, 0.5);
    TH1F *htriple3  = new TH1F("htriple3", "amplitude of a triple", 100, 0, 0.5);

    // Histograms for time differences for the different values of CDF
    TH1F *htime_12    = new TH1F("htime_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime_13    = new TH1F("htime_13", "time diff ;Time [ps];Counts", 100, -0.2, 0.2);
    TH1F *htime_23    = new TH1F("htime_23", "time diff ;Time [ps];Counts", 100, -0.2, 0.2);

    TH1F *htime20_12    = new TH1F("htime20_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);//unico gragfico che sto considereando
    TH1F *htime20_13    = new TH1F("htime20_13", "time diff ;Time [ps];Counts", 100, -100, 100);
    TH1F *htime20_23    = new TH1F("htime20_23", "time diff ;Time [ps];Counts", 100, -0.2, 0.2);

    TH1F *htime30_12    = new TH1F("htime30_12", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime30_13    = new TH1F("htime30_13", "time diff ;Time [ps];Counts", 100, -2000, 2000);
    TH1F *htime30_23    = new TH1F("htime30_23", "time diff ;Time [ps];Counts", 100, -2000, 2000);

    // Baseline histogram
    TH1F *hcharge1  = new TH1F("hcharge1", "charge1", 1000, 0, 10);
    TH1F *hcharge2  = new TH1F("hcharge2", "charge2", 1000, 0, 10);
    TH1F *hrisetime1  = new TH1F("hrisetime1", "risetime1", 1000, -2000, 2000);
    TH1F *hrisetime2  = new TH1F("hrisetime2", "risetime2", 1000, -2000, 2000);
    TH1F *hbaseline1  = new TH1F("hbaseline", "baseline", 1000, 0.9, 1.1);

    TH1F *hbaseline2  = new TH1F("hbaseline2", "baseline2", 1000, 0.9, 1.1);
    TH1F *hpeaktime1  = new TH1F("hpeaktime1", "peaktime1", 1000, -20000, 20000);
    TH1F *hpeaktime2  = new TH1F("hpeaktime2", "peaktime2", 1000, -20000, 20000);

    TH1F *conteggixcanale1  = new TH1F("conteggixcanale1", "conteggixcanale1", 1000, 0, 64);
    TH1F *conteggixcanale2  = new TH1F("conteggixcanale2", "conteggixcanale2", 1000, 0, 64);
    TH1F *conteggixcanale3  = new TH1F("conteggixcanale3", "conteggixcanale3", 100, 0, 7);



    ///Amplitude vs integral histogram
    TGraph *hintegral1 = new TGraph();
    TGraph *hintegral2 = new TGraph();
    TGraph *hintegral3 = new TGraph();
    //TGraph2D *timevsamplitude =new TGraph2D(); no more 2d plot
    TGraph *timevsintegral =new TGraph();
    //TGraph2D *timevsamplitude2 =new TGraph2D();no more 2d plot
    TGraph *timevsintegral2 =new TGraph();
    TGraph *timevsamplitude3 = new TGraph();
     TGraph *timevsamplitude =new TGraph();
      TGraph *timevsamplitude2 =new TGraph();

    TProfile *prof = new TProfile("prof", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);
    TProfile *prof2 = new TProfile("prof2", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);
    TProfile *prof3 = new TProfile("prof3", "Tempo medio per bin di ampiezza;Ampiezza;Tempo [ps]", 50, 0, 1);

    TGraph *cfdvschannel1 = new TGraph();

    // -----------------------
    //        Read tree
    // -----------------------
    Long64_t nentries = tree->GetEntries();
    bool cut[nentries]={false};
    int nGood, idx;

    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        nGood = 0;
        idx = -1;

        for (int j = 0; j < npulses; j++) // Loop inside the event
        {
            if (pulses_bad_pulse[j])
                continue; // Skip events where the fit didn't work

            nGood++;
            idx = j;
            hAmpAll->Fill(pulses_amplitude[j]);
        }

        if (nGood == 1 && idx >= 0)
        {
            hAmp1Hit->Fill(pulses_amplitude[idx]);
            hQ->Fill(pulses_integral[idx]);
            //hTime->Fill(pulses_time_cfd20[idx]);
        }

        // ---- Triple hits for all three FEBs ----
        // AS: Replaced assumption about pulse index ordering with explicit Board[] mapping
        // also see the hitfeb = [1, 1, anything] (otherwise we lose statistics because of smaller geometrical acceptance of 3rd detector)
        //i want hitfeb[2]==1
       // if (HitFeb[0] == 1 && HitFeb[1] == 1 && nGood == 2)
        if (HitFeb[0] == 1 && HitFeb[1] == 1 && nGood == 3 && HitFeb[2]==1)
        //if (HitFeb[0] == 1 && HitFeb[1] == 1 && nGood >= 2) //condizione che mi garantisce hit sulle prime due e whatever sull'ultim detector'
        //if (HitFeb[0] == 1 && HitFeb[1] == 1 && HitFeb[2] == 1 && nGood == 3) //triple per check
        {
            // Arrays indexed by FEB number (0,1,2)
            double ampFEB[3]  = {-1.0, -1.0, -1.0};
            double cfd10[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd20[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd30[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd50[3]   = {-9999.0, -9999.0, -9999.0};
            double cfd60[3]   = {-9999.0, -9999.0, -9999.0};
            double base[3]    = {-9999.0, -9999.0, -9999.0};
            double peak_time[3]  = {-9999.0, -9999.0, -9999.0};
            double risetime[3] = {-9999.0, -9999.0, -9999.0};

            // Fill according to Board[]
            for (int j = 0; j < npulses; j++)
            {
                if (pulses_bad_pulse[j])
                    continue;
                int det = Board[j]; // Which detector this pulse belongs to (either 0, 1, 2)

                if (det < 0 || det > 2)
                {
                    cout << "WARNING: unexpected Board index = " << det << " (event " << i << ", pulse " << j << ")" << endl;
                    continue;
                }

                //if(pulses_amplitude[j]>0 && pulses_integral[j]/pulses_amplitude[j] < 2000 )
                //{
                ampFEB[det] = pulses_amplitude[j];
                cfd10[det]  = pulses_time_cfd10[j];
                cfd20[det]  = pulses_time_cfd20[j];
                cfd30[det]  = pulses_time_cfd30[j];
                integral[det]= pulses_integral[j];
                cfd50[det]=pulses_time_cfd50[j];
                cfd60[det]=pulses_time_cfd60[j];
                base[det]=Baseline[j];
                peak_time[det]= pulses_peak_time[j];
                risetime[det]= pulses_rise_time[j];
                hTime20->Fill(cfd20[0]);// questi tre prima erano riempiti con j non corretti

                hTime50->Fill(cfd50[0]);//
                hTime60->Fill(cfd60[0]);
                hTime60_2->Fill(cfd60[1]);
                hcharge1->Fill(integral[0]);
                hcharge2->Fill(integral[1]);
                hrisetime1->Fill(risetime[0]);
                hrisetime2->Fill(risetime[1]);



                //cut[i]=true;}
            }

            // Skip se un detector non ha dato segnale buono
            if (ampFEB[0] <= 0 || ampFEB[1] <= 0 || ampFEB[2] <=0)
                continue;
            //&& (Channel[0]<16 || Channel[0]>32)
            //&& cfd60[0]<3600 &&
            //&& cfd30[0]<3100 && cfd30[0]>2900 cut on cfd 30%  && integral[0]/ampFEB[0]>8 && integral[1]/ampFEB[1]>8    ///// && ampFEB[0]> 0.2 && cfd60[0]<3350 && cfd60[0]>3000ss
           //if(Channel[1]<16 || Channel[1]>23) nonnfunziona: && (ampFEB[0]<0.03 || ampFEB[0]>0.03) &&(ampFEB[1]<0.03 || ampFEB[1]>0.0





           if(Channel[0]==12 && Channel[1]==12 && Channel[2]==5){
            if(peak_time[1]<0)
                continue;

                // Fill amplitude histograms           // AS: why this x 1000 scaling? -> To get values in mV
                htriple1->Fill(ampFEB[0]);
                htriple2->Fill(ampFEB[1]);
                htriple3->Fill(ampFEB[2]);

                if(Channel[0]<8)
                    peso=0;
                if(Channel[0]>7 && Channel[0]<16)
                    peso=1;
                if(Channel[0]>15 && Channel[0]<24)
                    peso=2;
                if(Channel[0]>23 && Channel[0]<32)
                    peso=3;
                if(Channel[0]>31 && Channel[0]<40)
                    peso=4;
                if(Channel[0]>39 && Channel[0]<48)
                    peso=5;
                if(Channel[0]>47 && Channel[0]<56)
                    peso=6;
                if(Channel[0]>55)
                    peso=7;

                if(Channel[1]<8)
                    peso2=0;
                if(Channel[1]>7 && Channel[1]<16)
                    peso2=1;
                if(Channel[1]>15 && Channel[1]<24)
                    peso2=2;
                if(Channel[1]>23 && Channel[1]<32)
                    peso2=3;
                if(Channel[1]>31 && Channel[1]<40)
                    peso2=4;
                if(Channel[1]>39 && Channel[1]<48)
                    peso2=5;
                if(Channel[1]>47 && Channel[1]<56)
                    peso2=6;
                if(Channel[1]>55)
                    peso2=7;

            //mean for every channel


                for(int i=0; i<64; i++){
                    if(Channel[0]==i){
                        tabella1[i].push_back(peak_time[0]);
                       // conteggi_canali1[i] += 1;

                    }
                }


                for(int i=0; i<64; i++){
                    if(Channel[1]==i){
                        tabella2[i].push_back(peak_time[1]);
                        //conteggi_canali2[i] += 1;

                    }
                }




hTime30->Fill(cfd30[0]);

                // CFD = 10%
                //time_12 = cfd60[0] - cfd60[1];
                //time_23 = cfd10[1] - cfd10[2];
                //time_13 = cfd10[0] - cfd10[2];

                //htime_12->Fill(time_12);
                //htime_23->Fill(time_23);
                //htime_13->Fill(time_13);

                // CFD = 20%
                //time20_12 = cfd20[0] - cfd20[1];
                //if(fabs(time20_12)>)continue;
                //time20_23 = cfd20[1] - cfd20[2];
                //time20_13 = cfd20[0] - cfd20[2];

                //htime20_12->Fill(time20_12);
                //htime20_23->Fill(time20_23);
                //htime20_13->Fill(time20_13);
                hpeaktime1->Fill(peak_time[0]);
                hpeaktime2->Fill(peak_time[1]);


                // CFD = 60%
                //double timedif1= peak_time[0]-cfd60[0];
                //double timedif2= peak_time[1]-cfd60[1];
                time30_12 = cfd30[0]-cfd30[1];
                time30_23 = cfd30[1] - cfd30[2];
                time30_13 = cfd30[0] - cfd30[2];



                htime30_12->Fill(time30_12);
                htime30_23->Fill(time30_23);
                htime30_13->Fill(time30_13);

                hbaseline1->Fill(base[0]);
                hbaseline2->Fill(base[1]);



                ///Integral vs amplitude distribution
                timevsamplitude->SetPoint(timevsamplitude->GetN(), ampFEB[0], cfd30[0]);
                timevsamplitude2->SetPoint(timevsamplitude2->GetN(), ampFEB[1], cfd30[1]);
                timevsamplitude3->SetPoint(timevsamplitude3->GetN(), ampFEB[2], cfd30[2]);

                prof->Fill(ampFEB[0], cfd60[0]);
                prof2->Fill(ampFEB[1], cfd60[1]);
                prof3->Fill(ampFEB[2], cfd60[2]);
                timevsintegral->SetPoint(timevsintegral->GetN(), integral[0], cfd60[0]);
                //if(integral[0]<40)
                    hintegral1->SetPoint(hintegral1->GetN(),ampFEB[0],integral[0]);
                //if(integral[1]<26)
                    hintegral2->SetPoint(hintegral2->GetN(),ampFEB[1],integral[1]);
                    hintegral3->SetPoint(hintegral3->GetN(),ampFEB[2],integral[2]);



                    conteggixcanale1->Fill(Channel[0]);
                    conteggixcanale2->Fill(Channel[1]);
                    conteggixcanale3->Fill(Channel[2]);














                if(ampFEB[0]<0)
                    cout<< "Amplitude < 0; event number: " << i << endl;
                cfdvschannel1->SetPoint(cfdvschannel1->GetN(), Channel[0], cfd60[0]);



        }
        }

        // Progress indicator for sanity :)
        if (i % 100000 == 0) cout << "Processed " << i << " / " << nentries << " events...\r" << flush;
    }




       for(int i=0; i<64; i++){
                        mean_1[i]=media(tabella1[i]);
                        //cout << "media canali primo riv:  "<< mean_1[i] <<endl;

                }



            for(int i=0; i<64; i++){
                        mean_2[i]=media(tabella2[i]);
                                                //cout << "media canali secondo riv:  "<< mean_2[i] <<endl;

                }
                  //for the new tree
    //TTree *copy = tree->CloneTree();
    TTree *tout = new TTree("channelMeans", "Medie per canale");
    tout->Branch("mean_1", mean_1, "mean_1[64]/D");
    tout->Branch("mean_2", mean_2, "mean_2[64]/D");


        tout->Fill();
    cout << "Until here code works fine" << endl;

    /// --- TF1 Polya function ---
    TF1 *fpolyaAmpl1 = new TF1("fpolyaAmpl1", "([0]/[1])*((([2]+1)^([2]+1)*(x/[1])^[2])/(TMath::Gamma([2]+1)))*exp(-([2]+1)*x/[1])", 0, 800);
    fpolyaAmpl1->SetTitle("Polya fit (ampl.)");
    fpolyaAmpl1->SetParName(0,"c");
    fpolyaAmpl1->SetParName(1,"#bar{Q}");
    fpolyaAmpl1->SetParName(2,"#theta");
    fpolyaAmpl1->SetLineColor(kGreen+2);
    //fpolyaAmpl1->SetRange(minax, amplMax[ci]);
    fpolyaAmpl1->SetNpx(10000);
    //fpolyaAmpl1->SetParameter(0,800);
    fpolyaAmpl1->SetParameter(1,180);
    //fpolyaAmpl1->SetParameter(2,30);
    //fpolyaAmpl1->SetParLimits(2,0.001,5000);
    fpolyaAmpl1->SetRange(0, 800);

    TF1 *fpolyaAmpl2 = (TF1*)fpolyaAmpl1->Clone("fpolyaAmpl2");
    TF1 *fpolyaAmpl3 = (TF1*)fpolyaAmpl1->Clone("fpolyaAmpl3");
    TF1 *fpolyacharge1 = (TF1*)fpolyaAmpl1->Clone("fpolyacharge1");
    TF1 *fpolyacharge2 = (TF1*)fpolyaAmpl1->Clone("fpolyacharge2");


    //fpolyaAmpl2->SetParameter(2,1);

    // --- Tre canvas separati ---
    // Canvas 1
    // AS: removed grids, I think the plots look nicer, but can be re-added if needed
    TCanvas *c1 = new TCanvas("c1", "htriple1", 900, 700);
    htriple1->SetLineColor(kAzure+1);
    htriple1->SetLineWidth(3);
    htriple1->SetFillColorAlpha(kAzure-4, 0.35);
    htriple1->SetTitle("Amplitudes of waveforms on the first FEB");
    htriple1->Draw("HIST");
    htriple1->Fit(fpolyaAmpl1,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
    c1->Update();
    c1->SaveAs("hist1.pdf");

    // Canvas 2
    TCanvas *c2 = new TCanvas("c2", "htriple2", 900, 700);
    htriple2->SetLineColor(kGreen+1);
    htriple2->SetLineWidth(3);
    htriple2->SetFillColorAlpha(kGreen-4, 0.35);
    htriple2->SetTitle("Amplitudes of waveforms on the second FEB");
    htriple2->Draw("HIST");
    htriple2->Fit(fpolyaAmpl2,"R"); // Fit su htriple2
    fpolyaAmpl2->Draw("SAME");
    c2->Update();
    c2->SaveAs("hist2.pdf");

    // Canvas 3
    TCanvas *c3 = new TCanvas("c3", "htriple3", 800, 600);
    htriple3->SetLineColor(kRed+1);
    htriple3->SetLineWidth(3);
    htriple3->SetFillColorAlpha(kRed-4, 0.35);
    htriple3->SetTitle("Amplitudes of waveforms on the third FEB");
    htriple3->Draw("HIST");
    htriple3->Fit(fpolyaAmpl3,"R"); // Fit su htriple3
    //fpolyaAmpl3->Draw("SAME");
    c3->Update();
    c3->SaveAs("hist3.pdf");


    TCanvas *ccharge1 = new TCanvas("ccharge1", "charge1", 900, 700);
    hcharge1->SetLineColor(kAzure+1);
    hcharge1->SetLineWidth(3);
    hcharge1->SetFillColorAlpha(kAzure-4, 0.35);
    hcharge1->SetTitle("Charge of waveforms on the first FEB");
    hcharge1->Draw("HIST");
    hcharge1->Fit(fpolyacharge1,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
    ccharge1->Update();



TCanvas *ccharge2 = new TCanvas("ccharge2", "charge2", 900, 700);
    hcharge2->SetLineColor(kAzure+1);
    hcharge2->SetLineWidth(3);
    hcharge2->SetFillColorAlpha(kAzure-4, 0.35);
    hcharge2->SetTitle("Charge of waveforms on the second FEB");
    hcharge2->Draw("HIST");
    hcharge2->Fit(fpolyacharge2,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
    ccharge2->Update();



TCanvas *crisetime1 = new TCanvas("crisetime1", "risetime1", 900, 700);
    hrisetime1->SetLineColor(kAzure+1);
    hrisetime1->SetLineWidth(3);
    hrisetime1->SetFillColorAlpha(kAzure-4, 0.35);
    hrisetime1->SetTitle("risetime of waveforms on the first FEB");
    hrisetime1->Draw("HIST");
    crisetime1->Update();




TCanvas *crisetime2 = new TCanvas("crisetime2", "risetime2", 900, 700);
    hrisetime2->SetLineColor(kAzure+1);
    hrisetime2->SetLineWidth(3);
    hrisetime2->SetFillColorAlpha(kAzure-4, 0.35);
    hrisetime2->SetTitle("risetime of waveforms on the second FEB");
    hrisetime2->Draw("HIST");
    crisetime2->Update();


TCanvas *cpeaketime1 = new TCanvas("cpeaktime1", "peaktime1", 900, 700);
    hpeaktime1->SetLineColor(kAzure+1);
    hpeaktime1->SetLineWidth(3);
    hpeaktime1->SetFillColorAlpha(kAzure-4, 0.35);
    hpeaktime1->SetTitle("peaktime of waveforms on the first FEB");
    hpeaktime1->Draw("HIST");
    cpeaketime1->Update();



TCanvas *cpeaketime2 = new TCanvas("cpeaktime2", "peaktime2", 900, 700);
    hpeaktime2->SetLineColor(kAzure+1);
    hpeaktime2->SetLineWidth(3);
    hpeaktime2->SetFillColorAlpha(kAzure-4, 0.35);
    hpeaktime2->SetTitle("peaktime of waveforms on the second FEB");
    hpeaktime2->Draw("HIST");
    cpeaketime2->Update();




















    // --- Canvas unico con tutti e tre ---
    TCanvas *cAll = new TCanvas("cAll", "All Histograms", 900, 700);

    htriple1->SetLineColor(kAzure+1);
    htriple1->SetLineWidth(3);
    htriple1->SetFillColorAlpha(kAzure-4, 0.35);
    htriple1->Draw("HIST");
    htriple1->Fit(fpolyaAmpl1,"R SAME");

    htriple2->SetLineColor(kGreen+1);
    htriple2->SetLineWidth(3);
    htriple2->SetFillColorAlpha(kGreen-4, 0.35);
    htriple2->Draw("HIST SAME");
    htriple2->Fit(fpolyaAmpl2,"R SAME");

    htriple3->SetLineColor(kRed+1);
    htriple3->SetLineWidth(3);
    htriple3->SetFillColorAlpha(kRed-4, 0.35);
    htriple3->Draw("HIST SAME");
    htriple3->Fit(fpolyaAmpl3,"R SAME");

    cAll->Update();
    cAll->SaveAs("histAll.pdf");

    ///amplitude vs integral plot
    cout<<"till here is fine"<<endl;

    TCanvas *cintegral = new TCanvas("cintegral", "amplitudevsintegral", 900, 700);

hintegral1->SetTitle("Amplitude vs Integral;amplitude(V);integral");
hintegral1->SetName("hintegral1");

hintegral1->SetMarkerStyle(20);
hintegral1->SetMarkerSize(1.0);
hintegral1->SetMarkerColor(kAzure+2);
hintegral1->SetLineColor(0);    // nessuna linea
hintegral1->SetLineWidth(0); // nessuna linea

double xmin = hintegral1->GetXaxis()->GetXmin();
double xmax = hintegral1->GetXaxis()->GetXmax();

// --- Configurazione della Funzione di Fit (TF1) ---
TF1* avsi = new TF1("avsi", "[0]+[1]*x", 0, 2);
//avsi->SetParameter(0, 0.0);    // intercetta iniziale
avsi->SetParLimits(1, 12, 26);   // pendenza iniziale stimata
avsi->SetNpx(5000);

// 1. ESEGUI IL FIT (Questa è la parte critica)
// L'opzione "Q" (quiet) può essere utile se il fit è rumoroso,
// ma per ora usiamo la chiamata standard.
hintegral1->Fit(avsi);

// 2. DISEGNA IL TGRAPH (disegna i punti e crea gli assi)
// Usa "AP" che è corretto per i punti e gli assi.
hintegral1->Draw("AP");

// 3. DISEGNA LA LINEA DI FIT SOPRA (aggiungi la linea di fit appena calcolata)
avsi->SetLineColor(kRed);
avsi->SetLineWidth(2);
// Disegna solo la linea di fit senza ricalcolare nulla
avsi->Draw("same");

// 4. Aggiorna il Canvas
cintegral->Update();
// cintegral->Draw(); // A volte serve anche questo, prova se Update non basta.


    cout << "till here is fine" << endl;

    TCanvas *cintegral2 = new TCanvas("cintegral2", "amplitudevsintegral", 900, 700);
    hintegral2->SetName("hintegral2");
    hintegral2->SetTitle("Amplitude vs Integral;amplitude;integral");
    // *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***
    hintegral2->SetMarkerStyle(20);
    hintegral2->SetMarkerSize(1.0);
    hintegral2->SetMarkerColor(kAzure+2);
    hintegral2->SetLineColor(0);    // nessuna linea
    hintegral2->SetLineWidth(0); // nessuna linea

    TF1* avsi2 = new TF1("avsi2", "[0] + [1]*x");
    avsi2->SetNpx(500);
    hintegral2->Fit(avsi2, "R");
    avsi2->SetLineColor(kRed);
    avsi2->SetLineWidth(2);

    // *** SOLO ORA DISEGNA ***
    hintegral2->Draw("AP");

    hintegral2->Fit(avsi, "R");  // "R" = fit solo nell'intervallo di x specificato
    hintegral2->Draw("AP");
    avsi2->Draw("same");          // adesso la linea del fit viene disegnata*/

    cintegral2->Update();






TCanvas *cconteggi = new TCanvas("cconteggi", "cconteggi", 900, 700);
    conteggixcanale1->SetLineColor(kAzure+1);
    conteggixcanale1->SetLineWidth(3);
    conteggixcanale1->SetFillColorAlpha(kAzure-4, 0.35);
    conteggixcanale1->SetTitle("eventi per canale");
   // htriple1->Draw("HIST");
    //htriple1->Fit(fpolyaAmpl1,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
    cconteggi->Update();
    //cconteggi->SaveAs("hist1.pdf");


TCanvas *cconteggi2 = new TCanvas("cconteggi2", "cconteggi2", 900, 700);
    conteggixcanale2->SetLineColor(kAzure+1);
    conteggixcanale2->SetLineWidth(3);
    conteggixcanale2->SetFillColorAlpha(kAzure-4, 0.35);
    conteggixcanale2->SetTitle("eventi per canale");
   // htriple1->Draw("HIST");
    //htriple1->Fit(fpolyaAmpl1,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
    cconteggi2->Update();
    //cconteggi->SaveAs("hist1.pdf");



TCanvas *cconteggi3 = new TCanvas("cconteggi3", "cconteggi3", 900, 700);
    conteggixcanale3->SetLineColor(kAzure+1);
    conteggixcanale3->SetLineWidth(3);
    conteggixcanale3->SetFillColorAlpha(kAzure-4, 0.35);
    conteggixcanale3->SetTitle("eventi per canale");
   // htriple1->Draw("HIST");
    //htriple1->Fit(fpolyaAmpl1,"R"); // Fit su htriple1
    //fpolyaAmpl1->Draw("SAME");
    cconteggi3->Update();
    //cconteggi->SaveAs("hist1.pdf");









    // ----------------------------------------------------------------------------------------------------

    // cout<<"time walk"<<endl;
    // TCanvas *ctimeamlitude = new TCanvas("ctimeaplitude", "amplitudevsintegral", 900, 700);
    // //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);
    // timevsamplitude->SetName("timeaplitude");
    // timevsamplitude->SetTitle("SAT vs amplitude; Amplitude; SAT[us]");
    //
    // // *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***
    // /*timevsamplitude->SetMarkerStyle(20);
    // timevsamplitude->SetMarkerSize(1.0);
    // timevsamplitude->SetMarkerColor(kAzure+2);
    //
    // timevsamplitude->SetLineColor(0);    // nessuna linea
    // timevsamplitude->SetLineWidth(0);    // nessuna linea*/
    //
    // gStyle->SetPalette(kViridis);
    // //timevsamplitude->Fit("timewalk", "R");
    // // *** SOLO ORA DISEGNA ***
    // timevsamplitude->Draw("COLZ");
    // ctimeamlitude->Update();
    //
    // TCanvas *ctimeamlitude2 = new TCanvas("ctimeaplitude2", "amplitudevsintegral", 900, 700);
    // //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);
    // timevsamplitude2->SetName("timeaplitude2");
    // timevsamplitude2->SetTitle("SAT vs amplitude; Amplitude; SAT[us]");
    //
    // // *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***
    // /*timevsamplitude->SetMarkerStyle(20);
    // timevsamplitude->SetMarkerSize(1.0);
    // timevsamplitude->SetMarkerColor(kAzure+2);
    // timevsamplitude->SetLineColor(0);    // nessuna linea
    // timevsamplitude->SetLineWidth(0);    // nessuna linea*/
    //
    // gStyle->SetPalette(kViridis);
    // //timevsamplitude->Fit("timewalk", "R");
    // // *** SOLO ORA DISEGNA ***
    // timevsamplitude2->Draw("COLZ");
    // ctimeamlitude2->Update();



TCanvas *ctimeamlitude = new TCanvas("ctimeaplitude", "amplitudevsintegral", 900, 700);



   //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);


    timevsamplitude->SetName("timeaplitude");
timevsamplitude->SetTitle("SAT vs amplitude; Amplitude; SAT");
// *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***

timevsamplitude->SetMarkerStyle(20);
timevsamplitude->SetMarkerSize(1.0);
timevsamplitude->SetMarkerColor(kAzure+2);
timevsamplitude->SetLineColor(0);    // nessuna linea
timevsamplitude->SetLineWidth(0);    // nessuna linea


//timevsamplitude->Fit("timewalk", "R");
// *** SOLO ORA DISEGNA ***
timevsamplitude->Draw("AP");
ctimeamlitude->Update();


TCanvas *ctimeamlitude2 = new TCanvas("ctimeaplitude2", "amplitudevsintegral", 900, 700);



   //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);


    timevsamplitude2->SetName("timeaplitude");
timevsamplitude2->SetTitle("SAT vs amplitude; Amplitude; SAT");
// *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***

timevsamplitude2->SetMarkerStyle(20);
timevsamplitude2->SetMarkerSize(1.0);
timevsamplitude2->SetMarkerColor(kAzure+2);
timevsamplitude2->SetLineColor(0);    // nessuna linea
timevsamplitude2->SetLineWidth(0);    // nessuna linea


//timevsamplitude->Fit("timewalk", "R");
// *** SOLO ORA DISEGNA ***
timevsamplitude2->Draw("AP");
ctimeamlitude2->Update();






TCanvas *cAmpVsTime = new TCanvas("cAmpVsTime", "Amplitude vs Time", 900, 700);
//cAmpVsTime->SetGrid();
//cAmpVsTime->SetTickx();
//cAmpVsTime->SetTicky();

  //TF1 *timewalk = new TF1("timewalk", "[0]+[1]/x", 0, 800);


    timevsamplitude3->SetName("timeaplitude");
timevsamplitude3->SetTitle("SAT vs amplitude; Amplitude; SAT");
// *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***

timevsamplitude3->SetMarkerStyle(20);
timevsamplitude3->SetMarkerSize(1.0);
timevsamplitude3->SetMarkerColor(kAzure+2);
timevsamplitude3->SetLineColor(0);    // nessuna linea
timevsamplitude3->SetLineWidth(0);    // nessuna linea


//timevsamplitude->Fit("timewalk", "R");
// *** SOLO ORA DISEGNA ***
timevsamplitude3->Draw("AP");
// Legenda
/*auto leg = new TLegend(0.15, 0.80, 0.45, 0.92);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(timevsamplitude3, "Measured points", "lp");
leg->Draw();*/

cAmpVsTime->Update();




    TCanvas *ctvsint = new TCanvas("ctvsint", "Time vs Integral", 900, 700);
    ctvsint ->SetGrid();

    timevsintegral->SetName("timevsintegral");

    // Stile linea/punti
    timevsintegral->SetLineColor(0);
    timevsintegral->SetLineWidth(0);
    timevsintegral->SetMarkerStyle(20);
    timevsintegral->SetMarkerSize(1.1);
    timevsintegral->SetMarkerColor(kBlue+1);

    // Titoli e assi
    timevsintegral->SetTitle("SAT vs integral;Integral;SAT [ps]");
    timevsintegral->GetXaxis()->SetTitleSize(0.045);
    timevsintegral->GetYaxis()->SetTitleSize(0.045);
    timevsintegral->GetXaxis()->SetLabelSize(0.04);
    timevsintegral->GetYaxis()->SetLabelSize(0.04);
    timevsintegral->Draw("AP");

    // Eventuale legenda
    /*auto leg = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg->AddEntry(timevsintegral, "Misure", "lp");
    leg->Draw();*/

    ctvsint->Update();

    TCanvas *cmediabin = new TCanvas("cmediabin", "Slewing Correction", 800, 600);
    //prof->Fit("pol3");
    TF1 *empirico = new TF1("empirico", "[0] + [1]/x", 0, 1);
    prof->Fit(empirico, "R");
    prof->Draw();
    cmediabin->Update();

    TCanvas *cmediabin2 = new TCanvas("cmediabin2", "Slewing Correction", 800, 600);
    //prof->Fit("pol3");
    TF1 *empirico2 = new TF1("empirico2", "[0] + [1]/x", 0, 1);
    prof2->Fit(empirico2, "R");
    prof2->Draw();
    cmediabin2->Update();

     TCanvas *cmediabin3 = new TCanvas("cmediabin2", "Slewing Correction", 800, 600);
    //prof->Fit("pol3");
    TF1 *empirico3 = new TF1("empirico3", "[0] + [1]/x", 0, 1);
    prof3->Fit(empirico3, "R");
    prof3->Draw();
    cmediabin3->Update();

    /*TCanvas *cintegral3 = new TCanvas("cintegral3", "amplitudevsintegral", 900, 700);
    hintegral3->SetName("hintegral1");
    hintegral3->SetTitle("Amplitude vs Integral;amplitude;integral");

    // *** IMPOSTAZIONE STILE PRIMA DEL DRAW ***
    hintegral3->SetMarkerStyle(20);
    hintegral3->SetMarkerSize(1.0);
    hintegral3->SetMarkerColor(kAzure+2);
    hintegral3->SetLineColor(0);    // nessuna linea
    hintegral3->SetLineWidth(0);    // nessuna linea
    // *** SOLO ORA DISEGNA ***
    hintegral3->Draw("AP");
    cintegral3->Update();*/

    //baseline
    TCanvas *cbaseline = new TCanvas("cbaseline", "baseline", 900, 700);
    //TCanvas *ctime12 = new TCanvas("ctime12", "htriple1", 900, 700);

    // 1. Definisci e Configura la Funzione Gaussiana (TF1)
    TF1 *f_gaus_b = new TF1("f_gaus", "gaus", 0.003, 0.009); // Definisci il nome, la formula, e il range iniziale del fit

    // 2. Imposta il Parametro della Media (Indice 1)
    // Il formato della gaus ROOT è: p0 * exp(-0.5 * ((x-p1)/p2)^2)
    // p0 = Ampiezza (Indice 0)
    // p1 = Media (Indice 1) <-- QUESTO E' QUELLO CHE VUOI IMPOSTARE
    // p2 = Sigma (Indice 2)
    //f_gaus->SetParameter(1, 300); // Esempio: imposta la media a 50.0
    //f_gaus->SetParameter(2, 80);
    // (Opzionale) Imposta i parametri iniziali per l'Ampiezza (0) e la Sigma (2) per un fit migliore
    // f_gaus->SetParameter(0, 500.0); // Ampiezza iniziale
    // f_gaus->SetParameter(2, 5.0);   // Sigma iniziale

    // 3. Esegui il Fit con la funzione predefinita
    // Nota: Passiamo l'oggetto TF1 al Fit, non solo il nome "gaus"
    hbaseline1->Fit(f_gaus_b, "R"); // "R" usa il range definito in TF1

    // 4. Configurazione e Disegno (il resto del tuo codice)
    hbaseline1->SetLineColor(kAzure+1);
    hbaseline1->SetLineWidth(3);
    hbaseline1->SetFillColorAlpha(kAzure-4, 0.35);
    hbaseline1->SetTitle("baseline fist det");
    hbaseline1->Draw("HIST");
    cbaseline->Update();

    // 5. Estrai i risultati dal TF1 dopo il fit (usa l'oggetto f_gaus)
    double amp_baseline   = f_gaus_b->GetParameter(0);
    double mean_baseline  = f_gaus_b->GetParameter(1); // Questo sarà il risultato del fit
    double sigma_baseline = f_gaus_b->GetParameter(2);





    //baseline
    TCanvas *cbaseline2 = new TCanvas("cbaseline2", "baseline2", 900, 700);
    //TCanvas *ctime12 = new TCanvas("ctime12", "htriple1", 900, 700);

    // 1. Definisci e Configura la Funzione Gaussiana (TF1)
    TF1 *f_gaus_b2 = new TF1("f_gaus", "gaus2", 0.003, 0.009); // Definisci il nome, la formula, e il range iniziale del fit

    // 2. Imposta il Parametro della Media (Indice 1)
    // Il formato della gaus ROOT è: p0 * exp(-0.5 * ((x-p1)/p2)^2)
    // p0 = Ampiezza (Indice 0)
    // p1 = Media (Indice 1) <-- QUESTO E' QUELLO CHE VUOI IMPOSTARE
    // p2 = Sigma (Indice 2)
    //f_gaus->SetParameter(1, 300); // Esempio: imposta la media a 50.0
    //f_gaus->SetParameter(2, 80);
    // (Opzionale) Imposta i parametri iniziali per l'Ampiezza (0) e la Sigma (2) per un fit migliore
    // f_gaus->SetParameter(0, 500.0); // Ampiezza iniziale
    // f_gaus->SetParameter(2, 5.0);   // Sigma iniziale

    // 3. Esegui il Fit con la funzione predefinita
    // Nota: Passiamo l'oggetto TF1 al Fit, non solo il nome "gaus"
    hbaseline1->Fit(f_gaus_b2, "R"); // "R" usa il range definito in TF1

    // 4. Configurazione e Disegno (il resto del tuo codice)
    hbaseline2->SetLineColor(kAzure+1);
    hbaseline2->SetLineWidth(3);
    hbaseline2->SetFillColorAlpha(kAzure-4, 0.35);
    hbaseline2->SetTitle("baseline second det");
    hbaseline2->Draw("HIST");
    cbaseline2->Update();

    // 5. Estrai i risultati dal TF1 dopo il fit (usa l'oggetto f_gaus)
    double amp_baseline2   = f_gaus_b2->GetParameter(0);
    double mean_baseline2 = f_gaus_b2->GetParameter(1); // Questo sarà il risultato del fit
    double sigma_baseline2 = f_gaus_b2->GetParameter(2);


    // ... estrazione degli errori
    /* TCanvas *ctime23 = new TCanvas("ctime23", "htriple1", 900, 700);
    htime20_23->SetLineColor(kAzure+1);
    htime20_23->SetLineWidth(3);
    htime20_23->SetFillColorAlpha(kAzure-4, 0.35);
    htime20_23->SetTitle("Amplitude");
    htime20_23->Draw("HIST");
    htime20_23->Fit("gaus"); // Fit su htime20_23
    ctime23->Update();

    TF1 *f2 = htime20_23->GetFunction("gaus");
    double amp_23   = f2->GetParameter(0);  ///taking gaus parameters
    double mean_23  = f2->GetParameter(1);
    double sigma_23 = f2->GetParameter(2);

    double e_amp_23   = f2->GetParError(0);
    double e_mean_23  = f2->GetParError(1);
    double e_sigma_23 = f2->GetParError(2);

    TCanvas *ctime13 = new TCanvas("ctime13", "htriple1", 900, 700);
    htime20_13->SetLineColor(kAzure+1);
    htime20_13->SetLineWidth(3);
    htime20_13->SetFillColorAlpha(kAzure-4, 0.35);
    htime20_13->SetTitle("Amplitude");
    htime20_13->Draw("HIST");
    htime20_13->Fit("gaus"); // Fit su htime20_13
    ctime13->Update();

    // Take gaussian from the correct histogram htime20_13
    TF1 *f3 = htime20_13->GetFunction("gaus");
    double amp_13   = f3->GetParameter(0);  // taking gaus parameters
    double mean_13  = f3->GetParameter(1);
    double sigma_13 = f3->GetParameter(2);

    double e_amp_13   = f3->GetParError(0);
    double e_mean_13  = f3->GetParError(1);
    double e_sigma_13 = f3->GetParError(2);*/

    TCanvas *ctimeperchannel = new TCanvas("ctimeperchannel", "Time vs Integral", 900, 700);
    ctvsint ->SetGrid();

    cfdvschannel1->SetName("cfdvschannel1");
    cfdvschannel1->SetLineColor(0);
    cfdvschannel1->SetLineWidth(0);
    cfdvschannel1->SetMarkerStyle(20);
    cfdvschannel1->SetMarkerSize(1.1);
    cfdvschannel1->SetMarkerColor(kBlue+1);
    cfdvschannel1->SetTitle("channel vs SAT;Channel;SAT [ps]");
    cfdvschannel1->GetXaxis()->SetTitleSize(0.045);
    cfdvschannel1->GetYaxis()->SetTitleSize(0.045);
    cfdvschannel1->GetXaxis()->SetLabelSize(0.04);
    cfdvschannel1->GetYaxis()->SetLabelSize(0.04);
    cfdvschannel1->Draw("AP");

    // Eventuale legenda
    /*auto leg = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg->AddEntry(timevsintegral, "Misure", "lp");
    leg->Draw();*/
    ctimeperchannel->Update();

    // ------------ Gaussian plot of the time difference ---------------

    TCanvas *ctime12 = new TCanvas("ctime12", "htriple1", 900, 700);
   //TCanvas *ctime12 = new TCanvas("ctime12", "htriple1", 900, 700);

    // 1. Definisci e Configura la Funzione Gaussiana (TF1)
    TF1 *f_gaus = new TF1("f_gaus", "gaus", -1000, 1000); // Definisci il nome, la formula, e il range iniziale del fit

    // 2. Imposta il Parametro della Media (Indice 1)
    // Il formato della gaus ROOT è: p0 * exp(-0.5 * ((x-p1)/p2)^2)
    // p0 = Ampiezza (Indice 0)
    // p1 = Media (Indice 1) <-- QUESTO E' QUELLO CHE VUOI IMPOSTARE
    // p2 = Sigma (Indice 2)
    //f_gaus->SetParameter(1, 300); // Esempio: imposta la media a 50.0
    //f_gaus->SetParameter(2, 80);
    // (Opzionale) Imposta i parametri iniziali per l'Ampiezza (0) e la Sigma (2) per un fit migliore
    // f_gaus->SetParameter(0, 500.0); // Ampiezza iniziale
    // f_gaus->SetParameter(2, 5.0);   // Sigma iniziale

    // 3. Esegui il Fit con la funzione predefinita
    // Nota: Passiamo l'oggetto TF1 al Fit, non solo il nome "gaus"
    htime30_12->Fit(f_gaus, "R"); // "R" usa il range definito in TF1

    // 4. Configurazione e Disegno (il resto del tuo codice)
    htime30_12->SetLineColor(kAzure+1);
    htime30_12->SetLineWidth(3);
    htime30_12->SetFillColorAlpha(kAzure-4, 0.35);
    htime30_12->SetTitle("time plot 20%");
    htime30_12->Draw("HIST");
    ctime12->Update();

    // 5. Estrai i risultati dal TF1 dopo il fit (usa l'oggetto f_gaus)
    double amp_12   = f_gaus->GetParameter(0);
    double mean_12  = f_gaus->GetParameter(1); // Questo sarà il risultato del fit
    double sigma_12 = f_gaus->GetParameter(2);


    double e_amp_12   = f_gaus->GetParError(0);
    double e_mean_12  = f_gaus->GetParError(1);
    double e_sigma_12 = f_gaus->GetParError(2);


    // ... estrazione degli errori
     TCanvas *ctime23 = new TCanvas("ctime23", "htriple1", 900, 700);
    htime30_23->SetLineColor(kAzure+1);
    htime30_23->SetLineWidth(3);
    htime30_23->SetFillColorAlpha(kAzure-4, 0.35);
    htime30_23->SetTitle("Amplitude");
    htime30_23->Draw("HIST");
    htime30_23->Fit("gaus"); // Fit su htime20_23
    ctime23->Update();

    TF1 *f2 = htime30_23->GetFunction("gaus");
    double amp_23   = f2->GetParameter(0);  ///taking gaus parameters
    double mean_23  = f2->GetParameter(1);
    double sigma_23 = f2->GetParameter(2);

    double e_amp_23   = f2->GetParError(0);
    double e_mean_23  = f2->GetParError(1);
    double e_sigma_23 = f2->GetParError(2);

    TCanvas *ctime13 = new TCanvas("ctime13", "htriple1", 900, 700);
    htime30_13->SetLineColor(kAzure+1);
    htime30_13->SetLineWidth(3);
    htime30_13->SetFillColorAlpha(kAzure-4, 0.35);
    htime30_13->SetTitle("Amplitude");
    htime30_13->Draw("HIST");
    htime30_13->Fit("gaus"); // Fit su htime20_13
    ctime13->Update();

    // Take gaussian from the correct histogram htime20_13
    TF1 *f3 = htime30_13->GetFunction("gaus");
    double amp_13   = f3->GetParameter(0);  // taking gaus parameters
    double mean_13  = f3->GetParameter(1);
    double sigma_13 = f3->GetParameter(2);

    double e_amp_13   = f3->GetParError(0);
    double e_mean_13  = f3->GetParError(1);
    double e_sigma_13 = f3->GetParError(2);

    // --------------- Time resolution (with error propagation) ---------------

    // Use: sigma_12^2 = sigma1^2 + sigma2^2, etc.
    // Inversion:
    // sigma1^2 = 1/2 * ( sigma_12^2 + sigma_13^2 - sigma_23^2 )
    // sigma2^2 = 1/2 * ( sigma_12^2 + sigma_23^2 - sigma_13^2 )
    // sigma3^2 = 1/2 * ( sigma_13^2 + sigma_23^2 - sigma_12^2 )

    double s12 = sigma_12;
    double s13 = sigma_13;
    double s23 = sigma_23;

    double es12 = e_sigma_12;
    double es13 = e_sigma_13;
    double es23 = e_sigma_23;

    double s1sq = 0.5 * ( s12*s12 + s13*s13 - s23*s23 );
    double s2sq = 0.5 * ( s12*s12 + s23*s23 - s13*s13 );
    double s3sq = 0.5 * ( s13*s13 + s23*s23 - s12*s12 );

    if (s1sq < 0) {
        cout << "WARNING: sigma1^2 < 0 (" << s1sq << "). Setting to 0 to avoid NaN." << endl;
        s1sq = 0;
    }
    if (s2sq < 0) {
        cout << "WARNING: sigma2^2 < 0 (" << s2sq << "). Setting to 0 to avoid NaN." << endl;
        s2sq = 0;
    }
    if (s3sq < 0) {
        cout << "WARNING: sigma3^2 < 0 (" << s3sq << "). Setting to 0 to avoid NaN." << endl;
        s3sq = 0;
    }

    double sigma_1 = TMath::Sqrt(s1sq);
    double sigma_2 = TMath::Sqrt(s2sq);
    double sigma_3 = TMath::Sqrt(s3sq);

    // Error propagation:
    // S1 = sigma1^2 = 1/2*(a + b - c), with a = s12^2, b = s13^2, c = s23^2
    // var(a) = (2*s12*es12)^2, etc.
    // var(S1) = (1/2)^2 * ( var(a) + var(b) + var(c) )
    // => var(S1) = s12^2*es12^2 + s13^2*es13^2 + s23^2*es23^2
    // then error on sigma1: e_sigma1 = sqrt(var(S1)) / (2*sigma1)

    double varS1 = s12*s12*es12*es12 + s13*s13*es13*es13 + s23*s23*es23*es23;
    double varS2 = s12*s12*es12*es12 + s23*s23*es23*es23 + s13*s13*es13*es13; // same form, ordering irrelevant
    double varS3 = s13*s13*es13*es13 + s23*s23*es23*es23 + s12*s12*es12*es12;

    double e_sigma_1 = 0.0;
    double e_sigma_2 = 0.0;
    double e_sigma_3 = 0.0;

    if (sigma_1 > 0)
        e_sigma_1 = TMath::Sqrt(varS1) / (2.0 * sigma_1);
    else
        e_sigma_1 = 0.0;

    if (sigma_2 > 0)
        e_sigma_2 = TMath::Sqrt(varS2) / (2.0 * sigma_2);
    else
        e_sigma_2 = 0.0;

    if (sigma_3 > 0)
        e_sigma_3 = TMath::Sqrt(varS3) / (2.0 * sigma_3);
    else
        e_sigma_3 = 0.0;

    cout.setf(std::ios::fixed); cout.precision(4);
    cout << "Results of time resolution with CFD = 20% (obtained from gaussian fits):\n"
         << "sigma12 (fit) = " << s12 << " +/- " << es12 << " ps\n"
         << "sigma23 (fit) = " << s23 << " +/- " << es23 << " ps\n"
         << "sigma13 (fit) = " << s13 << " +/- " << es13 << " ps\n\n"
         << "Derived per-detector resolutions (sigma):\n"
         << "sigma1 = " << sigma_1 << " +/- " << e_sigma_1 << " ps\n"
         << "sigma2 = " << sigma_2 << " +/- " << e_sigma_2 << " ps\n"
         << "sigma3 = " << sigma_3 << " +/- " << e_sigma_3 << " ps\n";

    // -----------------------
    //       Plot results
    // -----------------------
    TCanvas *c = new TCanvas("c", "Amplitude Comparison", 800, 600);
    hAmpAll->SetLineColor(kBlue);
    hAmp1Hit->SetLineColor(kRed);
    hAmpAll->Draw("HIST");
    hAmp1Hit->Draw("HIST SAME");

    auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->AddEntry(hAmpAll, "All events", "l");
    legend->AddEntry(hAmp1Hit, "Single-hit events", "l");
    legend->Draw();

    // -----------------------
    //     Save everything
    // -----------------------
    TFile *fout = new TFile("muon_run_only2_TProfile_run19_trple_channel12_channel5.root", "RECREATE");
    hAmpAll->Write();
    hAmp1Hit->Write();
    hQ->Write();
    hTime20->Write();
    hTime30->Write();
    hTime50->Write();
    hTime60->Write();
    hTime60_2->Write();
    //c->Write();
    c1->Write();
    htriple1->Write();
    htriple2->Write();
    htriple3->Write();
    htime_12->Write();
    //htime_23->Write();
    //htime_13->Write();
    htime20_12->Write();
    //htime20_23->Write();
    //htime20_13->Write();
    htime30_12->Write();
    htime30_23->Write();
    htime30_13->Write();
    hintegral1->Write();
    hintegral2->Write();
    hintegral3->Write();
    hcharge1->Write();
    hcharge2->Write();
    hrisetime1->Write();
    hrisetime2->Write();
    hpeaktime1->Write();
    hpeaktime2->Write();
    timevsamplitude->Write();
    timevsamplitude2->Write();
    timevsamplitude3->Write();
    hbaseline1->Write();
    hbaseline2->Write();
    timevsintegral->Write();
    prof->Write();
    prof2->Write();
    prof3->Write();
    cmediabin->Write();
    cfdvschannel1->Write();
    cintegral->Write();
    cintegral2->Write();
    tout->Write();
   // copy->Write();
    conteggixcanale1->Write();
    conteggixcanale2->Write();
    conteggixcanale3->Write();
    cAmpVsTime->Write();


    fout->Close(); // Close output file
    file->Close(); // Close input file
    cout << "Analysis complete. Results saved to muon_run_only2_TProfile_analysis.root" << endl;

    return 0;
}
