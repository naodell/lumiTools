#define vdmAnalyzer_cxx

#include "vdmAnalyzer.h"

int const NBUCKETS = 3564;
int const NORBITS  = 4096;
int const NMAXCHANNELS = 4;

unsigned const  time_i = 0.; //6.e6;
unsigned const  time_f = 1.5e7;

unsigned const  bunches[] = {1, 41, 81, 121, 161, 201, 241, 
    721, 761, 801, 841, 881, 921, 961,
    //1441, 1481, 1521, 1561, 1601, 1641, 1681, // April VdM
    1581, 1621, 1661, 1701, 1741, 1781, 1821, // July/Nov VdM
    2161, 2201, 2241, 2281, 2321, 2361, 2401, 
    2881 // Nov VdM
    //2881, 2921, 2961, 3001, 3041, 3081, 3121 // April/July VdM
};
unsigned const offsets[] = {2, 1, -1, 2};

void vdmAnalyzer::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    outFile = new TFile("vdmHistograms.root", "RECREATE");
    
    outFile->mkdir("Total");
    for (unsigned i = 0; i < sizeof(bunches)/sizeof(UInt_t); ++i) outFile->mkdir(Form("BCID_%i", bunches[i]));

    outFile->cd("Total");
    h1_EntryDuration    = new TH1F("h1_EntryDuration", "Time duration per entry", 10, 0., 1000.);
    h1_Total            = new TH1F("h1_Total", "Total hits from PLT FED", 10, 0., 1000.);

    h1_Rates[0]         = new TH1F("h1_Rates_OR", "BCM1F OR;BCID;N_{hits}", NBUCKETS, -0.5, float(NBUCKETS)-0.5);
    h1_Rates[1]         = new TH1F("h1_Rates_AND", "BCM1F AND ;BCID;N_{hits}", NBUCKETS, -0.5, float(NBUCKETS)-0.5);
    h1_Rates[2]         = new TH1F("h1_Rates_XOR1", "BCM1F XOR1;BCID;N_{hits}", NBUCKETS, -0.5, float(NBUCKETS)-0.5);
    h1_Rates[3]         = new TH1F("h1_Rates_XOR2", "BCM1F XOR2 ;BCID;N_{hits}", NBUCKETS, -0.5, float(NBUCKETS)-0.5);

    p1_RateVsTimeSum[0] = new TProfile("p1_RateVsTime_OR", "BCM1F OR rates vs time (all BCs);time (ms);rate (Hz)", (time_f - time_i)/2000, time_i, time_f);
    p1_RateVsTimeSum[1] = new TProfile("p1_RateVsTime_AND", "BCM1F AND rates vs time (all BCs);time (ms);rate (Hz)", (time_f - time_i)/2000, time_i, time_f);
    p1_RateVsTimeSum[2] = new TProfile("p1_RateVsTime_XOR1", "BCM1F XOR1 rates vs time (all BCs);time (ms);rate (Hz)", (time_f - time_i)/2000, time_i, time_f);
    p1_RateVsTimeSum[3] = new TProfile("p1_RateVsTime_XOR2", "BCM1F XOR2 rates vs time (all BCs);time (ms);rate (Hz)", (time_f - time_i)/2000, time_i, time_f);

    for (unsigned i = 0; i < sizeof(bunches)/sizeof(UInt_t); ++i) {
        outFile->cd(Form("BCID_%i", bunches[i]));
        p1_RateVsTime[i][0] = new TProfile("p1_RateVsTime_OR", Form("BCM1F OR rates vs time (BCID %i)", bunches[i]), (time_f - time_i)/2000, time_i, time_f);
        p1_RateVsTime[i][1] = new TProfile("p1_RateVsTime_AND", Form("BCM1F AND rates vs time (BCID %i)", bunches[i]), (time_f - time_i)/2000, time_i, time_f);
        p1_RateVsTime[i][2] = new TProfile("p1_RateVsTime_XOR1", Form("BCM1F XOR1 rates vs time (BCID %i)", bunches[i]), (time_f - time_i)/2000, time_i, time_f);
        p1_RateVsTime[i][3] = new TProfile("p1_RateVsTime_XOR2", Form("BCM1F XOR2 rates vs time (BCID %i)", bunches[i]), (time_f - time_i)/2000, time_i, time_f);
    }

    prevTime    = 0;    
    eventCount  = 0;
}

Bool_t vdmAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry);

    Time = Time + offset;
    prevTime = Time;

    if (unsigned(Time) < time_i || unsigned(Time) > time_f) return kTRUE;

    if (eventCount%1000 == 0 and eventCount != 0) {
        cout << eventCount << " entries processed..." << endl;
    }

    UInt_t timeDuration = 333; //(Time-prevTime);
    h1_EntryDuration->Fill(timeDuration);

    for (unsigned i = 0; i < unsigned(NMAXCHANNELS); ++i) {
        float ratePerOrbit = 0;
        UInt_t bunchCount = 0;

        for (unsigned j = 0; j < unsigned(NBUCKETS); ++j) {
            unsigned rate = (1000./timeDuration)*(Hist[i][j]);
            h1_Rates[i]->Fill(j, Hist[i][j]); 

            if (j == bunches[bunchCount] + offsets[i]) {
                ratePerOrbit += rate;
                p1_RateVsTime[bunchCount][i]->Fill(Time, rate);
                //cout << rate << ", " << Time << endl;
                ++bunchCount;
            }
        }
        
        p1_RateVsTimeSum[i]->Fill(Time, ratePerOrbit/bunchCount);
    }
    ++eventCount;

    return kTRUE;
}

void vdmAnalyzer::Terminate()
{
    //for (int i = 0; i < NMAXCHANNELS; ++i) h1_Rates[i]->Scale(1000./float(time_f - time_i));

    outFile->Write();
    outFile->Close();
}
