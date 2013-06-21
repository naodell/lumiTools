#define adcAnalyzer_cxx

#include "adcAnalyzer.h"
#include <TStyle.h>

// Declare global variables here //
const double adcBins[]      = {2., 4., 8., 12., 16., 24., 32., 40., 50., 64., 100.};
//const double threshold[]    = {2., 8/3.9, 8/3.9, 6/3.9, 6/3.9, 6/3.9, 8/3.9, 6/3.9};
const double threshold[]    = {3.7, 3.26, 3.26, 3.26, 3.04, 3.04, 4.56, 3.04};
const unsigned offset       = 3225;

// VDM time
//const TDatime time_i = "2012-04-16 18:33:45";
//const TDatime time_f = "2012-04-16 22:40:45";

// Fill 2536
const TDatime time_i = "2012-04-20 15:33:45";
const TDatime time_f = "2012-04-20 22:40:45";

void adcAnalyzer::Begin(TTree * /*tree*/)
{
    TString option = GetOption();
    triggerTime = 0;
    nEvent      = 0;

    // Tool initializations
    //spectralAnalyzer = new TSpectrum();

    // Output file
    outFile = new TFile("adcHistograms.root", "RECREATE");
    for (int i = 0; i < 8; ++i) outFile->mkdir(Form("ADC_CH%i", i), Form("ADC_CH%i", i)); 

    // Histogram initializations
    TH1::SetDefaultSumw2(kTRUE);
    h1_ADCSpectrumCombined	= new TH1D("h1_ADCSpectrumCombined", "ADC spectrum (all channels);ADC counts;N_{events}", 255, 0., 255.);
    h1_TriggerTime          = new TH1D("h1_TriggerTime", "Turn clock time", 100, 0., 100.);

    for (unsigned i = 0; i < 8; ++i) {
        outFile->cd(Form("ADC_CH%i", i));
        p1_RateVsTime[i]            = new TProfile(Form("p1_RateVsTime_%i", i), Form("Rate vs time (channel %i);time;rate", i), (time_f.Convert() - time_i.Convert())/2, time_i.Convert(), time_f.Convert()); 

        h1_ADCSpectrum[i]           = new TH1D(Form("h1_ADCSpectrum_%i", i), Form("ADC spectrum (channel %i);ADC counts;N_{events}", i), 130, -40.5, 90.5);
        h1_Baseline[i]              = new TH1D(Form("h1_Baseline_%i", i), Form("ADC baseline (channel %i);ADC counts;N_{events}", i), 255, 0.5, 255.5);
        h1_PeakValue[i]             = new TH1D(Form("h1_PeakValue_%i", i), Form("Signal peak (channel %i);ADC counts;N_{events}", i), 100, 0.5, 100.5);

        h1_TimeOverShoot[i]         = new TH1D(Form("h1_TimeOverShoot_%i", i), Form("Time in overshoot (channel %i);time (ns);N_{events}", i), 500, 0.5, 4000.5);
        h1_TimeOverShootPO[i]       = new TH1D(Form("h1_TimeOverShootPO_%i", i), Form("Time in overshoot per orbit (channel %i);time (ns);N_{events}", i), 750, 0.5, 12000.5);

        h1_TimeOverThresh[i]        = new TH1D(Form("h1_TimeOverThresh_%i", i), Form("Time over threshold (channel %i);time (ns); nEvents", i), 250, 0.5, 2000.5);
        h1_TimeOverThreshPO[i]      = new TH1D(Form("h1_TimeOverThreshPO_%i", i), Form("Time over threshold per orbit (channel %i);time (ns); nEvents", i), 625, 0.5, 10000.5);

        h1_HitsVsTime[i]            = new TH1D(Form("h1_HitsVsTime_%i", i), Form("ADC Channel %i;time (ns);N_{hits}", i), 11250, 0., 89000);
        h1_HitsPerOrbit[i]          = new TH1D(Form("h1_HitsPerOrbit_%i", i), Form("ADC Channel %i;N_{hits};N_{orbits}", i), 800, -0.5, 799.5);
        h1_EffectiveDeadtime[i]     = new TH1D(Form("h1_EffectiveDeadtime_%i", i), Form("ADC Channel %i;t_{OT+OS};N_{hits}", i), 600, 0.5, 6000.5);

        h2_PulseHeightVsTOT[i]      = new TH2D(Form("h2_PulseHeightVsTOT_%i", i), Form("pulse height vs. TOT (channel %i);time (ns); h_{pulse} (ADC counts)", i), 250, 0.5, 2000.5, 80, 0.5, 80.5);
        h2_OSHeightVsTOT[i]         = new TH2D(Form("h2_OSHeightVsTOT_%i", i), Form("overshoot height vs. TOT (channel %i);time (ns); h_{pulse} (ADC counts)", i), 250, 0.5, 2000.5, 40, 0.5, 40.5);
        h2_TOSVsTOT[i]              = new TH2D(Form("h2_TOSVsTOT_%i", i), Form("TOS vs. TOT (channel %i);TOT (ns); TOS", i), 250, 0.5, 2000.5, 500, 0.5, 4000.5);

        outFile->GetDirectory(Form("ADC_CH%i", i))->mkdir("Orbit scans", "Orbit Scans");

        for (unsigned j = 0; j < 20; ++j) {
            outFile->cd(Form("ADC_CH%i/Orbit scans", i));
            p1_ADCVsTime[i][j]      = new TProfile(Form("p1_ADCVsTime_%i_%i", i, j), Form("ADC Channel %i (acquisition %i);time (ns);ADC counts", i, 50*j), 22500, 0., 90000, 0., 255.);  
        }
    }
}

void adcAnalyzer::SlaveBegin(TTree * /*tree*/)
{
    TString option = GetOption();
}

Bool_t adcAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry);

    //if (time_sec < time_i.Convert() || time_sec > time_f.Convert()) return kTRUE;

    Int_t spectrum[eventCount], time[eventCount];
    Float_t dVdt[eventCount];

    // Counting variables
    Float_t channelBaseline = 0;
    Float_t peakValue       = 0;
    Float_t osValue         = 0;
    Float_t tot             = 0;
    Float_t overTime        = 0;
    Float_t overTimeTotal   = 0;
    Float_t underTime       = 0;
    Float_t underTimeTotal  = 0;
    UInt_t  nHits           = 0;
    UInt_t  crossingTime    = -1;
    bool    overThr         = false;
    bool    underThr        = false;
    bool    isRising        = false;
    bool    isDropping      = false;

    // Calculate baseline
    for (unsigned i = (offset - 500); i < offset; ++i) {
        channelBaseline += data[i];
    }
    channelBaseline /= 500;
    h1_Baseline[channel]->Fill(channelBaseline);
    //channelBaseline = h1_Baseline[channel]->GetMaximum();

    dVdt[0] = 0.;
    for (int i = 0; i < (int)eventCount; ++i) {
        time[i]     = 2*i;
        spectrum[i] = Int_t(data[i]);
        Float_t pulseHeight = channelBaseline - spectrum[i];
        if (i > 0) dVdt[i] = -1*(spectrum[i] - spectrum[i-1]);

        /*
        if (dVdt < 0) isDropping = true;
        else isDropping = false;

        if (dVdt > 0) isRising = true;
        else isRising = false;
        */

        if (pulseHeight > threshold[channel]) {
            if (!overThr) {

                // Record hits for TOA histogram
                if (entry%8 != 0 && i > (int)offset) {
                    //if (channel == 1) cout << time[i] << ", ";
                    crossingTime = time[i] - 2*offset;// - triggerTime);
                } else if (entry%8 == 0 && i < 1e4) {
                    crossingTime = time[i] - 2*offset;

                    triggerTime = 2*i;
                    h1_TriggerTime->Fill(triggerTime);
                }
                ++nHits;
            }

            overThr     = true;
            overTime    += 2;
            if (peakValue < pulseHeight) peakValue = pulseHeight;
        } else {
            if (overTime >= 14.) {
                h1_TimeOverThresh[channel]->Fill(overTime);
                h1_PeakValue[channel]->Fill(peakValue);
                h1_HitsVsTime[channel]->Fill(crossingTime);
                h2_PulseHeightVsTOT[channel]->Fill(overTime, peakValue);

                if (overTime > 200) tot = overTime;

                overTimeTotal += overTime;
            }

            overThr         = false;
            overTime        = 0;
            peakValue       = 0;
            crossingTime    = -1;
        }

        if (pulseHeight <= -threshold[channel]){
            underThr    = true;
            underTime   += 2.;

            if (osValue > pulseHeight) osValue = pulseHeight;
        } else {
            if (underTime >= 14) {
                h1_TimeOverShoot[channel]->Fill(underTime);
                underTimeTotal += underTime;

                if (tot > 0) {
                    h1_EffectiveDeadtime[channel]->Fill(tot+underTime);
                    h2_OSHeightVsTOT[channel]->Fill(tot, fabs(osValue));
                    h2_TOSVsTOT[channel]->Fill(tot, underTime);
                    tot = 0;
                }
            }
            underThr    = false;
            underTime   = 0;
        }

        h1_ADCSpectrum[channel]->Fill(pulseHeight);

        unsigned iScan = unsigned(floor(nEvent/100));
        if (nEvent%100 == 0 && iScan < 20) {
            p1_ADCVsTime[channel][iScan]->Fill(time[i], data[i]);
        }
    }

    //if (nHits > 0 && channel == 1) cout << endl;

    p1_RateVsTime[channel]->Fill(time_sec, nHits/89e-5);

    //if (entry < 8) {
    //    g1_dVdtVsTime[channel]  = new TGraph(eventCount, time, dVdt);
    //    outFile->GetDirectory(Form("ADC_CH%i", channel))->Append(g1_dVdtVsTime[channel]);
    //}
    
    h1_TimeOverThreshPO[channel]->Fill(overTimeTotal);
    h1_TimeOverShootPO[channel]->Fill(underTimeTotal);
    h1_HitsPerOrbit[channel]->Fill(nHits);

    if (entry%8 == 7) {
        ++nEvent;
        if (nEvent % 500 == 0 && nEvent != 0) cout << nEvent << " events processed..."<< endl;
    }

    return kTRUE;
}

void adcAnalyzer::SlaveTerminate()
{

}

void adcAnalyzer::Terminate()
{
    //for (unsigned i = 0; i < 8; ++i){
    //    h1_TimeOverThreshPO[i]->Scale(1/nEvent);
    //    h1_TimeOvershootPO[i]->Scale(1/nEvent);
    //    h1_HitsPerOrbit[i]->Scale(1/nEvent);
    //}

    outFile->Write();
    outFile->Close();
}
