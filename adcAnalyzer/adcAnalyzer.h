//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 26 07:52:46 2012 by ROOT version 5.27/06b
// from TTree /
// found on file: Memory Directory
//////////////////////////////////////////////////////////

#ifndef adcAnalyzer_h
#define adcAnalyzer_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include "TDatime.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
//#include "TSpectrum.h"
#include "TVirtualFFT.h"

class adcAnalyzer : public TSelector {
    private: 
		 TFile* 		outFile;

         // Event
		 TH1D*		h1_ADCSpectrumCombined;
		 TH1D*		h1_TriggerTime;

		 // Channel-by-Channel
		 TH1D*		h1_ADCSpectrum[8];
		 TH1D*		h1_Baseline[8];
		 TH1D*		h1_PeakValue[8];

		 TH1D*		h1_TimeOverShoot[8];
		 TH1D*		h1_TimeOverShootPO[8];
		 TH1D*		h1_TimeOverThresh[8];
		 TH1D*		h1_TimeOverThreshPO[8];
		 TH1D*		h1_EffectiveDeadtime[8];

         TH1D*      h1_HitsPerOrbit[8];
		 TH1D*		h1_HitsVsTime[8];
         TProfile*  p1_RateVsTime[8];
		 TProfile*  p1_ADCVsTime[8][20];

         TH2D*      h2_PulseHeightVsTOT[8];
         TH2D*      h2_OSHeightVsTOT[8];
         TH2D*      h2_TOSVsTOT[8];

		 TGraph*    g1_dVdtVsTime[8];
		 TGraph*    g1_ADCVsTime[8];
		 TProfile*	p1_BaselineVsTime[8];

         // For spectrum analysis
         //TSpectrum* spectralAnalyzer;

         // For counting events
         UInt_t     nEvent;
         UInt_t     nEntries;

         // Trigger
         Float_t    triggerTime;

    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types
        UInt_t time_sec;
        UInt_t time_usec;
        UInt_t acqStatus;
        UInt_t vmeStatus;
        UInt_t channel;
        UInt_t event;
        UInt_t eventHeader[4];
        Double_t baseMean;
        Double_t baseSigma;
        UInt_t eventCount;
        UChar_t data[0x200000];

        // List of branches
        TBranch* b_time_sec;
        TBranch* b_time_usec;
        TBranch* b_acqStatus;
        TBranch* b_vmeStatus;
        TBranch* b_channel;
        TBranch* b_event;
        TBranch* b_eventHeader;
        TBranch* b_baseMean;
        TBranch* b_baseSigma;
        TBranch* b_eventCount;
        TBranch* b_data;

        adcAnalyzer(TTree * /*tree*/ =0) { }

        virtual ~adcAnalyzer() { }
        virtual Int_t   Version() const {return 2;}
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) {return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;}
        virtual void    SetOption(const char *option) {fOption = option;}
        virtual void    SetObject(TObject *obj) {fObject = obj;}
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate();
        virtual void    Terminate();

        ClassDef(adcAnalyzer,0);
};

#endif

#ifdef adcAnalyzer_cxx
void adcAnalyzer::Init(TTree *tree)
{

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("time_sec",&time_sec, &b_time_sec);
    fChain->SetBranchAddress("time_usec",&time_usec, &b_time_usec);
    fChain->SetBranchAddress("ACQ_Status",&acqStatus, &b_acqStatus);
    fChain->SetBranchAddress("VME_Status",&vmeStatus, &b_vmeStatus);
    fChain->SetBranchAddress("channel",&channel, &b_channel);
    fChain->SetBranchAddress("Event",&event, &b_event);
    fChain->SetBranchAddress("EventHeader1",&eventHeader[0], &b_eventHeader);
    fChain->SetBranchAddress("EventHeader2",&eventHeader[1], &b_eventHeader);
    fChain->SetBranchAddress("EventHeader3",&eventHeader[2], &b_eventHeader);
    fChain->SetBranchAddress("EventHeader4",&eventHeader[3], &b_eventHeader);
    fChain->SetBranchAddress("BaseMean",&baseMean, &b_baseMean);
    fChain->SetBranchAddress("BaseSigma",&baseSigma, &b_baseSigma);
    fChain->SetBranchAddress("N",&eventCount, &b_eventCount);
    fChain->SetBranchAddress("data",data, &b_data);

}

Bool_t adcAnalyzer::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

#endif // #ifdef adcAnalyzer_cxx
