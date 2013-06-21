//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 23 14:47:03 2012 by ROOT version 5.32/01
// from TTree PLTHist/PLTHist
// found on file: data/VDM/Data_Histograms_20120416.161104.root
//////////////////////////////////////////////////////////

#ifndef vdmAnalyzer_h
#define vdmAnalyzer_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include <TROOT.h>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include "TH1F.h"
#include "TH2.h"
#include "TProfile.h"

class vdmAnalyzer : public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types
        Int_t           Time;
        Int_t           Total;
        Int_t           Hist[9][3564];

        // List of branches
        TBranch        *b_Time;   //!
        TBranch        *b_Total;   //!
        TBranch        *b_Hist;   //!

        // Output file with histograms, etc.
        TFile*          outFile;

        TH1F*           h1_EntryDuration;
        TH1F*           h1_Total;
        TH1F*           h1_Rates[6];

        TProfile*       p1_RateVsTimeSum[6];
        TProfile*       p1_RateVsTime[35][6];

        // Event counting, etc.
        UInt_t          eventCount;
        UInt_t          prevTime;
        UInt_t          offset;

        // Function declarations
        vdmAnalyzer(TTree * /*tree*/ =0) : fChain(0) { }
        virtual ~vdmAnalyzer() { }
        virtual Int_t   Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree) {TString option = GetOption();}
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate() {};
        virtual void    Terminate();

        ClassDef(vdmAnalyzer,0);
};

#endif

#ifdef vdmAnalyzer_cxx
void vdmAnalyzer::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("Time", &Time, &b_Time);
    fChain->SetBranchAddress("Total", &Total, &b_Total);
    fChain->SetBranchAddress("Hist", Hist, &b_Hist);

}

Bool_t vdmAnalyzer::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    cout << "\nOpening next file in chain...\n" << endl;
    offset = prevTime;

    return kTRUE;
}

#endif // #ifdef vdmAnalyzer_cxx
