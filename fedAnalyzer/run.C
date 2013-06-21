#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void run(int nEvents = 1e9, int entry0 = 0) {
    TChain *pltChain = new TChain("PLTHist");

    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.074454.root");
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.092422.root");
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.110349.root");
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.124317.root");
    
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.192115.root");
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.210042.root");
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120719.224010.root");

    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120416.181546.root");
    //pltChain->Add("../data/VDM/BCM1F/Data_Histogram_20120416.195513.root");

    pltChain->Add("../data/Data_Histogram_20121123.194021.root");
    pltChain->Add("../data/Data_Histogram_20121123.211948.root");
    pltChain->Add("../data/Data_Histogram_20121123.225916.root");

    TStopwatch timer;
    timer.Start();

    pltChain->Process("vdmAnalyzer.C+", "", nEvents, entry0);

    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() << endl;
    cout << "RealTime : " << timer.RealTime() << endl;
    cout << "\n";

}
