
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void run(int nEvents = 100000) {

    TChain* fChain = new TChain("T");

    ifstream sourceFiles("input.txt");
    char line[128];
    int  count = 0;

    while (sourceFiles >> line) {
        fChain->Add(line);      
        ++count;
    }

    cout << count << " files added!" << endl;
    sourceFiles.close();

    TStopwatch timer;
    timer.Start();

    fChain->Process("adcAnalyzer.C+", "", 8*nEvents, 0);

    timer.Stop();

    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() << " s (total), " << timer.CpuTime()/nEvents << " s (per event)" << endl;
    cout << "RealTime : " << timer.RealTime() << " s (total), " << timer.RealTime()/nEvents << " s (per event)" << endl;
    cout << "\n";
}

