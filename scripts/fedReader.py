#! /usr/bin/env python
import sys, math, csv, json, pprint
from array import array
import ROOT as r
import lumiTools as l

### Get Data/MC files ###
bcmData = []
bcmData.extend(l.read_csv_file('data_db/BunchByBunch/Fill2663_NoHeader1.csv'))

### Set time frame ###
timeStart = 0
timeStop  = int(len(bcmData)*1.1)
duration  = int(timeStop - timeStart)

### Declare histograms ###
outFile = r.TFile('histograms/fedHistograms.root', 'RECREATE')
outFile.mkdir('OR', 'OR')
#outFile.mkdir('XOR-', 'XOR-')
#outFile.mkdir('XOR+', 'XOR+')
#outFile.mkdir('AND', 'AND')

r.TH1.SetDefaultSumw2(r.kTRUE)
r.TProfile.SetDefaultSumw2(r.kTRUE)

h1_Intensities          = r.TH1D('h1_Intensities', 'I_{1}#times I_{2};I_{1}#times I_{2};N_{bx}', 50, 0.8, 1.8) 
h1_IntensityCorr        = r.TH1D('h1_IntensityCorr', 'I_{1}#times I_{2}/<I_{1}#times I_{2}>;BCID;I_{1}#times I_{2}/<I_{1}#times I_{2}>', 3564, 0.5, 3564.5) 
p1_IntensVsTime         = r.TProfile('p1_IntensVsTime', 'I_{1}#times I_{2};BCID;I_{1}#times I_{2}', 3564, 0.5, 3564) 

outFile.cd('OR')
h1_ProbabilityOr        = r.TH1D('h1_ProbabilityOr', 'Probability of BCM1F OR;BCID;p(hit)', 3564, 0.5, 3564.5) 
h1_ProbabilityOrSpec    = r.TH1D('h1_ProbabilityOrSpecific', 'Probability of BCM1F OR;BCID;p(hit)/I1*I2)', 3564, 0.5, 3564.5) 
h1_LumiOr               = r.TH1D('h1_LumiOr', 'BCM1F OR;BCID;Lumi', 3564, 0.5, 3564.5) 
h1_LumiOrSpec           = r.TH1D('h1_LumiOrSpecific', 'BCM1F OR;BCID;Lumi/I1*I2)', 3564, 0.5, 3564.5) 

h1_CorrectionsOr        = r.TH1D('h1_Corrections', 'Corrections BCM1F OR;BCID;correction', 3564, 0.5, 3564.5) 

bunches = 1380
time    = 607
orbit   = 11246
count   = 1220681 #time*orbit
avgProb = 0.3788

### Get average intensity ###
avgIntensity = 0.
for row in bcmData:
    bunch       = int(row['Bunch']) + 1
    intensity   = float(row['I1I2'])*1e-22

    if intensity > 0:
        avgIntensity += intensity
        p1_IntensVsTime.Fill(bunch, intensity)

avgIntensity /= bunches

### Fill histograms ###
for i, row in enumerate(bcmData):
    if i > 3560: break
    bunch       = int(row['Bunch']) + 1
    intensity   = float(row['I1I2'])*1e-22
    orProbLut   = float(bcmData[i+2]['LUTORCounts'])/float(count)
    hfLumi      = float(row['HFLumi'])
    intensCorr  = intensity/avgIntensity

    
    h1_ProbabilityOr.Fill(bunch, orProbLut)
    h1_ProbabilityOr.SetBinError(bunch, math.sqrt(orProbLut/count))
    h1_LumiOr.Fill(bunch, hfLumi)

    if intensity > 0:
        #print orProbLut, intensity
        h1_ProbabilityOrSpec.Fill(bunch, orProbLut/intensCorr)
        h1_ProbabilityOrSpec.SetBinError(bunch, math.sqrt(orProbLut/count)/intensCorr)
        h1_LumiOrSpec.Fill(bunch, hfLumi/intensCorr)

        h1_CorrectionsOr.Fill(bunch, avgProb/(orProbLut/intensCorr))
        h1_CorrectionsOr.SetBinError(bunch, math.sqrt(orProbLut/count))

        h1_IntensityCorr.Fill(bunch, intensCorr)
        h1_IntensityCorr.SetBinError(bunch, 0.01)


outFile.Write()
outFile.Close()
