#! /usr/bin/env python
import sys, math, csv, json, pprint, cProfile, os
from array import array
import ROOT as r
import CountingExperiment as c
import lumiTools as l

r.gStyle.SetOptStat(0)
algos = ['OR', 'XOR+', 'XOR-', 'AND']

### Set MC parameters here ###
doEff       = True
nOrbits     = 1000
timeStep    = 25

#probs       = [0.06, 0.06, 0.085, 0.083, 0.115, 0.095, 0.075, 0.095]
probs       = [0.06, 0.06, 0.085, 0.083, 0.115, 0.095, 0.075, 0.0]
probScale   = 1.

#fillPattern = l.get_fill_pattern('1000ns_50b_35_14_35')
#fillPattern = l.get_fill_pattern('50ns_840b_807_0_816_108bpi12inj')
fillPattern = l.get_fill_pattern('50ns_1380b_1331_0_1320_144bpi12inj')

bpo         = len(fillPattern)

### Initialize histograms, etc. ###
#outFile = r.TFile('histograms/MC/mc_test_1.root', 'RECREATE')
outFile = r.TFile('histograms/MC/test.root', 'RECREATE')
r.TH1.SetDefaultSumw2(r.kTRUE)

outFile.mkdir('single_channel')
outFile.mkdir('coincidences')

h1_hits = []
h1_probabilities = []
outFile.cd('single_channel')
for i in range(8):
    h1_hits.append(r.TH1D('h1_hits_{}'.format(i), 'BCM1F MC (ch {});BCID;hits'.format(i), bpo, 0.5, float(bpo+0.5)))
    h1_probabilities.append(r.TH1D('h1_probabilities_{}'.format(i), 'BCM1F MC (ch {});BCID;P(hit)'.format(i), bpo, 0.5, float(bpo+0.5)))

h1_coins = {}
outFile.cd('coincidences')
for algo in algos:
    h1_coins[algo] = r.TH1D('h1_coins_{}'.format(algo), 'BCM1F MC (algo {});BCID;P(hit)'.format(algo), bpo, 0.5, float(bpo+0.5))

### Run the experiment ###

print '\nStarting pseudoexperiments...\n'

lumiMC = c.CountingExperiment(timeStep = 25, nIter = bpo, pList = probs, sourceFile = 'histograms/adcHistograms.root', doEff = doEff)
lumiMC.set_fill_pattern(fillPattern)
lumiMC.initialize_efficiency_histograms()

for i in range(nOrbits):
    orbit = lumiMC.run_bunch_by_bunch()

    for j in range(bpo):
        for channel in range(8):
            if orbit[0][str(channel)][j]:
                h1_probabilities[channel].Fill(j+1, 1/float(nOrbits))
                h1_hits[channel].Fill(j+1)

        for algo in algos:
            if orbit[1][algo][j]: h1_coins[algo].Fill(j+1, 1/float(nOrbits))

    if (i+1)%500 is 0 and i is not 0: print '{} orbits generated...'.format(i+1)

#h1_coins['OR'].Scale(0.545/0.23)

outFile.Write()
outFile.Close()
