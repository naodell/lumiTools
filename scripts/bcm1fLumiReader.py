#! /usr/bin/env python
import sys, math, csv, json, pprint
from array import array
import ROOT as r
import lumiTools as l


#fillList = [2629, 2630, 2634, 2535, 2644, 2645, 2646, 2648, 2649, 2651, 2653, 2657, 2658]
fillList = [2520]

bcmData = []
#bcmData.extend(l.read_csv_file('data_db/Current_bcm.csv'))
bcmData.extend(l.read_csv_file('data_db/fill2520bcm.csv'))

#pixelData = []
#pixelData.extend(l.read_json_file('data/pixel_lumi_raw_2012.json'))

### Declare histograms ###

timeStart = 0
timeStop  = int(len(bcmData)*1.1)
duration  = int(timeStop - timeStart)

start     = r.TDatime('2012-05-15 12:28:44').Convert()
stop      = r.TDatime(bcmData[int(len(bcmData))-1]['DT']).Convert()
#stop      = r.TDatime(bcmData[int(len(bcmData))-1]['DT_UNTIL']).Convert()

outFile     = r.TFile('histograms/lumiHistograms.root', 'RECREATE')
outFile.mkdir('OR', 'OR')
outFile.mkdir('XOR-', 'XOR-')
outFile.mkdir('XOR+', 'XOR+')
outFile.mkdir('AND', 'AND')
outFile.mkdir('HFLumi', 'HFLumi')
outFile.mkdir('FBCT', 'FBCT')
outFile.mkdir('PerFill', 'PerFill')

r.TH1.SetDefaultSumw2(r.kTRUE)
r.TProfile.SetDefaultSumw2(r.kTRUE)

outFile.cd('OR')
p1_RateOr           = r.TProfile('p1_RateOr', 'Rate of BCM1F OR;time (GMT);rate (Hz)', duration, timeStart, timeStop) 
p1_MuOr             = r.TProfile('p1_MuOr', '#mu_{OR};time (GMT);#mu_{OR}', duration, timeStart, timeStop) 
p1_LumiOr           = r.TProfile('p1_LumiOr', 'Lumi (OR);time (GMT);Lumi (10^{30} cm^{-2}s^{-1})', duration, start, stop) 
p1_LumiOrByHF       = r.TProfile('p1_LumiOrByHF', 'Lumi (OR);time (GMT);Lumi_{OR}/Lumi_{HF}', duration, start, stop) 
p1_RateOrRatio      = r.TProfile('p1_RateOrRatio', 'Ratio of NIM/LUT rates;time (GMT);R_{NIM}/R_{LUT}', duration, start, stop) 

outFile.cd('AND')
p1_RateAnd          = r.TProfile('p1_RateAnd', 'Rate of BCM1F AND;time (GMT);rate (Hz)', duration, timeStart, timeStop) 
p1_MuAnd            = r.TProfile('p1_MuAnd', '#mu_{AND};time (GMT);#mu_{AND}', duration, timeStart, timeStop) 
p1_MuAnd_Alt        = r.TProfile('p1_MuAnd_Alt', '#mu_{AND} (Alternative);time (GMT);#mu_{AND}', duration, timeStart, timeStop) 

outFile.cd('XOR-')
p1_RateXOr1         = r.TProfile('p1_RateXOr1', 'Rate of BCM1F XOR-;time (GMT);rate (Hz)', duration, timeStart, timeStop) 
p1_MuXOr1           = r.TProfile('p1_MuXOr1', '#mu_{XOR1};time (GMT);#mu_{XOR1}', duration, timeStart, timeStop) 
p1_RateXOr1ByHF     = r.TProfile('p1_RateXOr1ByHF', 'Rate_{XOR-}/L_{HF};time (GMT);Rate_{XOR-}/L_{HF}', duration, timeStart, timeStop) 

outFile.cd('XOR+')
p1_RateXOr2         = r.TProfile('p1_RateXOr2', 'Rate of BCM1F XOR+;time (GMT);rate (Hz)', duration, timeStart, timeStop) 
p1_MuXOr2           = r.TProfile('p1_MuXOr2', '#mu_{XOR2};time (GMT);#mu_{XOR2}', duration, timeStart, timeStop) 
p1_RateXOr2ByHF     = r.TProfile('p1_RateXOr2ByHF', 'Rate_{XOR+}/L_{HF};time (GMT);Rate_{XOR+}/L_{HF}', duration, timeStart, timeStop) 


outFile.cd('HFLumi')
h2_HFLumiVsMuOr         = r.TH2D('h2_HFLumiVsMuOr', 'HF Lumi Vs #mu_{OR};HF Lumi (10^{30} cm^{-2}s^{-1}); #mu_{OR}', 100, 0., 10., 100, 0.1, 1.)  
h2_HFLumiVsMuAnd        = r.TH2D('h2_HFLumiVsMuAnd', 'HF Lumi Vs #mu_{AND};HF Lumi (10^{30} cm^{-2}s^{-1}); #mu_{AND}', 100, 0., 10., 100, 0., 0.1) 

p1_HFLumi               = r.TProfile('p1_HFLumi', 'HF Lumi;time (GMT); 10^{30} cm^{-2}s^{-1}', duration, start, stop) 
p1_HFLumiVsMuAnd        = r.TProfile('p1_HFLumiVsMuAnd', 'HF Lumi Vs #mu_{AND};HF Lumi (10^{30} cm^{-2}s^{-1}); #mu_{AND}', 100, 0., 10.) 
p1_HFLumiVsMuAnd_Alt    = r.TProfile('p1_HFLumiVsMuAnd_Alt', 'HF Lumi Vs #mu_{AND};HF Lumi (10^{30} cm^{-2}s^{-1}); #mu_{AND}', 100, 0., 10.) 

outFile.cd('FBCT')
p1_BeamIntensity            = r.TProfile('p1_BeamIntensity', 'I_{b1} x I_{b2};time (GMT); a.u.', duration, timeStart, timeStop) 
p1_BeamIntensityVsMuOr      = r.TProfile('p1_BeamIntensityVsMuOr', 'I_{b1} x I_{b2} Vs #mu_{OR};I_{b1} x I_{b2}; #mu_{OR}', 100, 0.4, 0.8)  
p1_BeamIntensityVsMuAnd     = r.TProfile('p1_BeamIntensityVsMuAnd', 'I_{b1} x I_{b2} Vs #mu_{AND};I_{b1} x I_{b2}; #mu_{AND}', 100, 0.4, 0.8) 

outFile.cd('PerFill')
p1_OrCorrections    = r.TProfile('p1_OrCorrections', 'L_{HF}/L_{#mu_{OR}} (1380);#mu_{OR};correction', 1000, 0., 1.)
p1_HFLumiVsMuOrSum  = r.TProfile('p1_HFLumiVsMuOrSum', 'HF <L_{bunch}>;#mu_{OR};<L>/bunch (cm^{-2}s^{-1})', 1000, 0., 1.)
p1_LumiVsMuOrSum    = r.TProfile('p1_LumiVsMuOrSum', 'BCM1F <L_{bunch}>;#mu_{OR};<L>/bunch (cm^{-2}s^{-1})', 1000, 0., 1.)

p1_HFLumiVsMuOr = {}
p1_LumiVsMuOr   = {}
for fill in fillList:
    p1_HFLumiVsMuOr[fill]   = r.TProfile('p1_HFLumiVsMuOr_Fill{}'.format(fill), 'HF <L_{bunch}> (fill '+str(fill)+');#mu_{OR};<L> [cm^{-2}s^{-1}]', 1000, 0., 1.)
    p1_LumiVsMuOr[fill]     = r.TProfile('p1_LumiVsMuOr_Fill{}'.format(fill), 'BCM1F <L_{bunch}> (fill '+str(fill)+');#mu_{OR};<L> [cm^{-2}s^{-1}]', 1000, 0., 1.)

#outFile.cd('PixelLumi')
#p1_PixelLumi        = r.TProfile('p1_PixelLumi', 'Pixel lumi; time (min); 10^{30} cm^{-2}s^{-1}', pixelRange, pixelStart, pixelEnd)

time    = 10.
bunches = 0.
orbit   = 11246.

current = tuple([0., 0.])

### Fill histograms ###
for row in bcmData:
    bunches     = float(row['BUNCHES'])
    #bunches     = float(row['COLLIDING_BUNCHES'])
    fill        = float(row['FILL'])
    beamMode    = row['BEAMMODE']
    start       = time #r.TDatime(row['DT']).Convert()
    #stop        = r.TDatime(row['DT']).Convert()
    #start       = r.TDatime(row['DT_FROM']).Convert()
    #stop        = r.TDatime(row['DT_UNTIL']).Convert()

    if bunches == 0 or beamMode != 'STABLE_BEAMS' or fill not in fillList: continue
    time        += 1

    ### RATES ###
    #lutOrRate   = float(row['LUT_OR_AVG'])
    #orRate      = float(row['ORHZ_AVG']) 
    #xor1Rate    = float(row['XOR1HZ_AVG']) 
    #xor2Rate    = float(row['XOR2HZ_AVG']) 
    #andRate     = float(row['ANDHZ_AVG']) 

    lutOrRate   = float(row['LUT_ORHZ'])
    orRate      = float(row['ORHZ']) 
    xor1Rate    = float(row['XOR1HZ']) 
    xor2Rate    = float(row['XOR2HZ']) 
    andRate     = float(row['ANDHZ']) 

    #print '{}, {}'.format(lutOrRate, orRate)

    ### PROBABILITIES ###
    r0          = 1 - orRate/(11246*bunches)
    rXOr1       = xor1Rate/(11246*bunches) 
    rXOr2       = xor2Rate/(11246*bunches)
    rAnd        = andRate/(11246*bunches)

    ### MU-CORRECTIONS ###
    muOr        = -1*r.TMath.log(r0)
    muXOr1      = r.TMath.log(1 + rXOr1/r0)
    muXOr2      = r.TMath.log(1 + rXOr2/r0)
    muAnd       = -1*r.TMath.log(r0*(1 + rXOr1/r0)*(1 + rXOr2/r0))
    muAnd_Alt   = -1*r.TMath.log(1 - rAnd)

    ### Beam parameters ###
    hfLumi      = 0.9*float(row['HFLUMI'])*1e30 # lower HF lumi by 10%
    
    #if row['I1_COLLIDABLE'] in ['', '0'] or row['I2_COLLIDABLE'] in ['', '0']:
    #    intensity1  = current[0]
    #    intensity2  = current[1]
    #else:
    #    intensity1  = float(row['I1_COLLIDABLE'])
    #    intensity2  = float(row['I2_COLLIDABLE'])

    IxI   = float(row['I1XI2']) #intensity1*intensity2

    ### Lumi calibrations ###

    if hfLumi/bunches > 8e28 and muOr > 0:

        p1_RateOr.Fill(time, orRate)
        p1_RateXOr1.Fill(time, xor1Rate)
        p1_RateXOr2.Fill(time, xor2Rate)
        p1_RateAnd.Fill(time, andRate)

        p1_MuOr.Fill(time, muOr)
        p1_MuXOr1.Fill(time, muXOr1)
        p1_MuXOr2.Fill(time, muXOr2)
        p1_MuAnd.Fill(time, muAnd)
        p1_MuAnd_Alt.Fill(time, muAnd_Alt)

        p1_HFLumi.Fill(start, hfLumi)
        p1_BeamIntensity.Fill(time, IxI)
        p1_BeamIntensityVsMuOr.Fill(IxI, muOr) 
        p1_BeamIntensityVsMuAnd.Fill(IxI, muAnd)

        h2_HFLumiVsMuOr.Fill(hfLumi/bunches, muOr);
        h2_HFLumiVsMuAnd.Fill(hfLumi/bunches, muAnd);
        p1_HFLumiVsMuAnd.Fill(hfLumi/bunches, muAnd);
        p1_HFLumiVsMuAnd_Alt.Fill(hfLumi/bunches, muAnd_Alt);

        lumiOr = bunches*11246*muOr/(2940e-30)
        correction = 1.#0.3408*r.TMath.ATan(9.421*(muOr - 0.4268)) + 1.851

        p1_LumiOr.Fill(start, lumiOr*correction);
        p1_LumiOrByHF.Fill(start, lumiOr/hfLumi);
        p1_RateXOr1ByHF.Fill(start, xor1Rate/hfLumi);
        p1_RateXOr2ByHF.Fill(start, xor2Rate/hfLumi);

        if lutOrRate > 0: p1_RateOrRatio.Fill(start, orRate/lutOrRate)

        #print 'time: {}, BCM1F lumi: {}'.format(row['DT_FROM'], lumiOr*correction) 

        p1_OrCorrections.Fill(muOr, hfLumi/(11245*bunches*muOr/2940e-30))
        p1_HFLumiVsMuOrSum.Fill(muOr, hfLumi);#/bunches);
        p1_LumiVsMuOrSum.Fill(muOr, 11245*muOr/2940e-30);
        p1_LumiVsMuOr[fill].Fill(muOr, lumiOr);
        p1_HFLumiVsMuOr[fill].Fill(muOr, hfLumi);
        
p1_OrCorrections.Rebin(5)

l.format_histogram_time(p1_LumiOr)
l.format_histogram_time(p1_HFLumi)
l.format_histogram_time(p1_LumiOrByHF)
l.format_histogram_time(p1_RateXOr1ByHF)
l.format_histogram_time(p1_RateXOr2ByHF)
l.format_histogram_time(p1_RateOrRatio)

outFile.Write()
outFile.Close()
