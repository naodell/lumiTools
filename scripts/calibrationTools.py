#! /usr/bin/env python
import sys, math, csv, json, pprint, subprocess
from array import array
import ROOT as r
import fitVDM as t
import lumiTools as l

class VDMCalibrator():
    '''
    Provides calibration for BCM1F lumi from VDM
    scan data.
    '''
    def __init__(self, fill, scanData, bunchData, corrections, offsets, algos, dataFormat = 'brm'):
        self._fill          = fill
        self._scanData      = scanData
        self._bunchData     = bunchData
        self._corrections   = corrections
        self._offsets       = offsets
        self._algos         = algos
        self._dataFormat    = dataFormat

    def calculate_sigma_mu(self, prob, nBunches, algo):
        '''
        Calculates errors for mu algorithms.
        '''

        r0 = 1 - prob['OR']

        if algo == 'OR':
            sigma    = sigmaR0/r0
            sigmaR0  = math.sqrt((2 - 3*r0 + r0*r0)/nBunches) 
            return sigma

        elif algo in ['XOR1', 'XOR2']:
            rXOR    = prob[algo] 
            sigma   = math.sqrt((r0*rXOR - 4*r0*pow(rXOR, 2) + pow(rXOR, 2))/(nBunches*r0*(r0 + rXOR)))
            return sigma

        elif algo == 'AND':
            rXOR1 = prob['XOR1']
            rXOR2 = prob['XOR2']
            sigma = math.sqrt(abs(rXOR1*rXOR2 + pow(r0, 2)*(1 - rXOR1 - rXOR2) - r0*rXOR1*rXOR2 - pow(r0, 3))/(nBunches*r0*(r0 + rXOR1)*(r0 + rXOR2)))
            return sigma

    def get_scan_arrays_brm(self, inputHistograms, scan, bunch, scanPlane, algos):
        '''
        Returns an array for a single bunch in a single 
        scan.  Array contains mu values versus separation 
        values and their associated errors.
        '''

        muVdm       = l.dict_map(algos, [{} for i in range(len(algos))]) 
        muErr       = l.dict_map(algos, [{} for i in range(len(algos))]) 
        delta       = []
        deltaErr    = []

        preSeparation   = 0.
        dipTime         = 0.
        orbit           = 11246

        for row in self._scanData[int(scan)-1]:

            if row['BUNCHES'] == '' or row['IPSCAN'] != '5' or row['SCAN_FLAG'] == 0: continue

            dipTime     = r.TDatime(row['DT']).Convert(r.kTRUE)
            separation  = float(row['SEP'])
            plane       = row['PLANE']
            time        = dipTime - self._offsets['PLT']

            #print time, dipTime, row['DT']

            if plane is not scanPlane or separation == 0.: continue

            ### Beam parameters ###
            delta.append(separation)

            ### GET RATES ###
            prob = {}

            prob['OR']     = inputHistograms['OR'].GetBinContent(inputHistograms['OR'].FindBin(time*1e3))/orbit
            prob['XOR1']   = inputHistograms['XOR1'].GetBinContent(inputHistograms['XOR1'].FindBin(time*1e3))/orbit
            prob['XOR2']   = inputHistograms['XOR2'].GetBinContent(inputHistograms['XOR2'].FindBin(time*1e3))/orbit
            prob['AND']    = inputHistograms['AND'].GetBinContent(inputHistograms['AND'].FindBin(time*1e3))/orbit
            prob['0']      = 1 - prob['OR'] 

            if prob['OR'] == 0: continue

            #print '{}, {}'.format(prob['OR'], time)

            ### MU-CORRECTIONS ###
            mu = {}
            mu['OR']        = -1*r.TMath.log(prob['0'])
            mu['XOR1']      = r.TMath.log(1 + prob['XOR1']/prob['0'])
            mu['XOR2']      = r.TMath.log(1 + prob['XOR2']/prob['0'])
            mu['AND']       = -1*r.TMath.log(prob['0']*(1 + prob['XOR1']/prob['0'])*(1 + prob['XOR2']/prob['0']))
            mu['AND_ALT']   = -1*r.TMath.log(1 - prob['AND'])


            if bunch is not 'Total':
                beamTime    = int((dipTime - self._offsets['intensity'])/60.)

                #print beamTime, dipTime, self._offsets['intensity']

                intensity1  = 1. #float(self._bunchData[0][beamTime][int(bunch)])*1e-11
                intensity2  = 1. #float(self._bunchData[1][beamTime][int(bunch)])*1e-11
                IxI = intensity1*intensity2

            for algo in algos:
                if prob[algo] == 0: continue

                if separation != preSeparation:
                    muVdm[algo][separation] = 0
                    muErr[algo][separation] = 0

                error = self.calculate_sigma_mu(prob, 11246., algo)

                muVdm[algo][separation] += (mu[algo]/IxI)/pow(error/IxI, 2)
                muErr[algo][separation] += pow(IxI/error, 2)

            preSeparation = separation


        # Prepare arrays with scan data for use in fits

        delta = sorted(list(set(delta)))

        for algo in algos:
            for separation in delta:
                muVdm[algo][separation] = muVdm[algo][separation]/muErr[algo][separation]    
                muErr[algo][separation] = math.sqrt(1/muErr[algo][separation])
                deltaErr.append(1e-7)

        return (muVdm, muErr, delta, deltaErr)

    def get_scan_arrays_lumi(self, inputHistograms, scan, bunch, scanPlane, algos):
        '''
        Uses CSV files provided by lumi group.  The files contain
        times for the beginning and end of each scan separation.
        Returns an array for a single bunch in a single scan.  
        Array contains mu values versus separation values and
        their associated errors.
        '''

        muVdm       = l.dict_map(algos, [{} for i in range(len(algos))]) 
        muErr       = l.dict_map(algos, [{} for i in range(len(algos))]) 
        delta       = []
        deltaErr    = []

        preSeparation   = 0.
        dipTime         = 0.
        orbit           = 11246.

        for row in self._scanData[int(scan)-1]:

            time_i      = int(row['TIME_I']) - self._offsets['PLT']
            time_f      = int(row['TIME_F']) - self._offsets['PLT']
            separation  = float(row['SEP'])
            plane       = row['PLANE']


            if plane is not scanPlane or separation == 0.: continue

            ### Beam parameters ###
            delta.append(separation)

            for time in range(time_i, time_f, 2):

                ### PROBABILITIES ###
                prob = {}

                prob['OR']     = inputHistograms['OR'].GetBinContent(inputHistograms['OR'].FindBin(time*1e3))/(2*orbit)
                prob['XOR1']   = inputHistograms['XOR1'].GetBinContent(inputHistograms['XOR1'].FindBin(time*1e3))/(2*orbit)
                prob['XOR2']   = inputHistograms['XOR2'].GetBinContent(inputHistograms['XOR2'].FindBin(time*1e3))/(2*orbit)
                prob['AND']    = inputHistograms['AND'].GetBinContent(inputHistograms['AND'].FindBin(time*1e3))/(2*orbit)
                prob['0']      = 1 - prob['OR'] 

                if prob['OR'] == 0: continue

                #print '{}, {}'.format(prob['OR'], time)

                ### MU-CORRECTIONS ###
                mu = {}
                mu['OR']        = -1*r.TMath.log(prob['0'])
                mu['XOR1']      = r.TMath.log(1 + prob['XOR1']/prob['0'])
                mu['XOR2']      = r.TMath.log(1 + prob['XOR2']/prob['0'])
                mu['AND']       = -1*r.TMath.log(prob['0']*(1 + prob['XOR1']/prob['0'])*(1 + prob['XOR2']/prob['0']))
                mu['AND_ALT']   = -1*r.TMath.log(1 - prob['AND'])


                if bunch is not 'Total':
                    iRow = int((int(row['TIME_I']) + time - self._offsets['intensity'])/60.)

                    #print iRow, time, int(row['TIME_I']), self._offsets['intensity']

                    intensity1  = float(self._bunchData[0][iRow][int(bunch)])*1e-11
                    intensity2  = float(self._bunchData[1][iRow][int(bunch)])*1e-11
                    IxI = intensity1*intensity2

                for algo in algos:
                    if prob[algo] == 0: continue

                    if separation != preSeparation:
                        muVdm[algo][separation] = 0
                        muErr[algo][separation] = 0

                    error = self.calculate_sigma_mu(prob, 11246., algo)

                    muVdm[algo][separation] += (mu[algo]/IxI)/pow(error/IxI, 2)
                    muErr[algo][separation] += pow(IxI/error, 2)


        # Prepare arrays with scan data for use in fits

        delta = sorted(list(set(delta)))

        for algo in algos:
            for separation in delta:
                muVdm[algo][separation] = muVdm[algo][separation]/muErr[algo][separation]    
                muErr[algo][separation] = math.sqrt(1/muErr[algo][separation])
                deltaErr.append(1e-7)

        return (muVdm, muErr, delta, deltaErr)


    def multi_scan(self, scans, bcids, pltData, doCorrFit, savePath):
        '''
        Analyzes multiple scans
        '''
    
        scanResults     = []
        scanResults2D   = []
        for scan in scans:
            fitResults      = l.dict_map(['OR', 'AND', 'XOR1', 'XOR2'], [{}, {}, {}, {}])
            fitResults2D    = l.dict_map(['OR', 'AND', 'XOR1', 'XOR2'], [{}, {}, {}, {}])
            resultFile      = open(savePath+'/fill'+str(self._fill)+'/scan'+scan+'_fitSummary.txt', 'w')


            for bcid in bcids:
                bunch = str(bcid)

                histograms = {
                    'OR':pltData.GetDirectory('BCID_'+bunch).Get('p1_RateVsTime_OR'), 
                    'AND': pltData.GetDirectory('BCID_'+bunch).Get('p1_RateVsTime_AND'), 
                    'XOR1': pltData.GetDirectory('BCID_'+bunch).Get('p1_RateVsTime_XOR1'), 
                    'XOR2': pltData.GetDirectory('BCID_'+bunch).Get('p1_RateVsTime_XOR2')
                }

                fitResult   = {'X':[], 'Y':[]}
                fitGraphs   = {}
                mu2D        = {'X':[], 'Y':[], 'Z':[]}
                mu2DErr     = {'X':[], 'Y':[], 'Z':[]}

                for plane in ['X', 'Y']:
                    if self._dataFormat == 'brm':
                        fitInput = self.get_scan_arrays_brm(histograms, scan, bunch, plane, self._algos)
                    elif self._dataFormat == 'lumi':
                        fitInput = self.get_scan_arrays_lumi(histograms, scan, bunch, plane, self._algos)

                    muOr  = (
                        [value for (key, value) in sorted(zip(fitInput[0]['OR'].keys(), fitInput[0]['OR'].values()))],
                        [value for (key, value) in sorted(zip(fitInput[1]['OR'].keys(), fitInput[1]['OR'].values()))]
                    )

                    vdmGraph = r.TGraphErrors(len(fitInput[2]), array('f', fitInput[2]), array('f', muOr[0]), array('f', fitInput[3]), array('f', muOr[1]))
                    vdmGraph.SetTitle(plane+'-scan BCID '+bunch+' (BCM1F);#Delta '+plane+' (mm);')
                    fitResult[plane] = t.Fit(vdmGraph, 'doubleGaussian', True, self._fill, '1', savePath+'/fill'+str(self._fill)+'/scan'+scan+'/vdmScan_1D_OR_BCM1F_BCID-'+bunch+'_'+plane)

                    # Collect data for correlated fits
                    fitGraphs[plane] = vdmGraph

                    if plane == 'X':
                        mu2D['X'].extend(fitInput[2])
                        mu2DErr['X'].extend(fitInput[3])
                        mu2D['Y'].extend([0. for n in range(len(fitInput[2]))])
                        mu2DErr['Y'].extend([0. for n in range(len(fitInput[2]))])
                    elif plane == 'Y':
                        mu2D['X'].extend([0. for n in range(len(fitInput[2]))])
                        mu2DErr['X'].extend([0. for n in range(len(fitInput[2]))])
                        mu2D['Y'].extend(fitInput[2])
                        mu2DErr['Y'].extend(fitInput[3])

                    mu2D['Z'].extend(muOr[0])
                    mu2DErr['Z'].extend(muOr[1])


                graph2D = r.TGraph2DErrors(len(mu2D['X']), \
                                    array('d', mu2D['X']), array('d', mu2D['Y']), array('d', mu2D['Z']),\
                                    array('d', mu2DErr['X']), array('d', mu2DErr['Y']), array('d', mu2DErr['Z']))

                fitResult2D = t.Fit2D(graph2D, fitGraphs['X'], fitGraphs['Y'], True, self._fill, fitResult['X'], fitResult['Y'], savePath+'/fill'+str(self._fill)+'/scan'+scan+'/vdmScan_2D_OR_BCM1F_BCID-'+bunch)
                graph2D.Clear()

                sigmaEff = r.TMath.Pi()*fitResult['X'][0]*fitResult['Y'][0]*(fitResult['X'][2] + fitResult['Y'][2])/1e30
                sigmaErr = (r.TMath.Pi()/1e30)*math.sqrt((pow(fitResult['X'][0]*fitResult['Y'][1], 2) \
                            + pow(fitResult['X'][1]*fitResult['Y'][0], 2))*pow(fitResult['X'][2] + fitResult['Y'][2], 2) \
                            + pow(fitResult['X'][0]*fitResult['X'][1], 2)*(fitResult['X'][3]*fitResult['X'][3] + fitResult['Y'][3]*fitResult['Y'][3]))

                fitResults['OR'][bunch]     = (fitResult['X'], fitResult['Y'], (sigmaEff*1e30, sigmaErr*1e30))
                fitResults2D['OR'][bunch]   = (fitResult2D, (sigmaEff*1e30, sigmaErr*1e30))

            subprocess.call('pdftk ' + savePath + '/fill' + str(self._fill) + '/scan' + scan + '/vdmScan_1D*.pdf output ' + savePath + '/fill' + str(self._fill) + '/scan' + scan + '_1D.pdf', shell = True)
            subprocess.call('pdftk ' + savePath + '/fill' + str(self._fill) + '/scan' + scan + '/vdmScan_2D*.pdf output ' + savePath + '/fill' + str(self._fill) + '/scan' + scan + '_2D.pdf', shell = True)
            l.print_vdm_tables('OR', bcids, fitResults, resultFile, False, False) 

            scanResults.append(fitResults)
            scanResults2D.append(fitResults2D)

        return scanResults, scanResults2D
