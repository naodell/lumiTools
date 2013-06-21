#! /usr/bin/env python
import sys, math, csv, json, pprint
from array import array
import ROOT as r
import fitVDM as t
import lumiTools as l
import calibrationTools as c

### Configuration ###

orbit       = 11246.
bunchHeader = ['DT']+range(1, 3565)
savePath    = 'plots'
#savePath    = '~/work/plots/bcm1f_lumi/VDM'


corrections = {'lengthX': 1, 'lengthY': 1, 'intensity1': 1, 'intensity2':1}
algos       = ['OR']#, 'XOR1', 'XOR2', 'AND']
#fills       = [2855, 2856]
fills       = [3316]

bcids = [1, 41, 81, 121, 161, 201, 241,
        721, 761, 801, 841, 881, 921, 961,
        1581, 1621, 1661, 1701, 1741, 1781, 1821,
        #1441, 1481, 1521, 1561, 1601, 1641, 1681, # April VdM (I think...)
        2161, 2201, 2241, 2281, 2321, 2361, 2401,
        #2881, 2921, 2961, 3001, 3041, 3081, 3121]#, 'Total'] # July VdM
        2881]

#bcids = [1, 'Total']

scanResults = []

#-----------#
# Fill 2855 #
#-----------#

if 2855 in fills:

### Scan configuration and corrections ###

    scans = ['1', '2'] #, '3'] 

### Input data files ###

    bcmData = [
        l.read_csv_file('data/July/VDMScanXY1.csv'),
        l.read_csv_file('data/July/VDMScanXY2.csv'),
        l.read_csv_file('data/July/VDMScanXY3.csv')
        ]

    bunchData = [
        l.read_csv_file('data/July/FBCT_2855_Timber/TIMBER_DATA_Fill_2855_BCTFR_A6R4_B1.csv', bunchHeader),
        l.read_csv_file('data/July/FBCT_2855_Timber/TIMBER_DATA_Fill_2855_BCTFR_A6R4_B2.csv', bunchHeader)
        ]

    pltData = r.TFile('histograms/vdmHistograms_2011-07-19_2855.root', 'OPEN')

    ### Initialize calibrator for fill 2855 ###

    calibrator = c.VDMCalibrator(2855, bcmData, bunchData, corrections, {'PLT':1342721696 - 1798, 'intensity':1342685132}, algos, dataFormat = 'brm')
    scanResults.extend(calibrator.multi_scan(scans, bcids, pltData, False, savePath)[0])

    g_Sig = l.draw_summary_plots(bcids[:35], 'OR', scanResults, '2855', savePath+'/fill2855')
    #l.overlay_widths(g_Sig, bcids[:35], True, True, savePath)


#-----------#
# Fill 2856 #
#-----------#

if 2856 in fills:

    ### Scan configuration and corrections ###

    scans = ['1', '2'] 

    ### Input data files ###

    bcmData = [
        l.read_csv_file('data/July/VDMScanXY4.csv')#,
        #l.read_csv_file('data/July/VDMScanXY5.csv')
        ]

    bunchData = [
        l.read_csv_file('data/July/FBCT_2856_Timber/TIMBER_DATA_Fill_2856_BCTFR_A6R4_B1.csv', bunchHeader),
        l.read_csv_file('data/July/FBCT_2856_Timber/TIMBER_DATA_Fill_2856_BCTFR_A6R4_B2.csv', bunchHeader)
        ]

    pltData = r.TFile('histograms/vdmHistograms_2011-07-19_2856.root', 'OPEN')

    ### Initialize calibrator for fill 2856 ###
    calibrator = c.VDMCalibrator(2856, bcmData, bunchData, corrections, {'PLT':1342724785 - 6310, 'intensity':1342713881}, algos, dataFormat = 'brm')
    scanResults.extend(calibrator.multi_scan(scans, bcids, pltData, savePath)[0])

    g_Sig = l.draw_summary_plots(bcids[:35], 'OR', scanResults, '2856', savePath+'/fill2856')

#-----------#
# Fill 3316 #
#-----------#

if 3316 in fills:

    ### Scan configuration and corrections ###
    scans = ['1', '2']#, '3', '4', '5'] 

    ### Input data files ###
    scanData = [
        l.read_csv_file('data/November/ScanPointsNov2012_Scan1.txt'),
        l.read_csv_file('data/November/ScanPointsNov2012_Scan2.txt'),
        l.read_csv_file('data/November/ScanPointsNov2012_Scan3.txt'),
        l.read_csv_file('data/November/ScanPointsNov2012_Scan4.txt'),
        l.read_csv_file('data/November/ScanPointsNov2012_Scan5.txt')
        ]

    bunchData = [
        l.read_csv_file('data/November/TIMBER_DATA_Fill_3316_BCTFR_A6R4_B1.csv', bunchHeader),
        l.read_csv_file('data/November/TIMBER_DATA_Fill_3316_BCTFR_A6R4_B2.csv', bunchHeader)
        ]

    pltData = r.TFile('histograms/vdmHistograms_2012-11-23_3316.root', 'OPEN')

    ### Initialize calibrator for fill 2856 ###
    calibrator = c.VDMCalibrator(3316, scanData, bunchData, corrections, {'PLT':1353700837 - 1218, 'intensity':1353692441}, algos, dataFormat = 'lumi')
    scanResults.extend(calibrator.multi_scan(scans, bcids, pltData, True, savePath)[0])

    g_Sig = l.draw_summary_plots(bcids[:35], 'OR', scanResults, '3316', savePath+'/fill3316')
    l.overlay_widths((g_Sig[0][0], g_Sig[1][0]), bcids[:35], True, True, savePath, '1')
    l.overlay_widths((g_Sig[0][1], g_Sig[1][1]), bcids[:35], True, True, savePath, '2')

