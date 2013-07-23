#! /usr/bin/env python
import sys, os, subprocess, math
from array import array
from multiprocessing import Process, Queue
import ROOT as r
import fitVDM as t
import simTools as sim
import lumiTools as l

r.gROOT.SetBatch()
doTest = False

'''
Simple toy MC to produce VdM scan distibutions.  First,
generates a 2-D double Gaussian with some simulated noise.
Then scan slices as in VdM scan to produce points to
feed to the vdmCalibrator.
'''


# Configuration:

# beamTypes variable is of the form (beam 1, beam 2); possible choices are 'SG' = single
# gaussian, 'DG' = double gaussian, and 'DGX' = double gaussian with
# off-diagonal widths. fitTypes variable can be 'singleGaussian', 'doubleGaussian', or
# 'skewGaussian'.

beamTypes   = ('DG','DG') 
fitTypes    = ['singleGaussian', 'doubleGaussian']#, 'skewGaussian']
do2D        = True
nThrows     = int(1e4)
nToys       = 50
nSPs        = 25
scanRange   = (0.3, 0.7) # Should be within [0, 1]
scanPoints  = [scanRange[0] + (scanRange[1] - scanRange[0])*i/float(nSPs) for i in range(nSPs)]

paramSuffix = 'TEST' #sim.get_current_time()
plotPath    = 'plots/{0[0]}_{0[1]}'.format(beamTypes)

if not os.path.exists(plotPath):
    os.system('mkdir -p '+plotPath)


r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)
r.gRandom.SetSeed(0)

# Command line arguments
if len(sys.argv) > 1:
    nThrows = int(sys.argv[1])

# Prepare output
canvas      = r.TCanvas('canvas', 'canvas', 700, 600)

'''
Use Gaussian beam shapes allowing for correlations in
the x and y plane.  For throwing toys we need to keep the 
the distributions well within x,y in (0, 1).  VdM scan can
be simulated by displacing the beam centroid wrt each other
'''

# Define functions for transverse beam distributions
f2_SG   = r.TF2('SingleGauss', 'exp(-(x-[0])**2/(2*[2]**2) - [5]*(x-[0])*(y-[1])/(2*[3]**2) - (y-[1])**2/(2*[4]**2))') 

f2_DG   = r.TF2('DoubleGauss', '[8]*exp(-(x-[4])**2/(2*([0]*[1]/(1+[8]*([1]-1)))**2) \
                                        - (y-[6])**2/(2*([2]*[3]/(1+[8]*([3]-1)))**2)) \
                                        + (1 - [8])*exp(-(x-[5])**2/(2*([0]/(1+[8]*([1]-1)))**2) \
                                        - (y-[7])**2/(2*([2]/(1+[8]*([3]-1)))**2))')

f2_DGX  = r.TF2('DoubleGaussX', '[8]*exp(-(x-[4])**2/(2*[0]**2) + (x-[4])*(y-[6])/(2*[9]**2) - (y-[6])**2/(2*[2]**2)) \
                                        + (1-[8])*exp(-(x-[5])**2/(2*([1]*[0])**2) + (x-[5])*(y-[7])/(2*([10]*[9])**2) - (y-[7])**2/(2*([3]*[2])**2))') 

# Set paramaeters by hand; begin with the beams
# being head on.  NB: since throwing toys from
# a 2D hist only returns values in [0,1], we need
# to contain the beams (or at least, the overlap)
# to well within this range.

##############
### BEAM 1 ###
##############

if beamTypes[0] == 'SG':
    f2_Beam1 = r.TF2('f2_Beam1', 'SingleGauss')
    f2_Beam1.SetParNames('x_{0}', 'y_{0}', '#sigma_{x}','#sigma_{xy}','#sigma_{y}', 'xFactor')

    f2_Beam1.SetParameter('#sigma_{x}', 0.04)
    f2_Beam1.SetParameter('#sigma_{xy}', 0.1)
    f2_Beam1.SetParameter('#sigma_{y}', 0.075)
    f2_Beam1.SetParameter('xFactor', -0.) # Allows mediation of strength of cross-term (0 is off, of course)

    f2_Beam1.SetParameter('x_{0}', 0.5)
    f2_Beam1.SetParameter('y_{0}', 0.499)

elif beamTypes[0] == 'DG':
    f2_Beam1 = r.TF2('f2_Beam1', 'DoubleGauss')
    f2_Beam1.SetParNames('#Sigma_{x}', '#sigma_{1,x}/#sigma{2,x}', '#Sigma_{y}', '#sigma_{1,y}/#sigma{2,y}',\
                         'x_{1}', 'x_{2}', 'y_{1}', 'y_{2}', 'Fraction')

    f2_Beam1.SetParameter('#Sigma_{x}', 0.07)
    f2_Beam1.SetParameter('#sigma_{1,x}/#sigma{2,x}', 0.7)
    f2_Beam1.SetParameter('#Sigma_{y}', 0.07)
    f2_Beam1.SetParameter('#sigma_{1,y}/#sigma{2,y}', 0.3)

    f2_Beam1.SetParameter('x_{1}', 0.5)
    f2_Beam1.SetParameter('x_{2}', 0.5)
    f2_Beam1.SetParameter('y_{1}', 0.5)
    f2_Beam1.SetParameter('y_{2}', 0.5)
    f2_Beam1.SetParameter('Fraction', 0.7)

elif beamTypes[0] == 'DGX':
    f2_Beam1 = r.TF2('f2_Beam1', 'DoubleGaussX')
    f2_Beam1.SetParNames('#sigma_{x}', 'S_{x}', '#sigma_{y}','S_{y}','x_{1}','x_{2}','y_{1}','y_{2}','Fraction', '#sigma_{xy}', 'S_{xy}')

    f2_Beam1.SetParameter('#sigma_{x}', 0.04)
    f2_Beam1.SetParameter('S_{x}', 0.5)
    f2_Beam1.SetParameter('#sigma_{y}', 0.05)
    f2_Beam1.SetParameter('S_{y}', 0.5)
    f2_Beam1.SetParameter('#sigma_{xy}', 0.05)
    f2_Beam1.SetParameter('S_{xy}', 1.)

    f2_Beam1.SetParameter('x_{1}', 0.5)
    f2_Beam1.SetParameter('x_{2}', 0.5)
    f2_Beam1.SetParameter('y_{1}', 0.5)
    f2_Beam1.SetParameter('y_{2}', 0.5)
    f2_Beam1.SetParameter('Fraction', 1.)
    

##############
### BEAM 2 ###
##############

if beamTypes[1] == 'SG':
    f2_Beam2 = r.TF2('f2_Beam2', 'SingleGauss') 
    f2_Beam2.SetParNames('x_{0}', 'y_{0}', '#sigma_{x}','#sigma_{xy}','#sigma_{y}', 'xFactor')

    f2_Beam2.SetParameter('#sigma_{x}', 0.05)
    f2_Beam2.SetParameter('#sigma_{xy}', 0.06)
    f2_Beam2.SetParameter('#sigma_{y}', 0.065)
    f2_Beam2.SetParameter('xFactor', 0.) # Allows mediation of strength of cross-term (0 is off, of course)

    f2_Beam2.SetParameter('x_{0}', 0.5)
    f2_Beam2.SetParameter('y_{0}', 0.501)


elif beamTypes[1] == 'DG':
    f2_Beam2 = r.TF2('f2_Beam2', 'DoubleGauss')
    f2_Beam2.SetParNames('#Sigma_{x}', '#sigma_{1,x}/#sigma{2,x}', '#Sigma_{y}', '#sigma_{1,y}/#sigma{2,y}', \
                         'x_{1}', 'x_{2}', 'y_{1}', 'y_{2}', 'Fraction')

    f2_Beam2.SetParameter('#Sigma_{x}', 0.08)
    f2_Beam2.SetParameter('#sigma_{1,x}/#sigma{2,x}', 0.75)
    f2_Beam2.SetParameter('#Sigma_{y}', 0.08)
    f2_Beam2.SetParameter('#sigma_{1,y}/#sigma{2,y}', 0.65)

    f2_Beam2.SetParameter('x_{1}', 0.5)
    f2_Beam2.SetParameter('x_{2}', 0.5)
    f2_Beam2.SetParameter('y_{1}', 0.5)
    f2_Beam2.SetParameter('y_{2}', 0.5)
    f2_Beam2.SetParameter('Fraction', 0.7)

elif beamTypes[1] == 'DGX':
    f2_Beam2 = r.TF2('f2_Beam2', 'DoubleGaussX')
    f2_Beam2.SetParNames('#sigma_{x}', 'S_{x}', '#sigma_{y}','S_{y}','x_{1}','x_{2}','y_{1}','y_{2}','Fraction', '#sigma_{xy}', 'S_{xy}')

    f2_Beam2.SetParameter('#sigma_{x}', 0.05)
    f2_Beam2.SetParameter('S_{x}', 0.5)
    f2_Beam2.SetParameter('#sigma_{y}', 0.04)
    f2_Beam2.SetParameter('S_{y}', 0.5)
    f2_Beam2.SetParameter('#sigma_{xy}', 0.05)
    f2_Beam2.SetParameter('S_{xy}', 1.)

    f2_Beam2.SetParameter('x_{1}', 0.5)
    f2_Beam2.SetParameter('x_{2}', 0.5)
    f2_Beam2.SetParameter('y_{1}', 0.5)
    f2_Beam2.SetParameter('y_{2}', 0.5)
    f2_Beam2.SetParameter('Fraction', 1.)


### Draw beam profiles
f2_Beam1.SetTitle('Beam 1 PDF;;')
f2_Beam1.Draw('surf1')
canvas.Print(plotPath + '/beam1_shape.pdf')

f2_Beam2.SetTitle('Beam 2 PDF;;')
f2_Beam2.Draw('surf1')
canvas.Print(plotPath + '/beam2_shape.pdf')

### "Convolve" beam 1 and 2 functions to get overlap truth
f2_Overlap = r.TF2('f2_Overlap', 'f2_Beam1*f2_Beam2')
f2_Overlap.SetTitle('Beam overlap (head-on);;')
f2_Overlap.Draw('surf1')
canvas.Print('{0}/overlap_shape.pdf'.format(plotPath))


sigma = {}
if beamTypes[0] == 'SG':
    sigX = [f2_Beam1.GetParameter('#sigma_{x}'),f2_Beam1.GetParameter('#sigma_{x}')]
    sigY = [f2_Beam1.GetParameter('#sigma_{y}'),f2_Beam1.GetParameter('#sigma_{y}')]
    sigma['X'] = sigX[0]
    sigma['Y'] = sigY[0]


### Carry out simulation
print '\nStarting VdM simulation with {0} throws per scan point'.format(nThrows)

### multiprocessing test area ###
if doTest:
    def mc_test(queue, simbot, beam1, beam2, nThrows):
        queue.put(simbot.single_scan_simulator(beam1, beam2, nThrows))

    resultQueue = Queue()
    jobs        = [Process(target = mc_test, args = (resultQueue, simulator, f2_Beam1, f2_Beam2, nThrows,)) for i in range(nToys)]
    for job in jobs: job.start()
    for job in jobs: job.join()
    results     = [resultQueue.get() for i in range(nToys)]

    for result in results:
        print result[1]

    exit()
##################################

simulator = sim.SimTools(beamTypes, scanPoints, canvas, plotPath, paramSuffix)

#truth, rates, sigRates, sigDelta, beamSpot, sigBeamSpot, beamWidth, sigBeamWidth = simulator.single_scan_simulator(f2_Beam1, f2_Beam2, nThrows)

results = []
for i in range(nToys):
    results.append(simulator.single_scan_simulator(f2_Beam1, f2_Beam2, nThrows))

print 'Scan finished!\n'

### Make directory structure ###
if not os.path.exists('{0}/fits'.format(plotPath)):
    os.makedirs('{0}/fits'.format(plotPath))
else:
    if os.listdir('{0}/fits'.format(plotPath)):
        subprocess.call('rm {0}/fits/*pdf'.format(plotPath), shell = True)

if not os.path.exists('{0}/scanPoints'.format(plotPath)):
    os.makedirs('{0}/scanPoints'.format(plotPath))
else:
    if os.listdir('{0}/fits'.format(plotPath)):
        subprocess.call('rm {0}/scanPoints/*pdf'.format(plotPath), shell = True)

if not os.path.exists('{0}/biasPlots'.format(plotPath)):
    os.makedirs('{0}/biasPlots'.format(plotPath))
else:
    if os.listdir('{0}/biasPlots'.format(plotPath)):
        subprocess.call('rm {0}/biasPlots/*pdf'.format(plotPath), shell = True)


biases = []

# Print Carry out fits
for i in range(nToys):

    (truth, rates, sigRates, sigDelta, beamSpot, sigBeamSpot, beamWidth, sigBeamWidth) = results[i] 

    fitGraphs   = {}
    fitResult   = {'X':{}, 'Y':{}}

    diffPoints = [x - scanPoints[-(j+1)] for j,x in enumerate(scanPoints)]

    # Do 1D VDM fit 
    for plane in ['X', 'Y']:

        offset = 0.5

        # Fit simulated rates
        g_fit = r.TGraphErrors(nSPs, array('f', diffPoints), array('f', rates[plane]), array('f', sigDelta[plane]), array('f', sigRates[plane]))
        g_fit.SetTitle('VdM scan ' + plane + ' sim;#Delta' + plane + ' (a.u.);')
        g_fit.Draw('AP')

        canvas.Print(plotPath + '/scanPoints/1D_{0}_{1}_{2}.pdf'.format(plane, paramSuffix, i+1))

        fitGraphs[plane] = g_fit

        for fitType in fitTypes:
            fitResult[plane][fitType] = t.Fit(g_fit.Clone(), fitType, True, 'SIM', '---', '{0}/fits/fit1D_{1}_{2}_{3}_{4}'.format(plotPath, plane, fitType, paramSuffix, i+1))

        g_fit.Clear()
        canvas.Clear()

        simulator.draw_beamspot_plots([diffPoints, sigDelta[plane]], [beamSpot[plane], sigBeamSpot[plane]], [beamWidth[plane], sigBeamWidth[plane]], plane)
        
    # Do 2D VDM fit
    graph2D = []
    if do2D and 'doubleGaussian' in fitTypes:
        graph2D.append(r.TGraph2DErrors(2*nSPs, array('d', diffPoints + [-0.00001 for n in range(nSPs)]),
                                           array('d', [0.00001 for n in range(nSPs)] + diffPoints), 
                                           array('d', rates['X'] + rates['Y']), array('d', 2*sigDelta['X']), 
                                           array('d', 2*sigDelta['Y']), array('d', sigRates['X'] + sigRates['Y'])))

        graph2D[0].SetMarkerStyle(21)
        graph2D[0].SetTitle('Simulated scan points;#Delta X;#Delta Y')
        graph2D[0].Draw('P0')

        canvas.Print(plotPath + '/scanPoints/2D_{0}_{1}.pdf'.format(paramSuffix, i+1))

        fitResult2D = t.Fit2D(graph2D[0].Clone(), fitGraphs['X'], fitGraphs['Y'], True, 'SIM', \
                              fitResult['X']['doubleGaussian'], fitResult['Y']['doubleGaussian'], \
                              '{0}/fits/fit2D_{2}_{3}'.format(plotPath, '', paramSuffix, i+1))

        fitResult['X']['2D'] = [fitResult2D[0][0]*1000, fitResult2D[2]]
        fitResult['Y']['2D'] = [fitResult2D[0][2]*1000, fitResult2D[3]]

    biases.append(simulator.draw_bias_plots(fitResult, [rates, sigRates], truth, fitTypes + ['2D'], i))

styles = {'singleGaussian':[r.kRed, 23], 'doubleGaussian':[r.kCyan+2, 22], 'skewGaussian':[r.kGreen-2, 24], '2D':[r.kViolet, 25]}

if do2D:
    fitTypes.append('2D')

for plane in ['X', 'Y']:
    biasHists = {}
    for fitType in fitTypes:
        hist = r.TH1F('h_bias_{0}'.format(fitType), 'fit biases (scan {0});bias;nToys'.format(plane), 50, -0.1, 0.1)
        hist.SetFillColor(0)
        hist.SetFillStyle(0)
        hist.SetLineColor(styles[fitType][0])
        hist.SetLineWidth(2)

        biasHists[fitType] = hist

    hMax = 0
    for bias in biases:
        for fitType in fitTypes:
            biasHists[fitType].Fill(bias[plane][fitType])
            
            histMax = biasHists[fitType].GetMaximum()
            if hMax < histMax:
                hMax = histMax
        
    legend = r.TLegend(0.6,0.78,0.95,0.94)
    legend.SetFillColor(0)
    legend.SetTextSize(0.045)

    for i,fitType in enumerate(fitTypes):
        legend.AddEntry(biasHists[fitType], '{0}'.format(fitType))
        if i is 0:
            biasHists[fitType].SetMaximum(1.1*hMax)
            biasHists[fitType].Draw('')
        else:
            biasHists[fitType].Draw('SAME')
    legend.Draw()

    canvas.Print('{0}/bias_{1}.pdf'.format(plotPath, plane))


if os.listdir('{0}/fits'.format(plotPath)):
    for fitType in fitTypes:
        if fitType == '2D':
            subprocess.call('pdftk {0}/fits/fit2D*.pdf output {0}/fit2D_{1}.pdf'.format(plotPath, paramSuffix), shell = True)
            subprocess.call('pdftk {0}/scanPoints/2D*.pdf output {0}/scanPoints2D_{1}.pdf'.format(plotPath, paramSuffix), shell = True)
        else:
            subprocess.call('pdftk {0}/fits/fit1D*{1}*.pdf output {0}/fit1D_{1}.pdf'.format(plotPath, fitType), shell = True)
            subprocess.call('pdftk {0}/scanPoints/1D*.pdf output {0}/scanPoints1D_{1}.pdf'.format(plotPath, paramSuffix), shell = True)

        subprocess.call('pdftk {0}/biasPlots/X_*.pdf output {0}/biasPlots_X.pdf'.format(plotPath, paramSuffix), shell = True)
        subprocess.call('pdftk {0}/biasPlots/Y_*.pdf output {0}/biasPlots_Y.pdf'.format(plotPath, paramSuffix), shell = True)


