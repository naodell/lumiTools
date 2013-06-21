#! /usr/bin/env python
import sys, os, subprocess, math
from array import array
import ROOT as r
import fitVDM as t
import simTools as sim
import lumiTools as l

r.gROOT.SetBatch()

'''
Simple toy MC to produce VdM scan distibutions.  First,
generates a 2-D double Gaussian with some simulated noise.
Then scan slices as in VdM scan to produce points to
feed to the vdmCalibrator.
'''


# Configuration
beamType    = 'DG'
fitTypes    = ['singleGaussian', 'doubleGaussian']#, 'skewGaussian']
do2D        = True
nToys       = int(1e4)
nSPs        = 25
scanPoints  = [0.20 + 0.025*i for i in range(nSPs)]

paramSuffix = 'TEST' #sim.get_current_time()
plotPath    = 'plots/' + beamType

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)
r.gRandom.SetSeed(0)

# Command line arguments
if len(sys.argv) > 1:
    nToys = int(sys.argv[1])

# Prepare output
canvas      = r.TCanvas('canvas', 'canvas', 700, 600)

'''
Use Gaussian beam shapes allowing for correlations in
the x and y plane.  For throwing toys we need to keep the 
the distributions well within x,y in (0, 1).  VdM scan can
be simulated by displacing the beam centroid wrt each other
'''

# Define functions for transverse beam distributions
f2_SingleGaus = r.TF2('SingleGauss', 'exp(-(x-[0])**2/(2*[2]**2) - [5]*(x-[0])*(y-[1])/(2*[3]**2) - (y-[1])**2/(2*[4]**2))') 

if beamType == 'SG':
    f2_Beam1 = r.TF2('f2_Beam1', 'SingleGauss')
    f2_Beam1.SetParNames('x_{0}', 'y_{0}', '#sigma_{x}','#sigma_{xy}','#sigma_{y}', 'xFactor')

    f2_Beam2 = r.TF2('f2_Beam2', 'SingleGauss') 
    f2_Beam2.SetParNames('x_{0}', 'y_{0}', '#sigma_{x}','#sigma_{xy}','#sigma_{y}', 'xFactor')

    # Set paramaeters by hand; begin with the beams
    # being head on.  NB: since throwing toys from
    # a 2D hist only returns values in [0,1], we need
    # to contain the beams (or at least, the overlap)
    # to well within this range.

    f2_Beam1.SetParameter('x_{0}', 0.5)
    f2_Beam1.SetParameter('y_{0}', 0.49)
    f2_Beam1.SetParameter('#sigma_{x}', 0.07)
    f2_Beam1.SetParameter('#sigma_{xy}', 0.1)
    f2_Beam1.SetParameter('#sigma_{y}', 0.05)
    f2_Beam1.SetParameter('xFactor', -0.)

    f2_Beam2.SetParameter('x_{0}', 0.5)
    f2_Beam2.SetParameter('y_{0}', 0.51)
    f2_Beam2.SetParameter('#sigma_{x}', 0.05)
    f2_Beam2.SetParameter('#sigma_{xy}', 0.06)
    f2_Beam2.SetParameter('#sigma_{y}', 0.05)
    f2_Beam2.SetParameter('xFactor', 0.)


#f2_DoubleGaus = r.TF2('DoubleGauss', '[8]*exp(-(x-[4])**2/(2*[0]**2) - (y-[6])**2/(2*[2]**2)) \
#                                        + (1-[8])*exp(-(x-[5])**2/(2*([1]*[0])**2) - (y-[7])**2/(2*([3]*[2])**2))') 

f2_DoubleGaus = r.TF2('DoubleGauss', '[8]*exp(-(x-[4])**2/(2*([0]*[1]/(1+[8]*([1]-1)))**2) - (y-[6])**2/(2*([2]*[3]/(1+[8]*([3]-1)))**2)) \
                                        + (1 - [8])*exp(-(x-[5])**2/(2*([0]/(1+[8]*([1]-1)))**2) - (y-[7])**2/(2*([2]/(1+[8]*([3]-1)))**2))')

if beamType == 'DG':

    f2_Beam1 = r.TF2('f2_Beam1', 'DoubleGauss')
    f2_Beam1.SetParNames('#Sigma_{x}', 'sigma_{1,x}/sigma{2,x}', '#Sigma_{y}', 'sigma_{1,y}/sigma{2,y}', 'x_{1}', 'x_{2}', 'y_{1}', 'y_{2}', 'Fraction')

    f2_Beam1.SetParameter('#Sigma_{x}', 0.07)
    f2_Beam1.SetParameter('sigma_{1,x}/sigma{2,x}', 0.5)
    f2_Beam1.SetParameter('#Sigma_{y}', 0.07)
    f2_Beam1.SetParameter('sigma_{1,y}/sigma{2,y}', 0.5)

    f2_Beam1.SetParameter('x_{1}', 0.5)
    f2_Beam1.SetParameter('x_{2}', 0.5)
    f2_Beam1.SetParameter('y_{1}', 0.5)
    f2_Beam1.SetParameter('y_{2}', 0.5)
    f2_Beam1.SetParameter('Fraction', 0.7)

    f2_Beam2 = r.TF2('f2_Beam2', 'DoubleGauss')
    f2_Beam2.SetParNames('#Sigma_{x}', 'sigma_{1,x}/sigma{2,x}', '#Sigma_{y}', 'sigma_{1,y}/sigma{2,y}', 'x_{1}', 'x_{2}', 'y_{1}', 'y_{2}', 'Fraction')

    f2_Beam2.SetParameter('#Sigma_{x}', 0.05)
    f2_Beam2.SetParameter('sigma_{1,x}/sigma{2,x}', 0.6)
    f2_Beam2.SetParameter('#Sigma_{y}', 0.05)
    f2_Beam2.SetParameter('sigma_{1,y}/sigma{2,y}', 0.6)

    f2_Beam2.SetParameter('x_{1}', 0.5)
    f2_Beam2.SetParameter('x_{2}', 0.5)
    f2_Beam2.SetParameter('y_{1}', 0.5)
    f2_Beam2.SetParameter('y_{2}', 0.5)
    f2_Beam2.SetParameter('Fraction', 0.7)


f2_DoubleGaus = r.TF2('DoubleGaussX', '[8]*exp(-(x-[4])**2/(2*[0]**2) + (x-[4])*(y-[6])/(2*[9]**2) - (y-[6])**2/(2*[2]**2)) \
                                        + (1-[8])*exp(-(x-[5])**2/(2*([1]*[0])**2) + (x-[5])*(y-[7])/(2*([10]*[9])**2) - (y-[7])**2/(2*([3]*[2])**2))') 

if beamType == 'DGX':

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

sigma = {}
if beamType == 'SG':
    sigX = [f2_Beam1.GetParameter('#sigma_{x}'),f2_Beam1.GetParameter('#sigma_{x}')]
    sigY = [f2_Beam1.GetParameter('#sigma_{y}'),f2_Beam1.GetParameter('#sigma_{y}')]
    sigma['X'] = sigX[0]
    sigma['Y'] = sigY[0]

    #sigma['X'] = (sigX[0]*sigX[1]/math.sqrt(pow(sigX[0], 2) + pow(sigX[1], 2)))
    #sigma['Y'] = (sigY[0]*sigY[1]/math.sqrt(pow(sigY[0], 2) + pow(sigY[1], 2)))


### Carry out simulation
print '\nStarting VdM simulation with {0} throws per scan point'.format(nToys)
simulator = sim.SimTools(beamType, scanPoints, canvas, plotPath, paramSuffix)
rates, sigRates, sigDelta, beamSpot, sigBeamSpot, beamWidth = simulator.single_scan_simulator(f2_Beam1, f2_Beam2, nToys)
print 'Scan finished!\n'


### Convolve beam 1 and 2 functions to get overlap truth
f2_Overlap = r.TF2('f2_Overlap', 'f2_Beam1*f2_Beam2')
f2_Overlap.SetTitle('Beam overlap (head-on);;')
f2_Overlap.Draw('surf1')
canvas.Print('{0}/overlap_shape.pdf'.format(plotPath))

if not os.path.exists('{0}/fits'.format(plotPath)):
    os.makedirs('{0}/fits'.format(plotPath))
else:
    if os.listdir('{0}/fits'.format(plotPath)):
        subprocess.call('rm {0}/fits/*pdf'.format(plotPath), shell = True)

truth = simulator.get_vdm_truth(f2_Beam1, f2_Beam2)

# Print rates
fitGraphs   = {}
fitResult   = {'X':{}, 'Y':{}}

# Do 1D VDM fit 
for plane in ['X', 'Y']:

    offset = 0.5

    # Fit simulated rates
    g_fit = r.TGraphErrors(nSPs, array('f', [x - offset for x in scanPoints]), array('f', rates[plane]), array('f', sigDelta[plane]), array('f', sigRates[plane]))
    g_fit.SetTitle('VdM scan ' + plane + ' sim;#Delta' + plane + ' (a.u.);')
    g_fit.Draw('AP')

    canvas.Print(plotPath + '/scanPoints_{0}_{1}.pdf'.format(plane, paramSuffix))

    fitGraphs[plane] = g_fit

    for fitType in fitTypes:
        fitResult[plane][fitType] = t.Fit(g_fit.Clone(), fitType, True, 'SIM', '---', '{0}/fits/fit1D_{1}_{2}_{3}'.format(plotPath, plane, fitType, paramSuffix))

    g_fit.Clear()
    canvas.Clear()

    simulator.draw_beamspot_plots([[x - offset for x in scanPoints], sigDelta[plane]], [beamSpot[plane], sigBeamSpot[plane]], beamWidth[plane], plane)
    
# Do 2D VDM fit
if do2D and 'singleGaussian' in fitTypes:
    graph2D = r.TGraph2DErrors(2*nSPs, array('d', [x - offset for x in scanPoints] + [-0.00001 for n in range(nSPs)]),
                                       array('d', [0.00001 for n in range(nSPs)] + [x - offset for x in scanPoints]), 
                                       array('d', rates['X'] + rates['Y']), array('d', 2*sigDelta['X']), 
                                       array('d', 2*sigDelta['Y']), array('d', sigRates['X'] + sigRates['Y']))

    graph2D.SetMarkerStyle(21)
    graph2D.SetTitle('Simulated scan points;#Delta X;#Delta Y')
    graph2D.Draw('P0')

    canvas.Print(plotPath + '/scanPoints_2D_{0}.pdf'.format(paramSuffix))

    fitResult2D = t.Fit2D(graph2D.Clone(), fitGraphs['X'], fitGraphs['Y'], False, 'SIM', fitResult['X']['singleGaussian'], fitResult['Y']['singleGaussian'], '{0}/fits/fit2D_{1}_{2}'.format(plotPath, '', paramSuffix))

    #graph2D.Clear()

simulator.draw_bias_plots(fitResult, [rates, sigRates], truth, fitTypes)

if os.listdir('{0}/fits'.format(plotPath)):
    for fitType in fitTypes:
            subprocess.call('pdftk {0}/fits/fit1D*{1}*.pdf output {0}/fit1D_{1}.pdf'.format(plotPath, fitType), shell = True)

#    if do2D and 'singleGaussian' in fitTypes:
#        subprocess.call('pdftk {0}/fits/fit2D*.pdf output {0}/fit2D_{1}.pdf'.format(plotPath, paramSuffix), shell = True)
