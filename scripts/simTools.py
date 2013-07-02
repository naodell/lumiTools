import os, sys, math, time, datetime, subprocess, multiprocessing
from array import array
import lumiTools as l
import ROOT as r

graphStyles = {'singleGaussian':[r.kRed, 23], 'doubleGaussian':[r.kCyan+2, 22], 'skewGaussian':[r.kGreen-2, 24], '2D':[r.kViolet, 25]}

def get_current_time():
    ''' 
    Returns a string of the current time with
    the format  
    '''

    now = datetime.datetime.now()
    currentTime = '{0:02d}{1:02d}{2:02d}_{3:02d}{4:02d}{5:02d}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)

    return currentTime

def set_graph_style(graph, title = ';x;y',  color = 1, width = 2, mStyle = 21, mSize = 0.8):
    
    graph.SetTitle(title)
    graph.SetFillColor(0)
    graph.SetFillStyle(0)
    graph.SetLineColor(color)
    graph.SetLineWidth(width)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(mStyle)
    graph.SetMarkerSize(mSize)


class SimTools():

    def __init__(self, beamType, scanPoints, canvas, plotPath = '', suffix = ''):
        
        self._suffix        = suffix
        self._plotPath      = plotPath
        self._beamType      = beamType
        self._scanPoints    = scanPoints
        self._nScanPoints   = len(scanPoints)
        self._canvas        = canvas


    def single_scan_simulator(self, f2_Beam1, f2_Beam2, nToys):
        '''
        Carries out simulation of one complete scan
        '''

        self._canvas.SetLogz()

        h2_Beam1    = r.TH2D('h2_Beam1', 'beam 1;x;y', 100, 0., 1., 100, 0., 1.)
        h2_Beam2    = r.TH2D('h2_Beam2', 'beam 2;x;y', 100, 0., 1., 100, 0., 1.)
        h2_Overlap  = r.TH2D('h2_Overlap', 'beam1/2 overlap;x;y', 100, 0., 1., 100, 0., 1.)

        # Set histogram styles 
        h2_Beam1.SetMarkerColor(r.kBlue)
        h2_Beam1.SetFillColor(r.kBlue)
        h2_Beam1.SetLineColor(r.kBlue)
        h2_Beam2.SetMarkerColor(r.kRed)
        h2_Beam2.SetFillColor(r.kRed)
        h2_Beam2.SetLineColor(r.kRed)
        h2_Overlap.SetMarkerColor(r.kViolet)
        h2_Overlap.SetLineColor(r.kViolet)

        # Outputs for fits
        rates       = {'X':[], 'Y':[]}
        sigRates    = {'X':[], 'Y':[]}
        sigDelta    = {'X':[], 'Y':[]}
        beamWidth   = {'X':{'X':[], 'Y':[]}, 'Y':{'X':[], 'Y':[]}}
        beamSpot    = {'X':{'X':[], 'Y':[]}, 'Y':{'X':[], 'Y':[]}}
        sigBeamSpot = {'X':{'X':[], 'Y':[]}, 'Y':{'X':[], 'Y':[]}}

        for plane in ['X', 'Y']:

            print 'Doing scan in ' + plane + ' direction...'
            #self._canvas.Print(self._plotPath + '/vdmScan_MC_' + plane + '.pdf[')

            if os.path.exists(self._plotPath + '/vdmScan_MC_' + plane + '.gif'):
                subprocess.call('rm ' + self._plotPath + '/vdmScan_MC_' + plane + '.gif', shell=True)

            # center x offset for y scan and vice versa
            if plane == 'Y':
                if self._beamType == 'SG':
                    f2_Beam2.SetParameter('x_{0}', 0.5)
                elif self._beamType in ['DG', 'DGX']:
                    f2_Beam2.SetParameter('x_{1}', 0.5)
                    f2_Beam2.SetParameter('x_{2}', 0.5)

            if plane == 'X':
                if self._beamType == 'SG':
                    f2_Beam2.SetParameter('y_{0}', 0.5)
                elif self._beamType in ['DG', 'DGX']:
                    f2_Beam2.SetParameter('y_{1}', 0.5)
                    f2_Beam2.SetParameter('y_{2}', 0.5)
            
            for index,point in enumerate(self._scanPoints):

                #print 'Simulating scan point {} of {}'.format(index + 1, len(self._scanPoints)) 

                # Scan beam 2 w.r.t. beam 1
                if plane == 'X':
                    if self._beamType == 'SG':
                        f2_Beam2.SetParameter('x_{0}', point)
                    elif self._beamType in ['DG', 'DGX']:
                        f2_Beam2.SetParameter('x_{1}', point)
                        f2_Beam2.SetParameter('x_{2}', point)
                elif plane == 'Y':
                    if self._beamType == 'SG':
                        f2_Beam2.SetParameter('y_{0}', point)
                    elif self._beamType in ['DG', 'DGX']:
                        f2_Beam2.SetParameter('y_{1}', point)
                        f2_Beam2.SetParameter('y_{2}', point)


                # Increase beam width as a function of beam
                # separation (TO DO)
                #f2_Beam1.SetParameter('#sigma_{x}', 0.)

                # For throwing toys from input distribution
                x_coord = r.Double(0)
                y_coord = r.Double(0)

                for i in range(nToys):
                    f2_Beam1.GetRandom2(x_coord, y_coord)
                    h2_Beam1.Fill(x_coord, y_coord, 1./nToys)

                    f2_Beam2.GetRandom2(x_coord, y_coord)
                    h2_Beam2.Fill(x_coord, y_coord, 1./nToys)

                h2_Overlap.Multiply(h2_Beam1, h2_Beam2)

                h2_Overlap.SetMaximum(5*nToys)
                h2_Overlap.Draw('COLZ')
                h2_Beam1.Draw('HIST SAME')
                h2_Beam2.Draw('HIST SAME')

                #self._canvas.Print(self._plotPath + '/vdmScan_MC_' + plane + '.pdf')
                self._canvas.Print(self._plotPath + '/vdmScan_MC_' + plane + '.gif+10')

                # Normalize overlap histogram to represent 
                # probability given by overlap of beam shapes
                beamspot    = [r.Long(0), r.Long(0), r.Long(0)]
                h2_Overlap.GetMaximumBin(beamspot[0], beamspot[1], beamspot[2])
                beamspot[0], beamspot[1] = beamspot[0]/100., beamspot[1]/100.

                mpVal       = h2_Overlap.GetMaximumBin()
                mpAmp       = h2_Overlap.GetBinContent(mpVal)
                olIntegral  = h2_Overlap.Integral()

                # record rates and errors for fits
                rates[plane].append(olIntegral)
                sigRates[plane].append(math.sqrt(0.1*rates[plane][index]/nToys)) 
                sigDelta[plane].append(0.01)

                # beamspot! 
                beamWidth[plane]['X'].append(h2_Overlap.ProjectionX().GetRMS())
                beamWidth[plane]['Y'].append(h2_Overlap.ProjectionY().GetRMS())
                beamSpot[plane]['X'].append(beamspot[0])
                beamSpot[plane]['Y'].append(beamspot[1])
                sigBeamSpot[plane]['X'].append(h2_Overlap.GetRMS(1))
                sigBeamSpot[plane]['Y'].append(h2_Overlap.GetRMS(2))

                h2_Overlap.Reset()
                h2_Beam1.Reset()
                h2_Beam2.Reset()

            #self._canvas.Print(self._plotPath + '/vdmScan_MC_' + plane + '.pdf]')

        self._canvas.SetLogz(0)

        return rates, sigRates, sigDelta, beamSpot, sigBeamSpot, beamWidth

    def get_vdm_truth(self, f2_Beam1, f2_Beam2):

        rates       = {'X':[], 'Y':[]}

        for plane in ['X', 'Y']:

            # center x offset for y scan and vice versa
            if plane == 'Y':
                if self._beamType == 'SG':
                    f2_Beam2.SetParameter('x_{0}', 0.5)
                elif self._beamType in ['DG', 'DGX']:
                    f2_Beam2.SetParameter('x_{1}', 0.5)
                    f2_Beam2.SetParameter('x_{2}', 0.5)

            if plane == 'X':
                if self._beamType == 'SG':
                    f2_Beam2.SetParameter('y_{0}', 0.5)
                elif self._beamType in ['DG', 'DGX']:
                    f2_Beam2.SetParameter('y_{1}', 0.5)
                    f2_Beam2.SetParameter('y_{2}', 0.5)
            
            for index,point in enumerate(self._scanPoints):

                # Scan beam 2 w.r.t. beam 1
                if plane == 'X':
                    if self._beamType == 'SG':
                        f2_Beam2.SetParameter('x_{0}', point)
                    elif self._beamType in ['DG', 'DGX']:
                        f2_Beam2.SetParameter('x_{1}', point)
                        f2_Beam2.SetParameter('x_{2}', point)
                elif plane == 'Y':
                    if self._beamType == 'SG':
                        f2_Beam2.SetParameter('y_{0}', point)
                    elif self._beamType in ['DG', 'DGX']:
                        f2_Beam2.SetParameter('y_{1}', point)
                        f2_Beam2.SetParameter('y_{2}', point)

                # Get beam overlap truth (convolution of beam 1 and 2)
                f2_Overlap = r.TF2('f2_Overlap', 'f2_Beam1*f2_Beam2')
                
                rates[plane].append(f2_Overlap.Integral(0., 1., 0., 1.))

        return rates


    def draw_bias_plots(self, f_fit, simRates, truth, fitTypes):
        

        for plane in ['X', 'Y']:

            # Set offset fo scanpoints; possibly unnecessary
            offset = 0.5

            # Normalize truth to simulated overlap; should be okay given
            # arbitrary normalization.
            truthScale = sum(simRates[0][plane])/sum(truth[plane])

            scanPoints  = [x - offset for x in self._scanPoints]
            biasRates   = []
            biasFits    = dict(zip(fitTypes, [[] for i in fitTypes]))
            fitRates    = dict(zip(fitTypes, [[] for i in fitTypes]))

            for i,rate in enumerate(simRates[0][plane]):
                truth[plane][i] = truthScale*truth[plane][i]

                biasRates.append((truth[plane][i] - simRates[0][plane][i])/truth[plane][i])
                #print biasRates[i]

                for fitType in fitTypes:
                    fitRate = f_fit[plane][fitType][-1].Eval(scanPoints[i])
                    biasFits[fitType].append((truth[plane][i] - fitRate)/truth[plane][i])
                    fitRates[fitType].append(fitRate)

            # Make graph from generated scan points (truth)
            g_truth     = r.TGraph(self._nScanPoints, array('f', scanPoints), array('f', truth[plane]))
            sigmaTrue   = g_truth.GetRMS()*0.5
            peakTrue    = g_truth.GetHistogram().GetMaximum()
            meanTrue    = g_truth.GetHistogram().GetMean()
            set_graph_style(g_truth, 'VDM Simulation;#Delta {0};#mu'.format(plane), r.kBlack, 1, 20, 0.8)

            # Get width of truth
            r.gStyle.SetOptFit(0)

            if self._beamType == 'SG':
                f_truth = r.TF1('f_truth', 'gaus', 0., 1.)
                f_truth.SetParameters(1., 0.5, 0.05)
            elif self._beamType == 'DG':
                f_truth = r.TF1('f_truth','[2]*([3]*exp(-(x-[4])**2/(2*([0]*[1]/([3]*[1]+1-[3]))**2)) + (1-[3])*exp(-(x-[4])**2/(2*([0]/([3]*[1]+1-[3]))**2)))')

                f_truth.SetParNames('#Sigma','#sigma_{1}/#sigma_{2}','Amplitude','Fraction','#mu')

                f_truth.SetParameter(0, 0.5)
                f_truth.SetParameter(1, 0.5)
                f_truth.SetParameter(2, 0.001)
                f_truth.SetParameter(3, 0.7)

                f_truth.SetParLimits(0,0.5*sigmaTrue,2*sigmaTrue)
                f_truth.SetParLimits(1,0.1,10)
                f_truth.SetParLimits(2,0.95*peakTrue,1.05*peakTrue)
                f_truth.SetParLimits(3,0.,0.5)

                f_truth.FixParameter(4, meanTrue)

            #g_truth.Fit('f_truth')

            sigmaTruth = sigmaTrue*1000 #f_truth.GetParameter('#Sigma')*1000. 

            # Prepare graphs
            g_rates     = r.TGraphErrors(self._nScanPoints, array('f', scanPoints), array('f', simRates[0][plane]), \
                                         array('f', [0.01 for i in range(self._nScanPoints)]), array('f', simRates[1][plane]))
            g_biasRates = r.TGraph(self._nScanPoints, array('f', scanPoints), array('f', biasRates))

            set_graph_style(g_rates, 'VDM Simulation;#Delta {0};#mu'.format(plane), r.kBlue, 1, 21, 0.8)
            set_graph_style(g_biasRates, ';#Delta ' + plane + ';(#mu_{truth} - #mu_{#color[2]{fit}||#color[4]{rates}})/#mu_{truth}', r.kBlue, 1, 20, 0.8)

            g_fit       = {}
            g_biasFits  = {}
            for i,fitType in enumerate(fitTypes):
                g_fit[fitType]       = r.TGraph(self._nScanPoints, array('f', scanPoints), array('f', fitRates[fitType]))
                g_biasFits[fitType]  = r.TGraph(self._nScanPoints, array('f', scanPoints), array('f', biasFits[fitType]))
                set_graph_style(g_fit[fitType], 'VDM Simulation;#Delta X;#mu', graphStyles[fitType][0], 1, graphStyles[fitType][1], 0.8)
                set_graph_style(g_biasFits[fitType], ';#Delta X;', graphStyles[fitType][0], 1, graphStyles[fitType][1], 0.8)


            # Build legend
            legend = r.TLegend(0.64,0.60,0.95,0.94)
            legend.SetFillColor(0)
            legend.SetTextSize(0.045)

            legend.AddEntry(g_truth, 'MC truth ({0})'.format(self._beamType))
            #legend.AddEntry(f_truth, 'MC truth ({0})'.format(self._beamType))
            legend.AddEntry(g_rates, 'Simulated rates')
            for fitType in fitTypes:
                legend.AddEntry(g_fit[fitType], '{0} fit'.format(fitType))

            # text box with fit bias
            textBox = r.TPaveText(0.15, 0.61, 0.35, 0.89, 'NDC')
            textBox.SetFillColor(0)
            textBox.SetFillStyle(0)
            textBox.SetLineWidth(0)
            textBox.SetLineColor(0)
            textBox.SetTextSize(0.05)
            for fitType in fitTypes:
                textBox.AddText('#color[{0}]{{#Delta#Sigma/#Sigma_{{truth}} = {1:.3f}}}'.format(g_fit[fitType].GetLineColor(), abs(sigmaTruth-f_fit[plane][fitType][0])/sigmaTruth))

            
            # Set up canvas for displaying inputs and bias metrics
            pad1 = r.TPad('pad1', '', 0., 0.35, 1., 1., 0)
            pad2 = r.TPad('pad2', '', 0., 0., 1., 0.35, 0)

            pad1.SetBottomMargin(0.)
            pad1.Draw()

            pad2.SetTopMargin(0.0)
            pad2.SetBottomMargin(0.2)
            pad2.Draw()
            pad2.SetGridx()
            pad2.SetGridy()

            pad1.cd()

            g_truth.SetMinimum(-0.00005)
            g_truth.GetYaxis().SetTitleSize(0.06)
            g_truth.GetYaxis().SetTitleOffset(0.448)

            g_truth.Draw('ACP')
            #f_truth.Draw('SAME')
            g_rates.Draw('CP SAME')
            for fitType in fitTypes:
                g_fit[fitType].Draw('CP SAME')

            legend.Draw()
            textBox.Draw('SAME')


            pad2.cd()
            g_biasRates.GetYaxis().SetNdivisions(5);
            g_biasRates.GetYaxis().SetLabelSize(0.07);
            g_biasRates.GetYaxis().SetTitleSize(0.09);
            g_biasRates.GetYaxis().SetTitleOffset(0.44);
            g_biasRates.GetYaxis().CenterTitle();

            g_biasRates.GetXaxis().SetLabelSize(0.08);
            g_biasRates.GetXaxis().SetTitleSize(0.09);
            g_biasRates.GetXaxis().SetTitleOffset(0.90);

            g_biasRates.Draw('ACP')
            for fitType in fitTypes:
                g_biasFits[fitType].Draw('CP SAME')

            self._canvas.Print('{0}/bias_{1}_{2}.pdf'.format(self._plotPath, plane, self._suffix))

            pad1.Delete()
            pad2.Delete()

    def draw_beamspot_plots(self, scanPoints, beamspot, beamWidth, plane):
        '''
        For drawing beamspot plots
        '''

        pad1 = r.TPad('pad2','pad2',0,.51,1,1)
        pad1.SetBottomMargin(0.05)
        pad1.Draw()

        pad2 = r.TPad('pad1','pad1',0,0,1,.49)
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()

        pad1.SetGridy()
        pad2.SetGridy()

        g_bsX = l.make_graph([array('f', scanPoints[0]), array('f', scanPoints[1])], [array('f', beamspot[0]['X']), array('f', beamspot[1]['X'])],\
                                    'VdM scan ' + plane + ' beamspot sim;; #color[4]{X}/#color[2]{Y} (a.u.)', 22, r.kBlue)
        g_bsY = l.make_graph([array('f', scanPoints[0]), array('f', scanPoints[1])], [array('f', beamspot[0]['Y']), array('f', beamspot[1]['Y'])],\
                                    'VdM scan ' + plane + ' beamspot sim;; Y (a.u.)', 21, r.kRed)

        g_beamWidthX = l.make_graph(array('f', scanPoints[0]), array('f', beamWidth['X']), ';#Delta' + plane + ' (a.u.); RMS_{#color[4]{BSX}/#color[2]{BSY}} (a.u.)', 22, r.kBlue, doErrors = False)
        g_beamWidthY = l.make_graph(array('f', scanPoints[0]), array('f', beamWidth['Y']), ';#Delta' + plane + ' (a.u.); RMS_{BS} (a.u.)', 21, r.kRed, doErrors = False)

        pad1.cd()
        g_bsX.SetMaximum(1.1*max(beamspot[0]['X']+beamspot[0]['Y']))
        g_bsX.SetMinimum(0.8*min(beamspot[0]['X']+beamspot[0]['Y']))
        g_bsX.GetYaxis().CenterTitle()

        g_bsX.Draw('AP')
        g_bsY.Draw('P SAME')

        pad2.cd()
        g_beamWidthX.SetMaximum(1.1*max(beamWidth['X']+beamWidth['Y']))
        g_beamWidthX.SetMinimum(0.8*min(beamWidth['X']+beamWidth['Y']))
        g_beamWidthX.GetYaxis().CenterTitle()

        g_beamWidthX.Draw('ACP')
        g_beamWidthY.Draw('CP SAME')
        
        self._canvas.Print(self._plotPath + '/beamspot_{0}.pdf'.format(plane))

        pad1.Delete()
        pad2.Delete()
