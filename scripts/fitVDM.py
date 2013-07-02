import ROOT as r
import array
import math

def Residuals(g,f):

        res=r.TGraphErrors()
       
        n=g.GetN()
        x=g.GetX()
        y=g.GetY()
        xe=g.GetEX()
        ye=g.GetEY()

        for i in range(n):
                val = f.Eval(x[i])
                dev = (y[i]-val)/ye[i]

                res.SetPoint(i,x[i],dev)
                res.SetPointError(i, 0., 0.)

        return res

def Plot(graph,f,fitType,fill,scanNumber,plotName):
       
        #r.gROOT.SetBatch()
        #r.gROOT.SetStyle('Plain')
        r.gStyle.SetPalette(1)
        r.gStyle.SetOptFit(111)
        r.gStyle.SetOptStat(0)
        r.gStyle.SetTitleBorderSize(0)

        c   = r.TCanvas('c','c',600,700)

        p1  = r.TPad('p1','p1',0,0,1,.25)
        p1.SetTopMargin(0)
        p1.SetBottomMargin(0.2)
        p1.SetGridx()
        p1.SetGridy()
        p1.Draw()

        p2  = r.TPad('p2','p2',0,.25,1,1)
        p2.SetBottomMargin(0)
        p2.SetLogy()
        p2.SetGridy()
        p2.SetGridx()
        p2.Draw()
       
        peak    = f.GetParameter(2)
        sigma   = f.GetParameter(0)
        sigmaE  = f.GetParError(0)

        title = graph.GetTitle()
        xAxis = graph.GetXaxis().GetTitle()
        graph.SetTitle(title + '-plane Scan')
        graph.SetMarkerStyle(8)
        graph.SetLineColor(1)

        graph.GetYaxis().SetTitle('#mu / (I1 x I2) [a.u.]')
        #graph.GetYaxis().SetRangeUser(5e-5*peak,1.15*peak)
        graph.GetYaxis().SetRangeUser(5e-5*peak,10*peak)
        graph.GetXaxis().SetTitle('')
        graph.GetYaxis().SetLabelFont(63)
        graph.GetYaxis().SetLabelSize(16)
        graph.GetYaxis().SetTitleFont(63)
        graph.GetYaxis().SetTitleSize(18)
        graph.GetYaxis().SetTitleOffset(1.7)
        graph.GetXaxis().SetNdivisions(0)

        res = Residuals(graph,f)
        res.GetXaxis().SetLabelFont(63)
        res.GetXaxis().SetLabelSize(16)
        res.GetXaxis().SetTitleFont(63)
        res.GetXaxis().SetTitleSize(18)
        res.GetXaxis().SetTitleOffset(3)
        res.GetXaxis().SetTitle(xAxis)
        res.GetYaxis().SetLabelFont(63)
        res.GetYaxis().SetLabelSize(16)
        res.GetYaxis().SetTitleFont(63)
        res.GetYaxis().SetTitleSize(18)
        res.GetYaxis().SetTitleOffset(1.5)
        res.GetYaxis().SetNdivisions(5);
        res.GetYaxis().SetTitle('Residuals [#sigma]')

        p2.cd()
        graph.Draw('AP')

        if fitType is 'doubleGaussian':

                sigma2 = f.GetParameter(1)
                h = f.GetParameter(3)
                s2 = sigma/(h*sigma2+1-h)
                a1 = peak*h
                a2 = peak*(1-h)
                s1 = sigma*sigma2/(h*sigma2+1-h)

                f1 = r.TF1('f1','gaus',-5.,5.)
                f1.SetParameters(a1,f.GetParameter(4),s1)
                f1.SetLineColor(2)
                f2 = r.TF1('f2','gaus',-5.,5.)
                f2.SetParameters(a2,f.GetParameter(4),s2)
                f2.SetLineColor(3)
                f3 = r.TF1('f3','pol0',-5.,5)
                f3.SetParameter(0,f.GetParameter(5))
                f3.SetLineColor(4)

                f1.Draw('same')
                f2.Draw('same')
                f3.Draw('same')

        elif fitType is 'singleGaussian':

                f1 = r.TF1('f1','gaus',-5.,5.)
                f1.SetParameters(peak,f.GetParameter(1),sigma)
                f1.SetLineColor(2)
                f2 = r.TF1('f2','pol0',-5.,5)
                f2.SetParameter(0,f.GetParameter(3))
                f2.SetLineColor(4)

                f1.Draw('same')
                f2.Draw('same')

        pad = r.TPaveText(0.15,0.78,0.30,0.88,'NDC')
        pad.SetTextAlign(12) # left alignment
        pad.SetTextSize(0.03)
        pad.SetFillColor(0)
        pad.AddText('CMS Preliminary')
        pad.AddText('VdM Scan '+scanNumber+': Fill '+str(fill))
        pad.Draw('same')

        pad2 = r.TPaveText(0.4,0.2,0.6,0.3,'NDC')
        pad2.SetTextSize(0.04)
        pad2.SetFillColor(0)
        pad2.AddText('#Sigma = '+str(round(1e3*sigma,3))+ '#pm' + str(round(1e3*sigmaE,3))+ '[#mum]')
        pad2.Draw('same')
       
        p1.cd()
        res.SetMarkerStyle(8)
        res.Draw('AP')

        c.Print(plotName+'.pdf')


def Plot2D(graph,f,fill,f2D, fit, plotName):

    #print "fit ", f2D.GetParameter(0), f2D.GetParError(0), f2D.GetParName(0), f2D.GetNumberFreeParameters()

    #r.gROOT.SetBatch()	
    #r.gROOT.SetStyle("Plain")
    r.gStyle.SetPalette(1)
    r.gStyle.SetOptFit(0)
    r.gStyle.SetOptStat(0)
    r.gStyle.SetTitleBorderSize(0)

    c = r.TCanvas("c","c",600,700)

    p1 = r.TPad('p1','p1',0,0,1,.25)
    p1.SetTopMargin(0)
    p1.SetBottomMargin(0.2)
    p1.SetGridx()
    p1.SetGridy()
    p1.Draw()

    p2  = r.TPad('p2','p2',0,.25,1,1)
    p2.SetBottomMargin(0)
    p2.SetLogy()
    p2.SetGridy()
    p2.SetGridx()
    p2.Draw()

    peak = f.Eval(0)

    title = graph.GetTitle()
    graph.SetTitle(title) #[0].upper() + "-plane BCID " + title[1:])

    graph.GetYaxis().SetTitle('#mu / (I1 x I2) [a.u.]')
    graph.GetYaxis().SetRangeUser(5e-5*peak,100*peak)
    #graph.GetXaxis().SetTitle("#Delta [mm]")
    #graph.GetXaxis().SetTitle('')
    graph.SetMarkerStyle(8)
    graph.GetYaxis().SetLabelFont(63)
    graph.GetYaxis().SetLabelSize(16)
    graph.GetYaxis().SetTitleFont(63)
    graph.GetYaxis().SetTitleSize(18)
    graph.GetYaxis().SetTitleOffset(1.5)
    graph.GetXaxis().SetNdivisions(0)


    res = Residuals(graph,f)
    res.GetXaxis().SetLabelFont(63)
    res.GetXaxis().SetLabelSize(16)
    res.GetYaxis().SetLabelFont(63)
    res.GetYaxis().SetLabelSize(16)
    res.GetXaxis().SetTitleFont(63)
    res.GetXaxis().SetTitleSize(18)
    res.GetXaxis().SetTitleOffset(3)
    res.GetYaxis().SetTitleFont(63)
    res.GetYaxis().SetTitleSize(18)
    res.GetYaxis().SetTitleOffset(1.5)

    p2.cd()
    graph.Draw('AP')

    f.SetLineColor(r.kRed)
    f.Draw("same")

    pad = r.TPaveText(0.15,0.75,0.35,0.88,"NDC")
    pad.SetTextAlign(12) # left alignment
    pad.SetTextSize(0.04)
    pad.SetFillColor(0)
    pad.AddText("CMS Preliminary")
    pad.AddText("VdM Scan: Fill "+str(fill))
    pad.Draw("same")

    #pad2 = r.TPaveText(0.4,0.2,0.6,0.3,"NDC")
    #pad2.SetTextSize(0.04)
    #pad2.SetFillColor(0)
    #pad2.AddText("#Sigma = "+str(round(1e3*sigma,3))+ "#pm" + str(round(1e3*sigmaE,3))+ "[#mum]")
    #pad2.Draw("same")

    pad3 = r.TPaveText(0.6,0.65,0.99,0.99,"NDC")
    #pad3.SetLineStyle(0)
    #pad3.SetLineColor(r.kBlack)
    #pad3.SetLineWidth(2)
    pad3.SetTextAlign(12) # left alignment
    pad3.SetTextSize(0.025)
    pad3.SetFillColor(0)
    pad3.AddText("#chi^{2}/ndof " + str(round(fit.Chi2(),1))+"/"+str(fit.Ndf()))

    for i in range(0,f2D.GetNumberFreeParameters()):
        pad3.AddText(str(f2D.GetParName(i)) + "  " + str(round(f2D.GetParameter(i),5)) + " #pm " +  str(round(f2D.GetParError(i),5)))

    pad3.Draw("same")

    p1.cd()
    res.SetMarkerStyle(8)
    res.GetXaxis().SetTitle("#Delta [mm]")
    res.GetYaxis().SetTitle("Residuals [#sigma]")
    res.Draw("AP")

    c.Print(plotName + '.pdf')


def Fit(graph, fitType, doPlot, fill, scanNumber, plotName = 'test'):

        # pre guess the fit parameters
        ExpSigma    = graph.GetRMS()*0.5
        ExpPeak     = graph.GetHistogram().GetMaximum()
        ExpMean     = graph.GetHistogram().GetMean()

        # double gaussian formula with substition to effective width and widths ration
        # Sigma_eff = h*sigma1 + (1-h)*sigma2
        # Sigma2 = sigma1/sigma2
        # [0] -> [0]*[1]/([3]*[1]+1-[3])
        # [1] -> [0]/([3]*[1]+1-[3])

        fits = {}

        # Double Gaussian
        fits['doubleGaussian'] = r.TF1('f_doubleGaussian','[2]*([3]*exp(-(x-[4])**2/(2*([0]*[1]/([3]*[1]+1-[3]))**2)) + (1-[3])*exp(-(x-[4])**2/(2*([0]/([3]*[1]+1-[3]))**2))) + [5]')
        fits['doubleGaussian'].SetParNames('#Sigma','#sigma_{1}/#sigma_{2}','Amplitude','Fraction','Mean','Constant')
        fits['doubleGaussian'].SetParameters(ExpSigma,1.,ExpPeak,0.99,0.)
        fits['doubleGaussian'].SetParLimits(0,0.5*ExpSigma,2*ExpSigma)
        fits['doubleGaussian'].SetParLimits(1,0.1,10)
        fits['doubleGaussian'].SetParLimits(2,0.2*ExpPeak,5*ExpPeak)
        fits['doubleGaussian'].SetParLimits(3,0.,1.)
        fits['doubleGaussian'].SetParLimits(5,0., 1e-2*ExpPeak) 
        fits['doubleGaussian'].SetLineColor(1)

        # Skewed Gaussian
        fits['skewGaussian'] = r.TF1('f_skewGaussian','[4] + [2]*exp(-(x-[1])**2/(2*[0]**2*(1 + [3]*TMath::Sign(1., x-[1]))))')
        fits['skewGaussian'].SetParNames('#Sigma','Mean','Amplitude','#alpha', 'Constant')
        fits['skewGaussian'].SetParameters(ExpSigma, 0., ExpPeak, 0.5)
        fits['skewGaussian'].SetParLimits(0,0.5*ExpSigma,2*ExpSigma)
        fits['skewGaussian'].SetParLimits(2,0.2*ExpPeak,5.*ExpPeak)
        fits['skewGaussian'].SetParLimits(3,0.,0.99)
        fits['skewGaussian'].SetParLimits(4,0.,1e-2*ExpPeak)
        fits['skewGaussian'].SetLineColor(1)

        # Single Gaussian
        fits['singleGaussian'] = r.TF1('f_singleGaussian','[3] + [2]*exp(-(x-[1])**2/(2*[0]**2))')
        fits['singleGaussian'].SetParNames('#Sigma','Mean','Amplitude','Constant')
        fits['singleGaussian'].SetParameters(ExpSigma,0.,ExpPeak)
        fits['singleGaussian'].SetParLimits(0,0.5*ExpSigma,2*ExpSigma)
        fits['singleGaussian'].SetParLimits(2,0.2*ExpPeak,5*ExpPeak)
        fits['singleGaussian'].SetParLimits(3,0.,1e-2*ExpPeak)
        fits['singleGaussian'].SetLineColor(1)

        for j in range(5):
                fit = graph.Fit('f_'+fitType,'SQ')
                if fit.CovMatrixStatus() == 3 and fit.Chi2()/fit.Ndf() < 2: break

        #if fitType in ['doubleGaussian', 'skewGaussian']:
        #    if fit.CovMatrixStatus() != 3 or fits[fitType].GetParameter(3) > 0.975 or fits[fitType].GetParameter(3) < 0.025:  
        #        for j in range(5):
        #                fit = graph.Fit('f_singleGaussian','SQ')
        #                if fit.CovMatrixStatus() == 3 and fit.Chi2()/fit.Ndf() < 2: break
        #        fitType = 'singleGaussian'

        #print fit.Status(), fit.CovMatrixStatus()      

        peak    = fits[fitType].GetParameter(2)
        peakE   = fits[fitType].GetParError(2)
        sigma   = fits[fitType].GetParameter(0)*1000
        sigmaE  = fits[fitType].GetParError(0)*1000
        cov     = fit.CovMatrix(0,2)
        const   = fits[fitType].GetParameter('Constant')
        corr    = peak/fits[fitType](0)
       
        xmax    = r.TMath.MaxElement(graph.GetN(),graph.GetX())
       
        sigma   = (const*2*xmax/math.sqrt(2*math.pi) + sigma*peak)/(const+peak)
        peak    = const + peak

        if doPlot:
                Plot(graph,fits[fitType],fitType,fill,scanNumber,plotName)    

        return [sigma, sigmaE, peak, peakE, cov, corr, fits[fitType]]


def Fit2D(graph2D, graph1DX, graph1DY, doPlot, fill, Results1DX, Results1DY, plotName = 'test'):
    '''
    R(x,y) = A [  f  exp(-{x-x01}**2/{2*sigx**2})    * exp(-{y-y01}**2/{2*sigy**2})
             + (1-f) exp(-{x-x02}**2/{2*(Sx*sigx)**2}) * exp(-{y-y02}**2/{2*(Sy*sigy)**2})

    sigx    - [0],
    Sx      - [1]
    sigy    - [2]
    Sy      - [3]
    x01     - [4], x02 - [5], y01 - [6], y02 - [7]
    f       - [8]
    A       - [9]
    '''

    #print "in Fit2D", graph2D

    FitResults = {'X':Results1DX[6], 'Y':Results1DY[6]}

    sigmaEffX   = FitResults['X'].GetParameter('#Sigma')
    sigmaRatioX = FitResults['X'].GetParameter('#sigma_{1}/#sigma_{2}')
    fractionX   = FitResults['X'].GetParameter('Fraction')
    meanX       = FitResults['X'].GetParameter('Mean')

    sigmaEffY   = FitResults['Y'].GetParameter('#Sigma')
    sigmaRatioY = FitResults['Y'].GetParameter('#sigma_{1}/#sigma_{2}')
    fractionY   = FitResults['Y'].GetParameter('Fraction')
    meanY       = FitResults['Y'].GetParameter('Mean')

    ExpSigmaX   = graph1DX.GetRMS()*0.5
    ExpPeakX    = graph1DX.GetHistogram().GetMaximum()

    ExpSigmaY   = graph1DY.GetRMS()*0.5
    ExpPeakY    = graph1DY.GetHistogram().GetMaximum()

    ExpPeak     = (ExpPeakX + ExpPeakY)/2.

    # Define functions for fitting vdm profiles
    f2D = r.TF2('f2D','[9]*([8]*exp(-(x - [4])**2/(2*([0]*[1]/(1 + [8]*([1] - 1)))**2) \
                                        - (y - [6])**2/(2*([2]*[3]/(1 + [8]*([3] - 1)))**2)) \
                                        + (1 - [8])*exp(-(x-[5])**2/(2*([0]/(1 + [8]*([1] - 1)))**2) \
                                        - (y - [7])**2/(2*([2]/(1 + [8]*([3] - 1)))**2)))')

    f2D.SetParNames('#Sigma_{x}', '#sigma_{1,x}/#sigma_{2,x}', '#Sigma_{y}', '#sigma_{1,y}/#sigma_{2,y}',\
                         'x_{1}', 'x_{2}', 'y_{1}', 'y_{2}', 'Fraction', 'Amp')
    f2D.SetParameters(sigmaEffX, sigmaRatioX, sigmaEffY, sigmaRatioY, meanX, meanX, meanY, meanY, (fractionX + fractionY)/2., ExpPeak/2.)

    f2D.SetParLimits(0, 0.5*sigmaEffX, 2*sigmaEffX)
    f2D.SetParLimits(1, 0.01, 1.99)
    f2D.SetParLimits(2, 0.5*sigmaEffY, 2*sigmaEffY)
    f2D.SetParLimits(3, 0.01, 1.99)
    f2D.SetParLimits(4, -1.1*meanX, 1.1*meanX)
    f2D.SetParLimits(5, -1.1*meanX, 1.1*meanX)
    f2D.SetParLimits(6, -1.1*meanY, 1.1*meanY)
    f2D.SetParLimits(7, -1.1*meanY, 1.1*meanY)
    f2D.SetParLimits(8, 0, 1.)
    f2D.SetParLimits(9, 0.9*ExpPeak, 1.1*ExpPeak)
    
    for l in range(5):
        fit2D = graph2D.Fit("f2D","SQ")
        if fit2D.CovMatrixStatus() == 3 and fit2D.Chi2()/fit2D.Ndf() < 2: break
        
        #print "Covariance matrix", fit2D.GetCovarianceMatrix().Print()

    f2Dfunc = graph2D.FindObject("f2D")

    params = [f2D.GetParameter(i) for i in range(0,10)]

    xmax = r.TMath.MaxElement(graph1DX.GetN(),graph1DX.GetX())
    ymax = r.TMath.MaxElement(graph1DY.GetN(),graph1DY.GetX())

    f_ProjX = r.TF1('f_ProjX','[9]*([8]*exp(-(x - [4])**2/(2*([0]*[1]/(1 + [8]*([1] - 1)))**2) \
                                        - ([6])**2/(2*([2]*[3]/(1 + [8]*([3] - 1)))**2)) \
                                        + (1 - [8])*exp(-(x-[5])**2/(2*([0]/(1 + [8]*([1] - 1)))**2) \
                                        - ([7])**2/(2*([2]/(1 + [8]*([3] - 1)))**2)))', -0.5, 0.5)

    f_ProjY = r.TF1('f_ProjY','[9]*([8]*exp(-([4])**2/(2*([0]*[1]/(1 + [8]*([1] - 1)))**2) \
                                        - (x - [6])**2/(2*([2]*[3]/(1 + [8]*([3] - 1)))**2)) \
                                        + (1 - [8])*exp(-([5])**2/(2*([0]/(1 + [8]*([1] - 1)))**2) \
                                        - (x - [7])**2/(2*([2]/(1 + [8]*([3] - 1)))**2)))', -0.5, 0.5)

    #f_ProjY = r.TF1('f_ProjY','f_ProjGauss', -1.5 * ymax, 1.5* ymax)
    #f_ProjX = r.TF1('f_ProjX','f_ProjGauss', -1.5 * xmax, 1.5* xmax)

    for i in range(0,10):
        f_ProjX.SetParameter(i, params[i])
        f_ProjY.SetParameter(i, params[i])

    if doPlot:
        Plot2D(graph1DX, f_ProjX, fill, f2Dfunc, fit2D, plotName + '_X')	
        Plot2D(graph1DY, f_ProjY, fill, f2Dfunc, fit2D, plotName + '_Y')	

    fitChi2 = fit2D.Chi2()/fit2D.Ndf()

    xmin = array.array('d',[0.0])
    ymin = array.array('d',[0.0])

    f2Dhelper = r.TF2("f2Dhelper","(-1.)*f2D")
    f2Dhelper.GetMinimumXY(xmin,ymin)     

    return [params, xmin, ymin, f2D, f_ProjX, f_ProjY]
