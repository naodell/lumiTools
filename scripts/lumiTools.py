import sys, math, csv, json, pprint, re, os, commands, datetime
from array import array
import ROOT as r

### Defines some useful functions
def read_json_file(filename):
    inFile  = open(filename, 'rb')
    data    = json.loads(inFile)
    pprint(data)
    inFile.close()
    return data


def read_csv_file(filename, fieldNames=None):
    inFile      = open(filename, 'rb')
    dataIter    = csv.DictReader(inFile, fieldnames=fieldNames)
    data        = list(dataIter)
    inFile.close()
    return data


def set_style(graph, color = 1, width = 2):
    graph.SetLineColor(color)
    graph.SetLineWidth(width)


def average(inputList):
    return sum(inputList)/len(inputList)


def dict_map(keyList, valueList):
    return dict(zip(keyList, valueList))


def format_histogram_time(hist):
    hist.GetXaxis().SetTimeDisplay(1);
    hist.GetXaxis().SetLabelOffset(0.03);
    hist.GetXaxis().SetTimeFormat("#splitline{%H:%M:%S}{%d/%m/%y}"); 
    hist.GetXaxis().SetTimeOffset(3600,"gmt");

def ratio_error(numer, denom):
    '''
    calculates the error of a ratio
    '''
    err = math.sqrt(pow(numer[1]/denom[0], 2) + pow(numer[0]*denom[1]/pow(denom[0], 2), 2))
    return err
    
def make_graph(varX, varY, title, marker, color, yScale = (0.9, 1.1), doErrors = True):
    '''
    A wrapper for TGraphErrors
    '''

    if doErrors:
        graph = r.TGraphErrors(len(varX[0]), varX[0], varY[0], varX[1], varY[1])
    else:
        graph = r.TGraph(len(varX), varX, varY)

    graph.SetTitle(title)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerColor(color)
    graph.SetFillColor(0)
    graph.SetLineColor(color)
    graph.GetYaxis().SetTitleOffset(1.5)
    #graph.GetYaxis().SetRangeUser(min(varY[0])*yScale[0], max(varY[0])*yScale[1])

    return graph

def ratio_graph(g_numer, g_denom, title, color, marker, rangeY):
    '''
    Makes a ratio of two graphs with errors
    '''
    
    bins = (g_numer.GetN(), g_denom.GetN())
    numX = (g_numer.GetX(), g_numer.GetEX())
    numY = (g_numer.GetY(), g_numer.GetEY())
    denX = (g_denom.GetX(), g_denom.GetEX())
    denY = (g_denom.GetY(), g_denom.GetEY())
        
    ratio   = []
    rErr    = []
    xVal    = []
    xErr    = []
    numBins = [numX[0][j] for j in range(bins[0])]
    for i in range(bins[1]):
        if denY[0][i] < 1e-9: continue
        index = numBins.index(float(denX[0][i]))

        ratio.append(100*(denY[0][i] - numY[0][index])/denY[0][i])
        rErr.append(100*ratio_error((numY[0][index], numY[1][index]), (denY[0][i], denY[1][i])))
        xVal.append(denX[0][i])
        xErr.append(1e-3)

    g_ratio = make_graph((array('f', xVal), array('f', xErr)), (array('f', ratio), array('f', rErr)), title, marker, color)

    return g_ratio
        



def print_vdm_tables(algo, bcids, fitResults, resultFile, doTwiki=False, doLatex=True):
    '''Prints tables with fit results, if doTwiki is false then printed for Latex'''

    if doTwiki:
        header = '   | bunch crossing | &Sigma;<sub>X</sub> (&mu;m) | &Sigma;<sub>Y</sub> (&mu;m) | &mu;<sub>Y,max</sub>/I<sub>1</sub>xI<sub>x</sub> (10 <sup>-22</sup>  p) | &mu;<sub>X,max</sub>/I<sub>1</sub>xI<sub>x</sub> (10 <sup>-22</sup> p) | &sigma;<sub>vis</sub> (10<sup>-30</sup> cm<sup>-2</sup>) |\n'
        outString = '  | {0} | {1[0][0]:.2f} &pm; {1[0][1]:.2f} | {1[1][0]:.2f} &pm; {1[1][1]:.2f} | {1[0][2]:.5f} &pm; {1[0][3]:.5f} | {1[1][2]:.5f} &pm; {1[1][3]:.5f} | {1[2][0]:.2f}  &pm; {1[2][1]:.2f}| \n'
    elif doLatex:
        header = '\tbunch crossing & $\Sigma_{X}$ ($\mu$m) & $\Sigma_{Y}$ ($\mu$m) & $\mu_{X,max}/I_{1}\\timesI_{2} (10^{22} p)$ & $\mu_{Y,max}/I_{1}\\timesI_{2} (10^{22} p)$ & $\sigma_{vis}$ (10^{-30}$ cm$^{-2}) \\\ \hline \n'
        outString = '\t{0} & {1[0][0]:.2f} $\pm$ {1[0][1]:.2f} & {1[1][0]:.2f} $\pm$ {1[1][1]:.2f} & {1[0][2]:.5f} $\pm$ {1[0][3]:.5f} & {1[1][2]:.5f} $\pm$ {1[1][3]:.5f} & {1[2][0]:.2f} $\pm$ {1[2][1]:.2f} \\\ \hline \n'
    else:
        header = 'bunch crossing,sigmaX,sigErrX,sigmaY,sigErrY,muX,muErrX,muY,muErrY,sigmaVis,sigErrVis\n'
        outString = '{0},{1[0][0]:.2f},{1[0][1]:.2f},{1[1][0]:.2f},{1[1][1]:.2f},{1[0][2]:.5f},{1[0][3]:.5f},{1[1][2]:.5f},{1[1][3]:.5f},{1[2][0]:.2f},{1[2][1]:.2f}\n'


    if doLatex or doTwiki:
        resultFile.write('\n\t{}\n'.format(algo))
        resultFile.write(header)
    else:
        resultFile.write(header)

    avgResult = [[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0.]]
    for bcid in bcids[:35]:
        bunch = str(bcid)
        resultFile.write(outString.format(bunch, fitResults[algo][bunch]))

        avgResult[0][0] += fitResults[algo][bunch][0][0]/pow(fitResults[algo][bunch][0][1], 2)
        avgResult[0][1] += pow(1/fitResults[algo][bunch][0][1], 2)
        avgResult[0][2] += fitResults[algo][bunch][0][2]/pow(fitResults[algo][bunch][0][3], 2)
        avgResult[0][3] += pow(1/fitResults[algo][bunch][0][3], 2)           
        avgResult[1][0] += fitResults[algo][bunch][1][0]/pow(fitResults[algo][bunch][1][1], 2)
        avgResult[1][1] += pow(1/fitResults[algo][bunch][1][1], 2)           
        avgResult[1][2] += fitResults[algo][bunch][1][2]/pow(fitResults[algo][bunch][1][3], 2)
        avgResult[1][3] += pow(1/fitResults[algo][bunch][1][3], 2)           
        avgResult[2][0] += fitResults[algo][bunch][2][0]/pow(fitResults[algo][bunch][2][1], 2)
        avgResult[2][1] += pow(1/fitResults[algo][bunch][2][1], 2)

    avgResult[0][0] = avgResult[0][0]/avgResult[0][1]
    avgResult[0][1] = math.sqrt(1/avgResult[0][1])
    avgResult[0][2] = avgResult[0][2]/avgResult[0][3]
    avgResult[0][3] = math.sqrt(1/avgResult[0][3])
    avgResult[1][0] = avgResult[1][0]/avgResult[1][1]
    avgResult[1][1] = math.sqrt(1/avgResult[1][1])
    avgResult[1][2] = avgResult[1][2]/avgResult[1][3]
    avgResult[1][3] = math.sqrt(1/avgResult[1][3])
    avgResult[2][0] = avgResult[2][0]/avgResult[2][1]
    avgResult[2][1] = math.sqrt(1/avgResult[2][1])

    resultFile.write(outString.format('Average', avgResult))
    #print outString.format('Average', avgResult)

    resultFile.write('\n')

def draw_summary_plots(bcids, algo, scanResults, fill, savePath):
    '''
    Draws summary plots for bunch-by-bunch VDM results
    '''
    g_SigX      = []
    g_MuX       = []
    g_SigY      = []
    g_MuY       = []
    g_SigVis    = []

    for scan, scanResult in enumerate(scanResults):
        result = []
        sigma  = {'X':[], 'errX':[], 'Y':[], 'errY':[]}
        mu     = {'X':[], 'errX':[], 'Y':[], 'errY':[]}
        sigVis = {'val':[], 'err':[]}

        for bcid in bcids:
            bunch = str(bcid)

            sigma['X'].append(scanResult[algo][bunch][0][0])
            sigma['errX'].append(scanResult[algo][bunch][0][1])
            mu['X'].append(scanResult[algo][bunch][0][2])
            mu['errX'].append(scanResult[algo][bunch][0][3])

            sigma['Y'].append(scanResult[algo][bunch][1][0])
            sigma['errY'].append(scanResult[algo][bunch][1][1])
            mu['Y'].append(scanResult[algo][bunch][1][2])
            mu['errY'].append(scanResult[algo][bunch][1][3])

            sigVis['val'].append(scanResult[algo][bunch][2][0])
            sigVis['err'].append(scanResult[algo][bunch][2][1])

        varX = (array('f', bcids), array('f', [0.01 for i in range(len(bcids))]))

        g_SigX.append(make_graph(varX, (array('f', sigma['X']), array('f', sigma['errX'])),
                                ';Bunch Crossing;#Sigma_{X} (#mum)', 22-scan, 1+scan, (0.88, 1.10)))
        g_SigY.append(make_graph(varX, (array('f', sigma['Y']), array('f', sigma['errY'])), 
                                ';Bunch Crossing;#Sigma_{Y} (#mum)', 22-scan, 1+scan, (0.9, 1.15)))
        g_MuX.append(make_graph(varX, (array('f', mu['X']), array('f', mu['errX'])), 
                                ';Bunch Crossing;#mu_{X} (10^{-22})', 22-scan, 1+scan, (0.8, 1.10)))
        g_MuY.append(make_graph(varX, (array('f', mu['Y']), array('f', mu['errY'])), 
                                ';Bunch Crossing;#mu_{Y} (10^{-22})', 22-scan, 1+scan, (0.8, 1.10)))
        g_SigVis.append(make_graph(varX, (array('f', sigVis['val']), array('f', sigVis['err'])), 
                                ';Bunch Crossing;#sigma_{vis} [10^{-30} cm^{2}]', 22-scan, 1+scan, (0.8, 1.15))) 

    canvas = r.TCanvas('canvas', 'canvas', 800, 600)
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)

    #f_avg = r.TF1('f_avg', 'pol0', 0, 3564)
    #f_avg.SetParName(0, '<#sigma_{vis}>')
    canvas.Print(savePath+'/sumGraphs_'+fill+'.pdf[')

    for i, graphs in enumerate([g_SigX, g_MuX, g_SigY, g_MuY, g_SigVis]):
        graphs[0].Draw('AP')
        #graphs[0].Fit('f_avg', 'SQ')
        #f_avg.SetLineColor(r.kRed)
        #f_avg.Draw('same')

        for g in graphs[1:]:
            g.Draw('P same')
            g.Draw('P same')

        legend = r.TLegend(0.83,0.81,0.98,0.95)
        legend.SetFillColor(0)
        legend.SetTextSize(0.045)
        legend.AddEntry(graphs[0], 'scan 1')
        for j,g in enumerate(graphs[1:]):
            legend.AddEntry(g, 'scan {}'.format(j+2))
        legend.Draw()

        text = r.TPaveText()
        text = r.TPaveText(0.15,0.78,0.40,0.88,'NDC')
        text.SetTextAlign(12) # left alignment
        text.SetTextSize(0.03)
        text.SetBorderSize(0)
        text.SetFillStyle(0)
        text.AddText('CMS Preliminary 2012')
        #text.AddText('#sqrt s = 8 TeV')
        text.AddText('November VdM Scan (Fill '+fill+')')
        text.Draw('same')

        canvas.SaveAs(savePath+'/sumGraphs_'+fill+'_'+str(i+1)+'.pdf')
        #canvas.Print(savePath+'/sumGraphs_'+fill+'.pdf')

    canvas.Print(savePath+'/sumGraphs_'+fill+'.pdf]')

    return (g_SigX, g_SigY)


def get_fill_pattern(schemeName):
    '''
    Adopted from Jeroen's calc_pixel_afterglow.py script,
    <svn url>
    '''

    lpc_url = 'http://lpc.web.cern.ch/lpc/documents/FillPatterns/'
    scheme_file_name = '%s.txt' % schemeName
    scheme_file_url = '%s/%s' % (lpc_url, scheme_file_name)

    ##########
    # Get the filling scheme text file.
    # First see if it is available locally, otherwise get it from the
    # LPC web site.
    ##########

    if not os.path.isfile(scheme_file_name):
        print 'Trying to get filling scheme file from the LPC web site'
        cmd = 'curl %s -o %s' % (scheme_file_url, scheme_file_name)
        (status, output) = commands.getstatusoutput(cmd)
        if status != 0:
            print >> sys.stderr, \
                'ERROR Could not get filling scheme from LPC web site: %s' % \
                output
            sys.exit(1)
        else:
            print 'Success!'
    else:
        print 'Using locally stored filling scheme file'

    # Read and parse the filling scheme definition file.
    try:
        scheme_file = open(scheme_file_name, 'r')
        scheme_def_lines = scheme_file.readlines()
        scheme_file.close()
    except IOError:
        print >> sys.stderr, \
            'ERROR Could not read filling scheme definition file'
        sys.exit(1)

    regexp = re.compile('^# *BEAM([1,2]{1}) *BEAM([1,2]{1}) slots *$')
    beam_patterns = {1 : [], 2 : []}
    data_block = None
    processing_data_block = False

    for line in scheme_def_lines:

        if line.startswith('#'):
            if processing_data_block:
                data_block = None
                processing_data_block = False
            match = regexp.match(line)
            if match:
                # Found the start of either the beam1 or the beam2
                # data block.
                # DEBUG DEBUG DEBUG
                assert set(match.groups()) == set(['1', '2'])
                # DEBUG DEBUG DEBUG end
                data_block = int(match.group(1))
            continue

        if not data_block:
            continue

        processing_data_block = True

        bx_nr = int(line.split()[0])
        beam_patterns[data_block].append(bx_nr)

    beam1_pattern = beam_patterns[1]
    beam2_pattern = beam_patterns[2]

    ## DEBUG DEBUG DEBUG
    assert len(beam1_pattern) == len(beam2_pattern)
    ## DEBUG DEBUG DEBUG end

    collidingPattern = []
    for i in range(3564):
        collidingPattern.append(False)
        if i in beam1_pattern and i in beam2_pattern:
            collidingPattern[i] = True

    return collidingPattern

def overlay_widths(g_Sig, bcids, doHF, doPixel, savePath, scan):
    '''
    For overlaying widths from BCM1F measurement and comparing
    to HF and/or the pixels
    '''

    ### Pixel data ###
    pixelBX   = [1, 721, 1581, 2161, 2241] 
    pixelSigX = [31.90, 31.55, 31.78, 31.59, 31.93]
    pixelErrX = [0.19, 0.2, 0.19, 0.18, 0.18]
    pixelSigY = [26.32, 25.93, 26.20, 26.35, 26.80]
    pixelErrY = [0.34, 0.36, 0.35, 0.35, 0.36]

    for i in range(len(pixelBX)):
        pixelErrX[i] = 1e-2*pixelErrX[i]*pixelSigX[i] 
        pixelErrY[i] = 1e-2*pixelErrY[i]*pixelSigY[i]

    g_PixSigX = make_graph((array('f', pixelBX), array('f', [0.01 for i in range(5)])),
    (array('f', pixelSigX), array('f', pixelErrX)), 'Pixels', 21, r.kBlue)

    g_PixSigY = make_graph((array('f', pixelBX), array('f', [0.01 for i in range(5)])), 
    (array('f', pixelSigY), array('f', pixelErrY)), 'Pixels', 21, r.kBlue)


    ### HF data ###
    #hfData = read_csv_file('data/VDM/HF/forComparisonWithATLAS1_July.txt') 
    hfData = read_csv_file('data/November/VDM_Nov_2012/HF_results/forComparisonWithATLAS' + scan + '.txt') 
    hfSigX = []
    hfErrX = []
    hfSigY = []
    hfErrY = []

    for row in hfData:
        #hfSigX.append(float(row['sigmaX']))
        #hfErrX.append(float(row['sigmaXe']))
        #hfSigY.append(float(row['sigmaY']))
        #hfErrY.append(float(row['sigmaYe']))

        hfSigX.append(float(row['CSx']))
        hfErrX.append(float(row['errCSx']))
        hfSigY.append(float(row['CSy']))
        hfErrY.append(float(row['errCSy']))

    g_HFSigX = make_graph((array('f', bcids), array('f', [0.01 for i in range(len(bcids))])),
    (array('f', hfSigX), array('f', hfErrX)), ';Bunch Crossing;#Sigma_{x} [#mu m]', 23, r.kRed)

    g_HFSigY = make_graph((array('f', bcids), array('f', [0.01 for i in range(len(bcids))])), 
    (array('f', hfSigY), array('f', hfErrY)), ';Bunch Crossing;#Sigma_{y} [#mu m]', 23, r.kRed)

    g_RatioHFX  = ratio_graph(g_Sig[0], g_HFSigX, ';Bunch Crossing; Percent Difference #Sigma_{X}', r.kBlue, 23, (-1., 2.5))
    g_RatioPixX = ratio_graph(g_Sig[0], g_PixSigX, ';Bunch Crossing;', r.kRed, 21, (-2., 3.))
    g_RatioHFY  = ratio_graph(g_Sig[1], g_HFSigY, ';Bunch Crossing; Percent Difference #Sigma_{Y}', r.kBlue, 23, (-1., 2.5))
    g_RatioPixY = ratio_graph(g_Sig[1], g_PixSigY, ';Bunch Crossing;', r.kRed, 21, (-4., 4.))

    ### Set styles ###
    g_Sig[0].SetLineColor(r.kBlack)
    g_Sig[0].SetMarkerColor(r.kBlack)
    g_Sig[0].SetMarkerStyle(21)
    g_Sig[1].SetLineColor(r.kBlack)
    g_Sig[1].SetMarkerColor(r.kBlack)
    g_Sig[1].SetMarkerStyle(21)

    ### Make overlays ###

    canvas = r.TCanvas('canvas', 'canvas', 800, 600)
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    #pad1 = r.TPad('p1','p1',0,0,1,.25)
    #pad2 = r.TPad('p2','p2',0,.25,1,1)
    #pad1.SetGridx()
    #pad2.SetGridx()
    #pad2.SetGridy()

    text = r.TPaveText(0.15,0.78,0.40,0.88,'NDC')
    text.SetTextAlign(12) # left alignment
    text.SetTextSize(0.03)
    text.SetFillColor(0)
    #text.AddText('CMS Preliminary 2012')
    text.AddText('VdM Scan ' + scan + ': Fill 3316')

    legend = r.TLegend(0.69,0.75,0.97,0.94)
    legend.SetFillColor(0)
    legend.SetTextSize(0.045)
    legend.AddEntry(g_Sig[0], 'BCM1F')
    legend.AddEntry(g_HFSigX, 'HF')
    #legend.AddEntry(g_PixSigX, 'Pixel')
    #legend.AddEntry(g_RatioHFX, 'HF - BCM1F')
    #legend.AddEntry(g_RatioPixX, 'Pixel - BCM1F')

    #pad1.cd()
    g_HFSigX.Draw('AP')
    #g_PixSigX.Draw('P same')
    g_Sig[0].Draw('P same')
    #g_RatioHFX.Draw('AP')
    #g_RatioPixX.Draw('P same')
    legend.Draw()
    text.Draw('same')

    #pad2.cd()
    #g_RatioHFX = g_Sig[0].GetHistogram()#Divide(g_HFSigX) 
    ##g_RatioHFX.Draw()

    canvas.SaveAs(savePath+'/overlays_scan' + scan + '_X.pdf')

    #pad1.cd()
    g_HFSigY.Draw('AP')
    #g_PixSigY.Draw('P same')
    g_Sig[1].Draw('P same')
    #g_RatioHFY.Draw('AP')
    #g_RatioPixY.Draw('P same')
    legend.Draw()
    text.Draw('same')

    canvas.SaveAs(savePath+'/overlays_scan' + scan + '_Y.pdf')

