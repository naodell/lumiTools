import sys, math, csv
import ROOT as r
import lumiTools as l

class CountingExperiment():
    '''
    Class for conducting simple counting experiments for BCM1F, i.e., a simple MC generator
    '''
    def __init__(self, timeStep = 25, nIter = 72, sourceFile= 'histograms/adcHistograms.root', doEff = False, \
            channels = [str(i) for i in range(8)], \
            pList = [0.5 for i in range(8)], \
            sList = [45 for i in range(8)], \
            gList = [(10,12.5) for i in range(8)]
            ):

        # configuration parameters
        self._sourceFile    = r.TFile(sourceFile, 'OPEN')
        self._timeStep      = timeStep
        self._nObjects      = len(channels)
        self._nIter         = nIter
        self._doEff         = doEff

        # initialize generators, counters
        self._rndGen        = r.TRandom3(0)
        self._counter       = 0.

        # Specify the list of outcomes that are being tested (diamond, algo, ?)
        self._channels       = channels 
        self._algos         = ['OR', 'XOR+', 'XOR-', 'AND'] 

        # Channel status booleans
        self._isLive        = l.dict_map(self._channels, [True for i in range(len(channels))]) 
        self._isTriggered   = l.dict_map(self._channels, [[False, False] for i in range(len(channels))]) #[<hit>, <recorded hit>]
        self._isOvershoot   = l.dict_map(self._channels, [False for i in range(len(channels))]) 

        # Saves pulse characteristics for evolving (for now, just TOT)
        self._pulseStatus   = l.dict_map(self._channels, [[0, 0, 0] for i in range(len(channels))]) 

        # Input parameters for detector modelling
        self._fillPattern   = [True for i in range(nIter)]       # list of bools saying whether bunch is colliding or not
        self._probabilities = l.dict_map(self._channels, pList)  # diamond-by-diamond probabilities 
        self._gate          = l.dict_map(self._channels, gList)  # gate width and offset w.r.t. center of hit distribution  
        self._overshoot     = l.dict_map(self._channels, sList)  # saturationn parameters
        #self._tosParams     = l.dict_map(self._channels, osList)

        # Spectra for modelling detector effects
        self._pulseSpectrum  = l.dict_map(self._channels, [])     # 2D hist for producing pulse height and tot 
        self._osSpectrum     = l.dict_map(self._channels, [r.TF1('f_osSpectrum', 'exp([0]*(x - [1]))', -5000, 0) for i in range(len(channels))])
        self._tosSpectrum    = l.dict_map(self._channels, [r.TF1('f_tosSpectrum', '[0] + [1]*atan([2]*(x - [3]))', 0, 2000) for i in range(len(channels))])
        self._bunchSpectrum  = r.TF1('f_bunchSpectrum', 'gaus', 0, 25) # for modelling distribution of colliding hits in time per bunch 
        self._albedoSpectrum = r.TF1('f_albedoSpectrum', 'pol0', 25, 50) # for modelling distribution of albedo hits in time per bunch 

        # Outputs
        self._hits  = l.dict_map(self._channels, [[False for i in range(self._nIter)] for j in range(self._nObjects)])
        self._coins = l.dict_map(self._algos, [[False for i in range(self._nIter)] for j in range(len(self._algos))])


    def initialize_efficiency_histograms(self):
        '''
        Gets histograms that will be samples to produce pulse 
        parameters for diamond hits.
        '''
        for chan in self._channels:
            if chan is not '0':
                self._pulseSpectrum[chan] = self._sourceFile.GetDirectory('ADC_CH'+chan).Get('h2_PulseHeightVsTOT_'+chan).Clone()
            else:
                self._pulseSpectrum[chan] = self._sourceFile.GetDirectory('ADC_CH1').Get('h2_PulseHeightVsTOT_1').Clone()

            self._tosSpectrum[chan].SetParameters(-6144, 1./(1.66e-4), 1./88.28, -14.43)
            self._osSpectrum[chan].SetParameter(0, -0.0006)
        
        self._bunchSpectrum.SetParameters(1, 12.5, 6.)
        self._albedoSpectrum.SetParameters(1, 1)
        

    def set_fill_pattern(self, fillPattern):
        '''
        Takes a list of booleans that represent whether a
        bunch is colliding or not.  Expand this to take a
        fill pattern from the lpc site and smear the hit
        probability based on the variation in current
        '''

        self._fillPattern   = fillPattern
        self._nIter         = len(fillPattern)


    def single_bunch_generator(self, bcid): 
        '''
        Generates diamond-by-diamond hits based on hit probability 
        and whether diamond is live.  Updates isLive and isTriggered
        variables.
        '''
        for chan in self._channels:  
            self._isTriggered[chan] = [False, False]

            if self._isLive[chan] or not self._doEff:  
                if self._fillPattern[bcid]: 
                    '''
                    If channel is live and the crossing is a colliding one, 
                    simulate hits in time gate.  Assume hits are Gaussian, 
                    gate is square.  
                    '''

                    if self._rndGen.Rndm() <= self._probabilities[chan]:
                        self._isTriggered[chan][0] = True
                        hitTime = self._bunchSpectrum.GetRandom()
                        width   = self._gate[chan][0]
                        offset  = self._gate[chan][1]

                        # Get pulse characteristics
                        tosTmp = r.Double(0)
                        phTmp  = r.Double(0)
                        self._pulseSpectrum[chan].GetRandom2(tosTmp, phTmp)
                        self._pulseStatus[chan][0], self._pulseStatus[chan][1] = tosTmp, phTmp

                        #print tosTmp,phTmp

                        if hitTime > (offset - width/2) and hitTime < (offset + width/2):
                            if self._isOvershoot[chan]:
                                if self._pulseStatus[chan][1] > (6 + self._pulseStatus[chan][2]): 
                                    self._isTriggered[chan][1] = True
                            else:
                                self._isTriggered[chan][1] = True

                        self._pulseStatus[chan][0] += hitTime

                elif self._fillPattern[bcid-1]: 
                    '''
                    If channel is live and the previous crossing is a colliding one, 
                    simulate albedo hits.  Assume hits drop off exponentially,
                    but none are recorded as diamond hits (yet...). 
                    '''

                    if self._rndGen.Rndm() <= 0.03*self._probabilities[chan]:
                        # Get pulse characteristics
                        tosTmp = r.Double(0)
                        phTmp  = r.Double(0)
                        self._pulseSpectrum[chan].GetRandom2(tosTmp, phTmp)
                        self._pulseStatus[chan][0], self._pulseStatus[chan][1] = tosTmp, phTmp

                        self._isTriggered[chan][0] = True
                        hitTime = self._albedoSpectrum.GetRandom()

                        self._pulseStatus[chan][0] += hitTime


        # update pulse characteristics
        self.update_pulse()

        # For debugging
        # print self._isTriggered
        # print self._isLive
        # print self._pulseStatus['1']
        # print '\n'
                
        self._counter += self._timeStep


    def update_pulse(self):
        '''
        Regulates deadtime counters for individual chans.  This will
        account for the characteristic decay spectrum of the deadtime
        and will have to include something special for overshoot of baseline.
        '''

        for chan in self._channels:

            if self._isTriggered[chan][0]:
                overTime, pulseHeight = self._pulseStatus[chan][0], self._pulseStatus[chan][1] 
                overTime *= 1.

                if overTime > 110 and chan is '7': 
                    self._isOvershoot[chan] = True
                    osTime = self._tosSpectrum[chan].Eval(overTime+20)
                    self._pulseStatus[chan][0] = overTime + 0.20*osTime 

                    expTime = 2.5*osTime
                    self._osSpectrum[chan].SetParameter(1, 0.9*expTime)
                    self._pulseStatus[chan][2] = expTime
            else:
                self._pulseStatus[chan][1] = 0

            if self._pulseStatus[chan][0] > self._timeStep:
                self._pulseStatus[chan][0] -= self._timeStep 
                self._isLive[chan] = False
            else:
                self._pulseStatus[chan][0] = 0
                self._isLive[chan] = True

                if self._pulseStatus[chan][2] > self._timeStep:
                    self._pulseStatus[chan][2] -= self._timeStep
                else:
                    self._pulseStatus[chan][2] = 0


    def make_coincidences(self, bx):
        '''
        Form coincidences of hist on a per bunch crossing
        basis.  For now, +/-z OR, AND, XOR+, XOR-
        '''

        if (True in [self._isTriggered[str(i)][1] for i in range(4)]) \
            and (True in [self._isTriggered[str(i)][1] for i in range(4,8)]):
            self._coins['XOR-'][bx] \
            = self._coins['XOR+'][bx] = False
            self._coins['AND'][bx]  \
            = self._coins['OR'][bx]   = True

        elif (True in [self._isTriggered[str(i)][1] for i in range(4)]) \
            and (True not in [self._isTriggered[str(i)][1] for i in range(4,8)]):
            self._coins['AND'][bx] \
            = self._coins['XOR-'][bx] = False
            self._coins['XOR+'][bx] \
            = self._coins['OR'][bx]   = True

        elif (True not in [self._isTriggered[str(i)][1] for i in range(4)]) \
            and (True in [self._isTriggered[str(i)][1] for i in range(4,8)]):
            self._coins['AND'][bx] \
            = self._coins['XOR+'][bx] = False
            self._coins['XOR-'][bx] \
            = self._coins['OR'][bx]   = True

        else:
            self._coins['AND'][bx] \
            = self._coins['XOR+'][bx] \
            = self._coins['XOR-'][bx] \
            = self._coins['OR'][bx] = False 



    def run_bunch_by_bunch(self):
        '''An experiment for bunch-by-bunch efficiency estimates'''

        for i in range(self._nIter):
            self.single_bunch_generator(i)

            for chan in self._channels:
                if self._isTriggered[chan][1]:
                    self._hits[chan][i] = True
                else: 
                    self._hits[chan][i] = False

            self.make_coincidences(i)

        for chan in self._channels:
            self._isLive[chan] = True
            self._pulseStatus[chan][0] = 0
            self._pulseStatus[chan][1] = 0
            
        return (self._hits, self._coins)
