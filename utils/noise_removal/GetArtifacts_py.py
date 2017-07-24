#get artifacts, store results in file
import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from scipy.signal import firwin,firwin2,lfilter, kaiser
import matplotlib.pyplot as plt
import pandas as pd

def GetNoiseHF(eeg,eegFS,lenHF,freqHF,winLenVarHF,winLenSmoothHF,thHF):

        '''
        '''
        
        #normalize
        if max(eeg) != 0:
                eeg = np.array(eeg)/max(eeg)
              
        #filter HF
        #window = kaiser(51, beta=14)
        winLen = lenHF
        # to be in conformance with MATLAB ( check documentation of fir1 in MATLAB, it automatically adds 1 for even sized windows)
        if lenHF % 2 ==0:
                winLen = lenHF + 1
        h = firwin(winLen, freqHF/(eegFS/2), pass_zero=False)  
        h = np.array(h)* 1.0/sum(h)
        print (h[19:39])
        plt.plot(h)
        plt.show()
        h = [-1*a for a in h]

        midPoint = math.floor(lenHF/2)
        h[midPoint] = h[midPoint] + 1
        print (h)
        #h = firwin(lenHF,freqHF/(eegFS/2),'high')
        sF = lfilter(h,1,eeg)
        #print (len(sF))
        #print (len(sF[midPoint:]))

        sF = list(sF[midPoint:]) + [0]*(midPoint)
        
        #print ( sF[:20])
                
        ###get variance
        #sf = pd.Series(sF).olling_stdr(winLenVarHF).std()
        rstd = pd.rolling_std(pd.Series(sF), window = int(winLenVarHF))
        rvar = [x*x for x in rstd]
        #print ( rvar ) 
        
        #square
        sF = np.square(sF)
                
        #smooth
        h = (1/int(winLenSmoothHF)) * np.ones(int(winLenSmoothHF))
        sF = lfilter(h,1,sF)
        midPoint = math.floor(winLenSmoothHF/2)
        sF = list(sF[midPoint:]) + [0]*(midPoint)
                
        #1=good signal
        signalOK = [x for x in sF if x < thHF ]
        #print ( signalOK)


def getNoiseSat(eeg,maxLengthCrossing):
    
    '''
    Input: 
    Output: 
    '''

    try:    

        #get signal difference
        dif = np.diff(eeg)

        #print ( " dif = " + str(dif[:10]) )

        k = [1 if x ==0 else 0 for x in dif]
        
        #print (" k " + str(k))
        #add first sample (diff takes one sample out)
        k = [k[0]] + k
        #print (" k " + str(k[:30]))        
        #skip extremely short saturations (1 value is 0), find 010
        
        k[0] = 0
        k[-1] = 0

        kUP = [i for i, x in enumerate(k) if x == 1 and k[i-1] == 0 ]
        kUP = np.array(kUP) + 1
        kDOWN = [i for i, x in enumerate(k) if x == 0 and k[i-1] == 1 ]

        d = kDOWN-kUP

        print ( d[:20])   
        kd = [ i for i,x in enumerate(k) if x < 2 ]
        #print ( kd[:20])   

        for i, x in enumerate ( kd ):
                #print ( str(i) + " : " + str(len(k)) ) 
                #print (kUP[i])
                #print (kDOWN[i])
                print ( str(i) + " : " + str(len(kUP)) + " : " + str(len(kDOWN) ) ) 
                if i >= len(kUP) or i >= len(kDOWN):
                        break
                k[kUP[i]:kDOWN[i]] = [0]*(kDOWN[i] - kUP[i])
            
        print ( k[:20] ) 
        
        k = 1- np.array(k)
        print ( str(len(k) ) ) 

        #overshoots allowed
        if maxLengthCrossing != 0:                
                #first find overshoots and allow short ones
                k[0] = 1
                k[-1] = 1
    
                k[0] = 0
                k[-1] = 0
        
                kDOWN = [i for i, x in enumerate(k) if x == 1 and k[i-1] == 0 ]
                kDOWN = np.array(kDOWN) + 1
                
                kUP = [i for i, x in enumerate(k) if x == 0 and k[i-1] == 1 ]
        
                d = kUP-kDOWN
                
                kd = [x for x in d if x > maxLengthCrossing] # %only select large crossings
                kDOWN = [ kDOWN[i] for i,x in enumerate(kd) if x == 1 ] 
                kUP = [ kUP[i] for i,x in enumerate(kd) if x == 1 ] 
        
                #create new series and find long enough windows
                k = []
                for i,x in enumerate(kDOWN):
                                    
                        if i >= len(kUP) or i >= len(kDOWN):
                                break
                        k[kUP[i]:kDOWN[i]] = [0]*(kDOWN[i] - kUP[i])
                
        signalOK = np.array(k)

    except:
        traceback.print_exc(file=sys.stdout)    
        
    return signalOK      

def getArtifacts():
    
    '''

    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    '''

    try:    

        # input directory
        # drIn = '/mnt/kant/dinod/ANALYSIS/basma/group1/MAT/';
        # output directory
        # drOut = '/mnt/kant/dinod/ANALYSIS/basma/group1/ART_border025dif/';
        
        # if exist(drOut,'dir') == 0; mkdir(drOut); end;
        
        doSaturation = 1 # DO saturation test (set to 1)
        doSmallSig = 1 # DO slow signal test
        doHF = 1 # DO high frequency artifact test
        
        eegFS = 250 # sampling rate
        
        eegData = loadmat('EEG181.mat')
        
        eegData = eegData["eegData"]
                
        eeg = eegData[14]
        
        eeg = eeg[50 * eegFS - 1 : 80 * eegFS - 1 ]         
        
        # saturation 
        maxLengthCrossing = 5 # %maximum length of crossing in samples (allow spikes to overshoot)
        
        # slow signal
        thVar = 1e-4 # variance threshold - DECREASE for LOWER sensitivity
        thVarMin = 1e-2 # minimum variance - DECREASE for LOWER sensitivity
        thVarMax = 2 # maximum variance
        winLenVar = 0.05*eegFS # minimum variance window under threshold
        winLenVarSlide = 0.05*eegFS # sliding variance window
        
        # high frequency
        lenHF = 50 # %length of HF filter [samples]
        freqHF = 100 # %frequency of HF filter [Hz]
        winLenVarHF = 0.1*eegFS # %variance window [samples]
        winLenSmoothHF = 0.2*eegFS # %smoothing window after HF [samples]
        thHF = 1e-5 # %threshold for HF signal - INCREASE for LOWER sensitivity
        
        # windows computation
        minLen = 1*eegFS #   %minimum required length of good signal in seconds
        winSafe = 0.25*eegFS #  %safe window around saturations to exclude
                
        sigSat = []
        sigSmall = []
        sigHF = []
        
        #get saturations (non-changing signal)
        #if doSaturation == 1:
            #sigSat = getNoiseSat(eeg,maxLengthCrossing)
        
        # get HF noise
        if doHF == 1:
            sigHF = GetNoiseHF(eeg,eegFS,lenHF,freqHF,winLenVarHF,winLenSmoothHF,thHF)
        
        # merge all three signals
        
        #k = []
        
        #if doSaturation == 1:
            #k = [x for x,y in zip(k, sigSat) if y == 1] 
        
        #if doHF == 1:
            #k = [x for x,y in zip(k, sigHF) if y == 1] 

        # compute windows
        
        #[sigOKCH, winOKCH] = GetWindows(k,minLen,winSafe)
        
        ##pltsubplot(5,1,1);
        #plt.plot(eeg,'k');
        #plt.title('eeg');
        ##axis tight;
        
        ##subplot(5,1,2);
        #plt.plot(sigSat,'r');
        #plt.title('saturations');
        ##axis tight;
        
        ##subplot(5,1,4);
        #plt.plot(sigHF,'r');
        #plt.title('hf signal');
        ##axis tight;
        
        ##subplot(5,1,5); hold on;
        #plt.plot(k,'k*');
        #plt.plot(sigOKCH,'r');
        #plt.title('good signal');
        ##axis tight;
        ##%}
        
        ##store data
        #signalOK[chI] = sigOKCH
        #windowsOK[chI] = winOKCH
        
        #save data
        #save([drOut files(fI).name], 'signalOK','windowsOK');

    except:
        traceback.print_exc(file=sys.stdout)    
    return      

getArtifacts()