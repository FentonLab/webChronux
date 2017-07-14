#get artifacts, store results in file
import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
#from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt

def getNoiseSat(eeg,maxLengthCrossing):
    
    '''
    Input: 
    Output: 
    '''

    try:    

        #get signal difference
        dif = np.diff(eeg)

        print ( " dif = " + str(dif[:10]) )

        k = [1 if x ==0 else 0 for x in dif]
        
        #print (" k " + str(k))
        #add first sample (diff takes one sample out)
        k = [k[0]] + k
        print (" k " + str(k[:30]))        
        ##skip extremely short saturations (1 value is 0), find 010
        k[0] = 0
        k[-1] = 0

        kUP = [i for i, x in enumerate(k) if x == 1 and k[i-1] == 0 ]
        #kUP = kUP + 1;
        kDOWN = [i for i, x in enumerate(k) if x == 0 and k[i-1] == 1 ]

        #d = kDOWN-kUP
        #kd = find(d < 2) # allow two samples being same
        #for i in range (len(kd)):
            #ind = kd(i)
            #k(kUP(ind):kDOWN(ind)) = 0
        
        #k = ~k

        #overshoots allowed
        #if maxLengthCrossing != 0:                
            ##first find overshoots and allow short ones
            #k(1) = 1; k(end) = 1
            #kDOWN = find(k(1:end-1) == 1 & k(2:end) == 0)
            #kDOWN = kDOWN + 1
            #kUP = find(k(1:end-1) == 0 & k(2:end) == 1)
            #d = kUP - kDOWN + 1 # %length of overshoots
            #kd = d > maxLengthCrossing # %only select large crossings
            #kDOWN = kDOWN(kd)
            #kUP = kUP(kd)
        
            ##create new series and find long enough windows
            #k = ones(1,length(k))
            #for i in range ( len(kDOWN)):
                #st = kDOWN(i)
                #ed = kUP(i)
                #k(st:ed) = 0
                
        #signalOK = logical(k);

    except:
        traceback.print_exc(file=sys.stdout)    
        
    return       

def getArtifacts():
    
    '''
    compute phase lag variance
    PLV plotted at the end
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
        doHF = 0 # DO high frequency artifact test
        
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
        freqHF = 250 # %frequency of HF filter [Hz]
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
        if doSaturation == 1:
            sigSat = getNoiseSat(eeg,maxLengthCrossing)
        
        ## get HF noise
        #if doHF == 1:
            #sigHF = GetNoiseHF(eeg,eegFS,lenHF,freqHF,winLenVarHF,winLenSmoothHF,thHF)
        
        # merge all three signals
        
        k = []
        
        #if doSaturation == 1:
            #k = k & sigSat
        
        #if doHF == 1:
            #k = k & sigHF

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