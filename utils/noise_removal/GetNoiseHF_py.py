#get artifacts, store results in file

import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math

#from scipy.signal import kaiser, group_delay, hilbert, angle

import matplotlib.pyplot as plt

def GetNoiseHF(eeg,eegFS,lenHF,freqHF,winLenVarHF,winLenSmoothHF,thHF):

        #%{
        #%testing:
        #load('data samples/sat1.mat');
        #eeg = x;
        #lenHF = 50;               
        #freqHF = 250;             
        #winLenVarHF = eegFS/10;   
        #winLenSmoothHF = eegFS/5; 
        #thHF = 1e-6;              
        #%}
        
        #normalize
        eeg = eeg/max(eeg)
              
        #filter HF
        h = fir1(lenHF,freqHF/(eegFS/2),'high')
        sF = lfilter(h,1,eeg);
        sF = [sF(lenHF/2+1:end) ,zeros(1,lenHF/2)]
                
        #get variance
        sF = rvar(double(sF),winLenVarHF,length(sF))
                
        #square
        sF = sF**2
                
        #smooth
        h = (1/winLenSmoothHF)*ones(1,winLenSmoothHF)
        sF = filter(h,1,sF)
        sF = [sF(winLenSmoothHF/2+1:end) ,zeros(1,winLenSmoothHF/2)]
                
        #1=good signal
        signalOK = sF < thHF;
             
        #%{
        #subplot(3,1,1);
        #plot(eeg,'k');
        #subplot(3,1,2);
        #plot(sF,'r');
        #subplot(3,1,3);
        #plot(signalOK,'r');
        #%}
