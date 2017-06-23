import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
#from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt
import scipy.signal.welch

# get frequency grid
def getPSD(lowerFrequency, upperFrequency, samplingFrequency, padding):

    WINDOW_SIZE = 3
    
    try:    

        #eegData = loadmat('EEG181.mat')
        
        #eegFS = 250
        
        #eegData = eegData["eegData"]
        
        #eeg = eegData[14] [50*eegFS:80*eegFS]

        #eeg = eeg - np.mean(eeg, axis = 0)    
        
        #Twindow = 1 # length of subwindow
        #Nwindow = Twindow * eegFS # length of window in samples
        #Noverlap = round ( Nwindow /2) # 50% overlap
        #NFFT = 2**11 # %resolution
        
        #[~,f] = pwelch(rand(1,length(eeg)), Nwindow, Noverlap, NFFT, eegFS) # %get frequency vector (can be done other ways)
         
        f = welch(x, fs=1.0, window='hanning', noverlap=None, nfft=None, detrend='constant')        
        
        #pwelch
        Pxx = pwelch(eeg, Nwindow, Noverlap, NFFT, eegFS)
        Pxx = 20*log10(Pxx) #convert to dB
        
        plot(f,Pxx)

    except:
        traceback.print_exc(file=sys.stdout)    
        
    return   
