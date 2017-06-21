import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
#from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt

def getPLV():
    
    '''
    compute phase lag variance
    PLV plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    '''

    try:    

        eegData = loadmat('EEG181.mat')
        
        eegFS = 250
        
        eeg1 = eegData[:14]
        eeg2 = eegData[:15]        

        eeg1 = [ x[50*eegFS:80*eegFS] for x in eeg1] 
        eeg1 = [ x[50*eegFS:80*eegFS] for x in eeg2] 

        eeg1 = eeg1 - np.mean(eeg1, axis = 0)        
        eeg2 = eeg2 - np.mean(eeg2, axis = 0)        
        
        eegFS = 2000 # %sampling rate
            
        # create filters
        freqs = np.linspace( 5,5,80)
        
        filters = {} #cell(length(freqs),2);
        
        bw = 4 # %bandwidth
        nyq = eegFS/2 #;
        attendB = 40 # %attenuation
        attenHz = 2 # %transition band
        #gamma freqs
        
        for fI in range (len(freqs)):
        
            lcut = freqs[fI]-bw/2
            hcut = freqs[fI]+bw/2
            
            Fstop1 = (lcut - attenHz)  / nyq  
            Fpass1 = lcut  / nyq 
            Fpass2 = hcut / nyq
            Fstop2 = (hcut + attenHz) / nyq
            Astop1 = attendB
            Apass  = 1
            Astop2 = attendB
            h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2);
            Hd = kaiser(h, 'kaiserwin')
            b = Hd.Numerator
            
            filters[fI] = ()
            filters[fI](0) = b
            
            #group delay
            [a,f] = group_delay(b,1,nyq,eegFS)
            k = f >= lcut & f <= hcut
            gd = np.floor(np.mean(a[k]))
            filters[fI][1] = gd
        
        PLVs = np.zeros(len(freqs))
            
        #filters
        for fI in range ( len(filters) ) :
            b = filters[fI][0]
            gd = filters[fI][1]
        
            ph1 = lfilter(b,1,eeg1)
            ph1 = [ph1[gd+1:end], np.zeros(1,gd)]
            ph1 = np.angle(hilbert(ph1))
        
            ph2 = lfilter(b,1,eeg2);
            ph2 = [ph2[gd+1:end] ,zeros(1,gd)]
            ph2 = np.angle(hilbert(ph2))
        
            # phase diff
            phiD = ph1 - ph2    
        
            # PLV
            PLVs[fI] = (1/len(phiD))*(abs(sum(np.exp(1i*phiD)))) 
        
        plot(freqs,PLVs)
        
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return             