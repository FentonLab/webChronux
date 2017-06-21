import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt

def getPAC():

    '''
    compute phase amplitude coupling, as in Vinnick
    PAC plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    '''

    try:

        x = loadmat('EEG181.mat');
        
        eegFS = x['eegFS'][0][0] 
        
        eegData = x['eegData']
        
        eeg = eegData[:15]
        
        eeg = [ x[50*eegFS:80*eegFS] for x in eeg] 

        eeg = eeg - np.mean(eeg, axis = 0)
        
        # edges for phase
        edges = np.linspace(-math.pi,math.pi,21)
        
        eegFS = 2000 #sampling rate
        
        winLen = len(eeg) #window length
        
        # create alpha filter
        
        lcut = 9
        hcut = 13
        attendB = 40
        attenHz = 2
        nyq = eegFS/2
        
        Fstop1 = (lcut - attenHz)  / nyq  
        Fpass1 = lcut  / nyq 
        Fpass2 = hcut / nyq
        Fstop2 = (hcut + attenHz) / nyq
        Astop1 = attendB
        Apass  = 1
        Astop2 = attendB
        
        h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1,Fpass2, Fstop2, Astop1, Apass, Astop2)
     
        Hd = kaiser(h, 'kaiserwin')
        
        bAlpha = Hd.Numerator
        #group delay
        [a,f] = group_delay(bAlpha,1,nyq,eegFS)
        
        k = f >= lcut & f <= hcut
        
        gdAlpha = np.floor(mean(a[k]))
           
        #create multiple gamma filters
        fGamma = np.linspace (20, 100, 5) 

        filtersGamma = {} # cell(length(fGamma),2);

        bw = 20 # %bandwidth
        attendB = 40 # attenuation
        attenHz = 4 # transition band
        
        #gamma freqs
        
        for fI in range ( len(fGamma) ): 
        
            lcut = fGamma(fI)-bw/2
            hcut = fGamma(fI)+bw/2
            
            Fstop1 = (lcut - attenHz)  / nyq  
            Fpass1 = lcut  / nyq 
            Fpass2 = hcut / nyq
            
            Fstop2 = (hcut + attenHz) / nyq
            Astop1 = attendB
            Apass  = 1
            Astop2 = attendB
            
            h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1,Fpass2, Fstop2, Astop1, Apass, Astop2)
            Hd = kaiser(h, 'kaiserwin')
            
            b = Hd.Numerator

            filtersGamma[fI] = ()
            
            filtersGamma[fI][0] = b
            
            #group delay
            [a,f] = group_delay(b,1,nyq,eegFS)
            
            k = f >= lcut & f <= hcut
            gd = np.floor(np.mean(a[k]))
            filtersGamma[fI](1) = gd
        
        MIs = np.zeros(len(fGamma))
        
        # phase (alpha)
        ph = lfilter(bAlpha,1,eeg)
        ph = [ph[gdAlpha+1:end] , np.zeros(gdAlpha)]
        ph = np.angle(hilbert(ph))
        
        #amplitude for gamma range + MI
        for fI in range(len(fGamma)): 
            
            b = filtersGamma[fI][0]
            gd = filtersGamma[fI][1]
            
            #amplitude
            amp = lfilter(b,1,eeg)
            
            amp = [ amp [gd+1:end] , np.zeros(gd)]
            amp = abs(hilbert(amp))
            
            #compute raw modulation index
            MI = getMI(amp,ph,edges)
        
            MIs(fI) = MI;
        end
        
        plt.plot(fGamma,MIs)
        
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return    
