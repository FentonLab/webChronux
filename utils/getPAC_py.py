import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from scipy.signal import kaiser, group_delay, hilbert, kaiserord, firwin, lfilter
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

        eegData = loadmat('EEG181.mat')
        
        eegFS = 250 # sampling frequency
        
        eegData = eegData["eegData"]
                
        eeg = eegData[14]
        
        print ( eeg[:5])
        
        eeg = eeg[50 * eegFS - 1 : 80 * eegFS - 1 ]   
        
        eeg = eeg - np.mean(eeg, axis = 0)         
        
        # edges for phase
        edges = np.linspace(-math.pi,math.pi,21)
        
        print ( " eeg = " + str(eeg[:10]))
        print ( " edges = " + str(len(edges)))
        
        #eegFS = 2000 #sampling rate
        
        #winLen = len(eeg) #window length
        
        ## create alpha filter
        
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
        
        print ( " nyq = " + str(nyq) )
        print ( " Fstop1 = " + str(Fstop1) )
        print ( " Fpass1 = " + str(Fpass1) )
        print ( " Fpass2 = " + str(Fpass2) )
        print ( " Astop1 = " + str(Astop1) )
        
        ####################################
        width = 2.0/nyq # pass to stop transition width
        
        # The desired attenuation in the stop band, in dB.
        ripple_db = 40.0
        
        # Compute the order and Kaiser parameter for the FIR filter.
        N, beta = kaiserord(ripple_db, width)
        
        #print ('N = ',N, 'beta = kaiser param = ', beta)
        
        # The cutoff frequency of the filter.
        cutoff_hz = 9.0
        
        # Use firwin with a Kaiser window to create a lowpass FIR filter.
        hpftaps = firwin(N, cutoff_hz/nyq, window=('kaiser', beta))
        
        #----------------------------------------------------------
        # now create the taps for a high pass filter.
        # by multiplying tap coefficients by -1 and
        # add 1 to the centre tap ( must be even order filter)
        
        hpftaps = [-1*a for a in hpftaps]
        hpftaps[33] = hpftaps[33] + 1
        
        #----------------------------------------------------------
        # Now calculate the tap weights for a lowpass filter at say 15hz
        
        cutoff_hz = 13.0
        lpftaps = firwin(N, cutoff_hz/nyq, window=('kaiser', beta))
        
        # Subtract 1 from lpf centre tap for gain adjust for hpf + lpf
        lpftaps[33] = lpftaps[33] - 1
        
        #----------------------------------------------------------
        # Now add the lpf and hpf coefficients to form the bpf.
        
        taps = [sum(pair) for pair in zip(hpftaps, lpftaps)]
        
        print ( " num taps = " + str(len(taps)))
        print ( " taps [:10] = " + str(taps[:10]))
        
        #----------------------------------------------------------
        # Use lfilter to filter test signal with Bandpass filter.
        filtered_x = lfilter(taps, 1.0, eeg)   
        print ( " filtered_x = " + str(filtered_x[:10]))        
        ####################################
        
        #h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1,Fpass2, Fstop2, Astop1, Apass, Astop2)
     
        #Hd = kaiser(h, 'kaiserwin')
        
        #bAlpha = Hd.Numerator
        ##group delay
        #[a,f] = group_delay(bAlpha,1,nyq,eegFS)
        
        #k = f >= lcut & f <= hcut
        
        #gdAlpha = np.floor(mean(a[k]))
           
        ##create multiple gamma filters
        #fGamma = np.linspace (20, 100, 5) 

        #filtersGamma = {} # cell(length(fGamma),2);

        #bw = 20 # %bandwidth
        #attendB = 40 # attenuation
        #attenHz = 4 # transition band
        
        ##gamma freqs
        
        #for fI in range ( len(fGamma) ): 
        
            #lcut = fGamma(fI)-bw/2
            #hcut = fGamma(fI)+bw/2
            
            #Fstop1 = (lcut - attenHz)  / nyq  
            #Fpass1 = lcut  / nyq 
            #Fpass2 = hcut / nyq
            
            #Fstop2 = (hcut + attenHz) / nyq
            #Astop1 = attendB
            #Apass  = 1
            #Astop2 = attendB
            
            #h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1,Fpass2, Fstop2, Astop1, Apass, Astop2)
            #Hd = kaiser(h, 'kaiserwin')
            
            #b = Hd.Numerator

            #filtersGamma[fI] = ()
            
            #filtersGamma[fI][0] = b
            
            ##group delay
            #[a,f] = group_delay(b,1,nyq,eegFS)
            
            #k = f >= lcut & f <= hcut
            #gd = np.floor(np.mean(a[k]))
            #filtersGamma[fI](1) = gd
        
        #MIs = []
        
        ## phase (alpha)
        #ph = lfilter(bAlpha,1,eeg)
        #ph = [ph[gdAlpha+1:end] , np.zeros(gdAlpha)]
        #ph = np.angle(hilbert(ph))
        
        ##amplitude for gamma range + MI
        #for fI in range(len(fGamma)): 
            
            #b = filtersGamma[fI][0]
            #gd = filtersGamma[fI][1]
            
            ##amplitude
            #amp = lfilter(b,1,eeg)
            
            #amp = [ amp [gd+1:end] , np.zeros(gd)]
            #amp = abs(hilbert(amp))
            
            ##compute raw modulation index
            #MI = getMI(amp,ph,edges)
        
            #MIs(fI) = MI;
        
        #plt.plot(fGamma,MIs)
        
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return    

getPAC()
