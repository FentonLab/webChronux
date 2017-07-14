import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
import scipy.signal as signal
from scipy.signal import *
#from scipy.signal import mfreqz, impz, kaiser, group_delay, hilbert
import matplotlib.pyplot as plt

from numpy import cos, sin, pi, absolute, arange
from numpy.random import normal
from scipy.signal import kaiserord, lfilter, firwin, freqz
#from pylab import figure, clf, plot, xlabel, ylabel, xlim, ylim, title, grid, axes, show
from math import sqrt, log10

def getMI(amp, phase, edges):

    '''
    compute modulation index as shown in Tort 2010
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
         - amp: Amplitude
         - phase: Phase
         - edges: Edges
    Output: 
         - MI
    '''

    MI = ''

    try:
        
        # get histogram for phases
        h = []
        
        for hi in range ( len(edges) - 1 ):
            
            k = [ 1 if x >= edges[hi] and x < edges[hi+1] else 0 for x in phase]
            #print ( " k = " + str(k) ) 
            
            if sum(k) > 0: # only if there are some values, otherwise keep zero
                print (sum([x for i,x in enumerate(amp) if k[i] == 1 ]))
                
                h.append( np.mean([x for i,x in enumerate(amp) if k[i] == 1 ])) # mean of amplitudes at that phase bin
            else:
                h.append(0)
                
        #print ( " h = " + str(h) )
    
        ## fix last value
        #k = [x == y for x,y in zip ( phase, edges [end])] 
        #if not all([x == 0 for x in k]):
            #h[end] = h[end] + np.mean(amp[k])
        
        ## convert to probability
        h = h / sum(h)
                
        #print ( " h = " + str(h) )
        
        ## replace zeros by eps
        #k = h == 0
        #h[k] = eps
        
        ## calculate modulation index
        Hp = -1 * sum( [x*np.log(x) for x in h]) # entropy of the histogram
        #print ( " Hp = " + str(Hp) )

        D = np.log(len(h)) - Hp
        #print ( " D = " + str(D) )
 
        MI = D / np.log(len(h))
        #print ( " MI = " + str(MI) )
    
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return MI

def test():
        
    '''
    compute phase amplitude coupling, as in Vinnick
    PAC plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    '''

    try:
        
        eegData = loadmat('../EEG181.mat')
        
        eegFS = 250 # sampling frequency
        
        eegData = eegData["eegData"]
                
        eeg = eegData[14]
        
        print ( eeg[:5])
        
        eeg = eeg[50 * eegFS - 1 : 80 * eegFS - 1 ]   
        
        eeg = eeg - np.mean(eeg, axis = 0)         
        
        # edges for phase
        edges = np.linspace(-math.pi,math.pi,21) 
        
        #edges = np.array(list(edges).append(math.pi))
        
        lcut = 9.0 
        hcut = 13.0
        
        #----------------------------------------------------------
        # Set the sampling frequency to 100. Sample period is 0.01. 
        # A frequency of 1 will have 100 samples per cycle, and 
        # therefore have 4 cycles in the 400 samples.
        
        sample_rate = 250.0
        nyq = sample_rate / 2.0
        width = 2.0/nyq # pass to stop transition width
        
        # The desired attenuation in the stop band, in dB.
        ripple_db = 40.0
        
        # Compute the order and Kaiser parameter for the FIR filter.
        N, beta = kaiserord(ripple_db, width)
        
        #print ('N = ',N, 'beta = kaiser param = ', beta)
        
        # The cutoff frequency of the filter.
        cutoff_hz = lcut
        
        # Use firwin with a Kaiser window to create a lowpass FIR filter.
        hpftaps = firwin(N, cutoff_hz/nyq, window=('kaiser', beta))    
        
        #print (hpftaps[:10])
        
        #----------------------------------------------------------
        # now create the taps for a high pass filter.
        # by multiplying tap coefficients by -1 and
        # add 1 to the centre tap ( must be even order filter)
        
        hpftaps = [-1*a for a in hpftaps]
        print ( len(hpftaps))
        midPoint = int(np.round(len(hpftaps)/2))
        if midPoint % 2 != 0:
            midPoint = midPoint -1
        hpftaps[midPoint] = hpftaps[midPoint] + 1
        
        #----------------------------------------------------------
        # Now calculate the tap weights for a lowpass filter at say 15hz
        
        cutoff_hz = hcut
        lpftaps = firwin(N, cutoff_hz/nyq, window=('kaiser', beta))
        
        # Subtract 1 from lpf centre tap for gain adjust for hpf + lpf
        lpftaps[midPoint] = lpftaps[midPoint] - 1
        
        taps = [sum(pair) for pair in zip(hpftaps, lpftaps)]

        denom = [0]*len(taps)
        denom[0] = 1
        #denom[-1] = 1
        
        [a,f] = group_delay( [ taps , denom ] , int(nyq))  
        
        #print ( " num taps = " + str(len(taps)))
        #print ( " taps [:10] = " + str(taps[:10]))        
        
        #print ( " a = " + str(a) ) 
        #print ( " f = " + str(len(f)) ) 
        
        bAlpha = taps           
        
        bAlphax = np.array(a) * sample_rate/(2*math.pi)
        #print ( " bAlpha **** = " + str(bAlpha) )         
        
        k = [f[i] for i,m in enumerate(bAlphax) if m >= 9 and m <= 13]
        #print ( " k = " + str(k) ) 
        
        gdAlpha = math.floor(np.mean(k))
        #print ( "gdAlpha = " + str(gdAlpha ) )
        
        fGamma = np.arange(20,101,5)
        
        #print ( " fGamma = " + str(fGamma) )

        bw = 20 #bandwidth
        attendB = ripple_db #attenuation
        attenHz = 4 #transition band
        
        filtersGamma = {}

        for fI in range ( len(fGamma) ): 
        
            lcut = fGamma[fI]-bw/2
            hcut = fGamma[fI]+bw/2
            
            Fstop1 = (lcut - attenHz)  / nyq  
            Fpass1 = lcut  / nyq 
            Fpass2 = hcut / nyq
            
            Fstop2 = (hcut + attenHz) / nyq
            Astop1 = attendB
            Apass  = 1
            Astop2 = attendB
            
            #################
            hpftaps = firwin(N, lcut/nyq, window=('kaiser', beta))    
            hpftaps = [-1*a for a in hpftaps]

            midPoint = int(np.round(len(hpftaps)/2))
            if midPoint % 2 != 0:
                midPoint = midPoint -1

            hpftaps[midPoint] = hpftaps[midPoint] + 1

            lpftaps = firwin(N, hcut/nyq, window=('kaiser', beta))
            lpftaps[midPoint] = lpftaps[midPoint] - 1
            
            taps = [sum(pair) for pair in zip(hpftaps, lpftaps)]
    
            denom = [0]*len(taps)
            denom[0] = 1
            
            [a,f] = group_delay( [ taps , denom ] , int(nyq))  

            x = np.array(a) * sample_rate/(2*math.pi)

            k = [f[i] for i,m in enumerate(x) if m >= 9 and m <= 13]
            val = math.floor(np.mean(k))
            #################

            filtersGamma[fI] = []
            
            filtersGamma[fI].append (taps)
            filtersGamma[fI].append (val)

        MIs = []
        
        # phase (alpha)
        
        #print ( " bAlpha = " + str(bAlpha[:20] ) ) 
        #print ( " eeg = " + str(eeg[:20] ) ) 

        ph = lfilter(bAlpha,1,eeg)
        
        #print ( " ph before angle = " + str(ph[:10] ) ) 
        
        ph = np.append(ph[gdAlpha+1:], np.zeros(gdAlpha))
        ph = np.angle(hilbert(ph))
        
        #print ( " ph after angle = " + str(ph[:10] ) ) 
        
        ##amplitude for gamma range + MI
        for key, value in filtersGamma.items(): 
            
            #amplitude
            amp = lfilter(value[0],1,eeg)
            gd = value[1]
            
            amp = np.append(amp[gd+1:], np.zeros(gd))
            amp = abs(hilbert(amp))
            
            #compute raw modulation index
            #print (" amp = " + str(amp[:20]) ) 
            #print (" ph = " + str(ph[:20]) ) 
            #print (" edges = " + str(edges) ) 
            
            MI = getMI(amp,ph,edges)
            MIs.append(MI)
            #break
        
        print ( " MIs = " + str(MIs))
        plt.plot(fGamma,MIs)
        plt.show()

        #print ( [str(k) + ":" + str(v[1]) + "-" + str(v[0][:5]) for k,v in filtersGamma.items() ]  )
        
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return    

test()    
