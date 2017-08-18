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
                
        eeg1 = eegData[14]
        
        print ( eeg1[:5])
        
        eeg1 = eeg1[50 * eegFS - 1 : 80 * eegFS - 1 ]   
        
        eeg1 = eeg1 - np.mean(eeg1, axis = 0)         

        eeg2 = eegData[15]
        
        print ( eeg2[:5])
        
        eeg2 = eeg2[50 * eegFS - 1 : 80 * eegFS - 1 ]   
        
        eeg2 = eeg2 - np.mean(eeg2, axis = 0)         
        
        # edges for phase
        edges = np.linspace(-math.pi,math.pi,21) 
        
        x = np.linspace(0,200,5000);
        
        s1 = np.sin( x * pi / 180 )
        
        s2 = np.sin( x * pi / 18 )
        
        s2 = np.array ( [x*y for x,y in zip(s1,s2) ] ) 
        
        s3 = s1 + s2 
        
        eeg = s3
        
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

        winLen = N
        # to be in conformance with MATLAB ( check documentation of fir1 in MATLAB, it automatically adds 1 for even sized windows)
        if N % 2 ==0:
                winLen = N + 1        
        
        # Use firwin with a Kaiser window to create a lowpass FIR filter.
        hpftaps = firwin(winLen, cutoff_hz/nyq, window=('kaiser', beta), pass_zero=False)    
        
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

        lpftaps = firwin(winLen, cutoff_hz/nyq, window=('kaiser', beta), pass_zero=False)
        
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

        PLVs = []
        
        # phase (alpha)
        
        #print ( " bAlpha = " + str(bAlpha[:20] ) ) 
        #print ( " eeg = " + str(eeg[:20] ) ) 

        #ph = lfilter(bAlpha,1,eeg)
        
        ##print ( " ph before angle = " + str(ph[:10] ) ) 
        
        #ph = np.append(ph[gdAlpha+1:], np.zeros(gdAlpha))
        #ph = np.angle(hilbert(ph))
        
        #print ( " ph after angle = " + str(ph[:10] ) ) 
        
        ##amplitude for gamma range + MI
        for key, value in filtersGamma.items(): 
            
            #amplitude
            amp1 = lfilter(value[0],1,eeg1)
            amp2 = lfilter(value[0],1,eeg2)
            gd = value[1]
            
            amp1 = np.append(amp1[gd+1:], np.zeros(gd))
            amp2 = np.append(amp2[gd+1:], np.zeros(gd))

            ph1 = np.angle(hilbert(amp1))
            ph2 = np.angle(hilbert(amp2))
            
            phd = ph1 - ph2
            #compute raw modulation index
            #print (" amp = " + str(amp[:20]) ) 
            #print (" ph = " + str(ph[:20]) ) 
            #print (" edges = " + str(edges) ) 
            
            PLVs.append((1/len(phd))*(abs(sum(np.exp(1j*phd))))) 

            #break
        
        print ( " MIs = " + str(PLVs))
        plt.plot(fGamma,PLVs)
        plt.show()

        #print ( [str(k) + ":" + str(v[1]) + "-" + str(v[0][:5]) for k,v in filtersGamma.items() ]  )
        
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return    

test()    
