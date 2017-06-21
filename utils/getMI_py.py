import numpy as np
import os
import sys
import traceback

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

    try:
        
        # get histogram for phases
        h = np.zeros(len(edges)-1)
        
        for hi in range ( len(edges) - 1 ):
            
            k = phase >= edges[hi] & phase < edges[hi+1]
            
            if sum(k) > 0: # only if there are some values, otherwise keep zero
                h[hi] = np.mean(amp[k]) # mean of amplitudes at that phase bin
    
        # fix last value
        k = [x == y for x,y in zip ( phase, edges [end])] 
        if not all([x == 0 for x in k]):
            h[end] = h[end] + np.mean(amp[k])
        
        # convert to probability
        h = h / sum(h)
        
        # replace zeros by eps
        k = h == 0
        h[k] = eps
        
        # calculate modulation index
        Hp = -1 * sum( [x*np.log(x) for x in h]) # entropy of the histogram
        D = np.log(len(h)) - Hp
        MI = D / np.log(len(h))
    
    
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return MI


