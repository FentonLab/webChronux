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

def testkaiser():

    '''
    '''

    MI = ''

    try:

        lcut = 9.0
        hcut = 13.0
        sample_rate = 250.0
        nyq = sample_rate / 2.0
        width = 2.0/nyq # pass to stop transition width
        
        # The desired attenuation in the stop band, in dB.
        ripple_db = 40.0
        
        # Compute the order and Kaiser parameter for the FIR filter.
        N, beta = kaiserord(ripple_db, width)

        hpftaps = firwin(N, lcut/nyq, window=('kaiser', beta))    
        hpftaps = [-1*a for a in hpftaps]
        
        midPoint = int(np.round(len(hpftaps)/2))
        if midPoint % 2 != 0:
            midPoint = midPoint -1
        
        hpftaps[midPoint] = hpftaps[midPoint] + 1
        
        lpftaps = firwin(N, hcut/nyq, window=('kaiser', beta))
        lpftaps[midPoint] = lpftaps[midPoint] - 1
        
        taps = [sum(pair) for pair in zip(hpftaps, lpftaps)]
        
        plt.plot(taps)
        plt.show()
        
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return    

testkaiser()    
        