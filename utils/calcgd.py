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

def calcgd(gd):
        
    '''
    compute phase amplitude coupling, as in Vinnick
    PAC plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    '''

    try:
        
        a = np.zeros(len(gd))
        
        a[0] = 1
        a[-1] = 1

        c = np.convolve(gd, a)

        c = c.transpose()
        nc = len(c)
        #cr = c*(0:(nc-1));
         
        if isWholeUnitCircle:
            s=1
        else:
            s=2
         
        if len(n) == 1:
            w = (2 * math.pi / s*(0:n-1)/n)
            if s*n >= nc: # pad with zeros to get the n values needed
               # dividenowarn temporarily suppresses warnings to avoid "Divide by zero"
                try:
                    gd = (fft([cr zeros(1,s*n-nc)]),fft([c zeros(1,s*n-nc)]))
                except: 
                    pass
                gd = np.real(gd[1:n]) - np.ones[n]*(na-1)
            else: # find multiple of s*n points greater than nc
                nfact = s*ceil(nc/(s*n))
                mmax = n*nfact
                # dividenowarn temporarily suppresses warnings to avoid "Divide by zero"
                gd = dividenowarn(fft(cr,mmax), fft(c,mmax))
                gd = real(gd(1:nfact:mmax)) - ones(1,n)*(na-1)
            
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return    

calcgd()       
