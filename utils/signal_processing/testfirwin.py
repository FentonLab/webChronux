#get artifacts, store results in file
import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from scipy.signal import firwin,firwin2,lfilter, kaiser
import matplotlib.pyplot as plt
import pandas as pd

def testfirwin():
    
    '''
    Input: 
    Output: 
    '''

    try:   
        
        plt.figure(1)
        

        
        n = 61
        a = signal.firwin(n, cutoff = 0.3, window = "hamming")
        plt.subplot(311)
        plt.plot(a)
        
        #Frequency and phase response
        w, h = mfreqz(a)

        plt.subplot(312)

        plt.plot(w, abs(h))       

        impz(a)
        plt.subplot(313)

        plt.show()
    except:
        traceback.print_exc(file=sys.stdout)    
    return  

testfirwin()