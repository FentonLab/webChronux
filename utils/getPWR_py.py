import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
#from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt

def getPWR():
    
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
	wFactor = 8
	bands = np.linspace( 2,60)
	wavelets = getWaveletsNorm( bands, wFactor, eegFS )
	
	eeg = eegData[:15]        

	eeg = [ x[50*eegFS:80*eegFS] for x in eeg] 	
	
	pwrs = np.zeros(1,len(wavelets))
	
	# filter
	for filtI in range (len(wavelets)):
	
	    psi = wavelets[filtI]
	
	    #convolution of wavelet and signal
	    
	    c = np.convolve (eeg,psi)
	    
	    #fix start and end
	    N = round((length(psi)-1)/2)
	    
	    c = c[N:len(c)-N]
	    
	    if len(c) > shape(eeg)[0]:
		c = c[shape(eeg)[0])
		
	    power = (abs(c))**2
	
	    pwrs[filtI] = np.mean(power)
	
	plt.plot(bands,pwrs)
	
    except:
        traceback.print_exc(file=sys.stdout)    
        
    return   	