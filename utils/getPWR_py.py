import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from math import pi
#from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt
#import pdb
def getWaveletsNorm( bands, wFactor, eegFS ):

    #bands ["alpha"] = np.linspace(4,7)
    #bands ["beta"]  = np.linspace(9,12)
    #bands ["low_gamma"] = np.linspace(15,50,5)
    #bands ["high_gamma"] = np.linspace(70,90,5)

    wavelets = {} # wavelets in each band    
    
    for  index, band in enumerate ( bands ) :
        
        f = band
        sigmaF = f/wFactor    #practical setting for spectral bandwidth (Tallon-Baudry)
        sigmaT = 1/(sigmaF*pi*2) #wavelet duration
        t = np.arange ( -4*sigmaT*eegFS , 4*sigmaT*eegFS , 1)

        t = t/eegFS
        
        #print (" before " + str( t[-3:]) )
        #print (len(t) )

        #length of t has to be even
        if len(t) % 2 == 0 :
            #print ( " in mod ") 
            t = t[:-1]
        #print (" before " + str( t[-3:]) )

        S1 = np.exp(    (-1*(t ** 2)) / (2 * ( sigmaT ** 2 ) )    ) # gaussian curve
        S2 = np.exp(2*1j*pi*f*t) # sinewave
        
        #print ( " S1 = " + str(S1[-3:]))
        #print ( " S2 = " + str(S2[-3:]))
        
        A = (sigmaT * math.sqrt(pi))  ** -0.5  #normalization for total power = 1
        psi = A * np.array([abs((x*y).real) for x,y in zip ( S1,S2) ])
        #psi =  [x*y for x,y in zip ( S1,S2) ]

        #print (A)
        
        #print ( psi)
        
        #psi = real(psi); %to make output real
        s = sum(psi)/2 #divide by total energy so convolution results in one
        #print ( s)        
        psi = psi / s
        #print ( psi[:10])
        wavelets[band] = psi
        
        #print ( psi)
        #break

    return wavelets

def getPWR():
    try:
        eegData = loadmat('EEG181.mat')
        eegFS = 250
        wFactor = 8
        bands = np.arange( 2,61)
        #print (bands)
        wavelets = getWaveletsNorm( bands, wFactor, eegFS )
    except:
        traceback.print_exc(file=sys.stdout)
    return

getPWR()



#def getPWR():

    #'''
    #compute phase lag variance
    #PLV plotted at the end
    #Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    #Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    #Input:
    #Output:
    #'''

    #try:
        #eegData = loadmat('EEG181.mat')
        #eegFS = 250
        #wFactor = 8
        #bands = np.linspace( 2,60)
        ##wavelets = getWaveletsNorm( bands, wFactor, eegFS )

        ##eeg = eegData[:15]

        ##eeg = [ x[50*eegFS:80*eegFS] for x in eeg]

        ##pwrs = np.zeros(1,len(wavelets))

        ### filter
        ##for filtI in range (len(wavelets)):

            ##psi = wavelets[filtI]

            ###convolution of wavelet and signal

            ##c = np.convolve (eeg,psi)

            ###fix start and end
            ##N = round((length(psi)-1)/2)

            ##c = c[N:len(c)-N]

            ##if len(c) > shape(eeg)[0]:
                ##c = c[shape(eeg)[0])

            ##power = (abs(c))**2

            ##pwrs[filtI] = np.mean(power)

        ##plt.plot(bands,pwrs)

    #except:
        #traceback.print_exc(file=sys.stdout)

    #return