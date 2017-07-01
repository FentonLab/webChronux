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

        #length of t has to be even
        if len(t) % 2 == 0 :

            t = t[:-1]

        S1 = np.exp(    (-1*(t ** 2)) / (2 * ( sigmaT ** 2 ) )    ) # gaussian curve
        S2 = np.exp(2*1j*pi*f*t) # sinewave
        
        A = (sigmaT * math.sqrt(pi))  ** -0.5  #normalization for total power = 1
        psi = A * np.array([x*y for x,y in zip ( S1,S2) ])

        s = sum([abs(x.real) for x in psi])/2 #divide by total energy so convolution results in one

        psi = psi / s

        wavelets[band] = psi
        
    return wavelets

def getPWR():
    try:
        eegData = loadmat('EEG181.mat')
        
        eegData = eegData["eegData"]
        
        eegFS = 250
        wFactor = 8
        bands = np.arange( 2,61)
        #print (bands)
        wavelets = getWaveletsNorm( bands, wFactor, eegFS )
        #print ( (eegData) ) 
        eeg = eegData[14]        
        eeg = eeg[50*eegFS - 1 : 80*eegFS - 1 ]   
        
        #print ( eeg[:10])
        
        pwrs = [] 

        # filter
        for filtI, psi in wavelets.items():

            #convolution of wavelet and signal
            
            #print (" psi = " + str( psi[:10] ) )            

            c = np.convolve (eeg,psi)
            
            #print (len(c))

            #print (c[:10])
            
            #break

            #fix start and end
            N = round((len(psi)-1)/2)

            print ( "n = " + str(N))
            c = c[N-1:len(c)-N-1]
            print ( " c = " + str(c[:10]))

            if len(c) > len(eeg):
                c = c[:len(eeg)]

            power = (abs(c))**2

            #print ( " power = " + str(power[:10]) )
            #pwrs[filtI] = np.mean(power)
            
            #print ( " mean power = " + str(np.mean ( power ) ) ) 
                    
            pwrs.append(np.mean(power))

            #break
        plt.plot(bands,pwrs)
        plt.show()
                
    except:
        traceback.print_exc(file=sys.stdout)
    return

getPWR()