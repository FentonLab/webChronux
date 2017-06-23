import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
import matplotlib.pyplot as plt
from math import ceil
from multitaper import *
import scipy

# get frequency grid
def getGridIndices(lowerFrequency, upperFrequency, paddedNumDataPoints, samplingFrequency):

  try:

      frequencyResolution = float ( samplingFrequency ) / float ( paddedNumDataPoints )
      
      gridValues = np.arange ( 0, samplingFrequency , frequencyResolution )
      
      gridValues = gridValues[ :paddedNumDataPoints ]

      gridIndices = [index for index, x in enumerate (gridValues) if x>= lowerFrequency and x<= upperFrequency ]

      gridValues = [x for index, x in enumerate (gridValues) if x>= lowerFrequency and x<= upperFrequency ]

  except:
    traceback.print_exc(file=sys.stdout)

  return gridValues , gridIndices

def getWPLI():
    
    '''
    compute phase lag variance
    PLV plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    
    '''
    
    WINDOW_SIZE = 3 # window size in seconds
    
    try:    

        eegData = loadmat('EEG181.mat')
        
        eegFS = 250 # sampling frequency
        
        eegData = eegData["eegData"]
        
        eeg1 = eegData[14]
        eeg2 = eegData[15]        

        eeg1 = eeg1[10*eegFS - 1 : 100*eegFS - 1 ]  
        eeg2 = eeg2[10*eegFS - 1 : 100*eegFS - 1 ]  

        eeg1 = eeg1 - np.mean(eeg1, axis = 0)        
        eeg2 = eeg2 - np.mean(eeg2, axis = 0)        
        
        winLen = WINDOW_SIZE * eegFS
        
        tapers = [4,9]
        pad = 0
        #Fs = 250
        fpass = [0,125]
        timeBandWidth = 4
        numTapers = 9
        
        paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + pad ) )

        gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, eegFS)

        NF = ( paddedNumDataPoints / 2 ) + 1
 
        timeBandWidth = 4
        numTapers = 9
 
        [tapers, eigenValues] = dpss_windows(int(winLen),  timeBandWidth, numTapers)

        tapers = np.array(tapers)
        tapers = math.sqrt(eegFS) * tapers
        
        Nw = int(np.floor(len(eeg1)/winLen))

        S12x = []
        S1x = []
        S2x = []
        
        crossSpectrumWinDataList = []
        
        for wI in range ( Nw ) :
        
            start = (wI) * winLen
            end = start + winLen
        
            #print ( " start = " + str(start) + " end = " + str(end) ) 
        
            eeg1x = eeg1[start:end]
            eeg2x = eeg2[start:end]
            
            crossSpectrumChannelDataList = []

            for taperIndex, taper in enumerate ( tapers ) :
     
                taperData1 = [float(a)*float(b) for a,b in zip(eeg1x,taper)]
                
                fftData1 = scipy.fftpack.fft(taperData1,paddedNumDataPoints)
                
                fftData1 = [fftData1[x] for x in gridIndices]
                
                fftData1 = (1.0/float(eegFS) ) * np.array (fftData1)
     
                taperData2 = [float(a)*float(b) for a,b in zip(eeg2x,taper)]
                fftData2 = scipy.fftpack.fft(taperData2,paddedNumDataPoints)

                fftData2 = [fftData2[x] for x in gridIndices]

                fftData2 = (1.0/float(eegFS) ) * np.array (fftData2)
                
                crossSpectrumChannelData = np.array ([np.conj(x) * y for x,y in zip ( fftData1, fftData2) ])            

                crossSpectrumChannelDataList.append(crossSpectrumChannelData)
            
            crossSpectrumChannelDataAvg = np.mean(crossSpectrumChannelDataList, axis = 0)     
            
            crossSpectrumChannelDataAvg = [x.imag for x in crossSpectrumChannelDataAvg]
            
            crossSpectrumWinDataList.append(crossSpectrumChannelDataAvg)

        outsum = np.nansum(crossSpectrumWinDataList, axis = 0)
        outsumW = np.nansum(np.absolute(crossSpectrumWinDataList), axis = 0)
        outssq = np.nansum(np.square(crossSpectrumWinDataList), axis = 0)
        
        wpli     = (np.square(outsum) - outssq) / (np.square (outsumW) - outssq)

        plt.plot(gridValues,wpli)            
        
        plt.show()
            
    except:
      traceback.print_exc(file=sys.stdout)
    return
        
getWPLI()
        

