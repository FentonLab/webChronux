import pyedflib
import scipy
from scipy import *
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import traceback
from scipy.io import savemat, loadmat
from multitaper import *

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

def showFig(filePath):

  try:
  
    
    data = loadmat(filePath)["S1x"]
    
    print(np.shape(data))
    print (data[:20,:20])
    #plt.figure(1, figsize = (8.5,11))
    #plt.imshow(np.log10(np.array(data)).transpose())
    #plt.gca().invert_yaxis()    
    #plt.show()

    data = loadmat("/Users/smitra/projects/webChronux/utils/EEG181.mat")["eegData"]

    beginWin = 0
    endWin = 0
    samplingFrequency = 250
    upperFrequency = 100
    lowerFrequency = 0
    timeBandWidth = 4
    timeWindow = 3 # time window in seconds
    stepSize = 1 # in seconds
    
    numDataPoints =  timeWindow * samplingFrequency
    stepSizeWindow = stepSize * samplingFrequency
    padding = 0
  
    winLen = timeWindow * samplingFrequency  
    print ("win len = " + str(winLen))
    paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + padding ) )
  
    numTapers = 9
    [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )

    fpass = [lowerFrequency,upperFrequency]
  
    gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, samplingFrequency)
    
    dataMatrix = []
  
    spectrumChannelSumData = [0] * len(gridIndices) 

    spectrogramData = []

    channelData = data[14]

    #numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
    numWindows = int ( ( len ( channelData )) / ( numDataPoints  ) )

    print ( " numWindows = " + str(numWindows) )
    
    #theta (4-7 Hz), alpha (9-13 Hz), beta (15-25 Hz) and gamma (30-50 Hz)       

    for windowNum in range ( numWindows ) :

        beginWin = windowNum * numDataPoints
        endWin = beginWin + numDataPoints
        
        windowData = channelData [ beginWin : endWin]
        if len(windowData) == 0:
          break

        for taperIndex, taper in enumerate ( tapers ) :
          print ("window data = " + str(np.shape(windowData)))
          print ("taper = " + str(np.shape(taper)))
          taperData = [float(a)*float(b) for a,b in zip(windowData,taper)]
          fftData = scipy.fftpack.fft(taperData,paddedNumDataPoints)
          fftData = [fftData[x] for x in gridIndices]
          fftData = (1.0/float(samplingFrequency) ) * np.array (fftData)
          print (['({0.real:.4f} + {0.imag:.4f}i)'.format(x) for x in fftData[:20]])
          spectrumChannelData = np.array([abs(x*conj(x)) for x in fftData])
          #break
          spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)

        spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
        #plt.plot([log(x) for x in spectrumChannelAvgData])
        #plt.show()
        #break

        spectrogramData.append(list(spectrumChannelAvgData))
    print(np.shape(spectrogramData))
    spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
    spectrumPSD = np.array(spectrumPSD)/100
    #plt.plot([log(x) for x in spectrumPSD])
    #plt.show()    
    plt.figure(1, figsize = (8.5,11))
    plt.imshow(np.array(log(spectrogramData)).transpose())
    #plt.gca().invert_yaxis()
    plt.show()

  except:
    traceback.print_exc(file=sys.stdout)
    return

showFig("/Users/smitra/projects/webChronux/matspectrogram.mat")
print (" end 999 ")
