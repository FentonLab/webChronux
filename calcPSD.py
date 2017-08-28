import pyedflib
from scipy import *
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import traceback
import mne
from multitaper import *
import scipy
from scipy.io import savemat, loadmat
SELECTED_CHANNELS = ['F7','T5', 'F3','P3', 'F4','P4', 'F8','T6', 'Fp1','O1', 'Fp2','O2']

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

def calcAvgPSD(filePath):

  try:

    print ( " in analyze data ")
    f = pyedflib.EdfReader(filePath)
    n = f.signals_in_file
  
    signal_labels = f.getSignalLabels()
  
    beginWin = 0
    endWin = 0
    
    samplingFrequency = 250
    upperFrequency = 100
    lowerFrequency = 0
    
    timeBandWidth = 4
    timeWindow = 5 # time window in seconds
    stepSize = 2 # in seconds
    
    numDataPoints =  timeWindow * samplingFrequency
    padding = 1
  
    winLen = timeWindow * samplingFrequency  
    paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + padding ) )
  
    numTapers = 2 * timeBandWidth -1
  
    [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )
  
    #numTapers = len(tapers)
    numTapers = 9
    print ("numTapers="+ str(numTapers))

    fpass = [lowerFrequency,upperFrequency]
  
    gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, samplingFrequency)
    
    dataMatrix = []
  
    spectrumChannelSumData = [0] * ( upperFrequency - lowerFrequency )
    #spectrumChannelSumData = [0] * len(gridIndices) 
    
    for channelIndex in range(n):
  
      spectrogramData = []
  
      channelData = f.readSignal(channelIndex)
      
      channelLabel = signal_labels[channelIndex]
      
      # only process selected channels
      if channelLabel not in SELECTED_CHANNELS:
        continue
      print ("for channel " + channelLabel)
  
      numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize * samplingFrequency  ) )
  
      print ( " numWindows = " + str(numWindows) )
      
      #theta (4-7 Hz), alpha (9-13 Hz), beta (15-25 Hz) and gamma (30-50 Hz)       
      thetaPowerSpectra = []
      alphaPowerSpectra = []
      betaPowerSpectra = []
      gammaPowerSpectra = []
  
      for windowNum in range ( numWindows ) :
  
          beginWin = windowNum * samplingFrequency * stepSize
          endWin = beginWin + numDataPoints
         
          print (" window num = " + str(windowNum) )

          windowData = channelData [ beginWin : endWin]

          if len(windowData) == 0:
            break
  
          for taperIndex, taper in enumerate ( tapers ) :
  
            #print (" taperIndex = " + str(taperIndex))
            #print (" tapers = " + str(tapers))
  
            taperData = [float(a)*float(b) for a,b in zip(windowData,taper)]
            
            fftData = scipy.fftpack.fft(taperData,paddedNumDataPoints)
            #print ( " fftData before = " + str(fftData) ) 

            fftData = (1.0/float(samplingFrequency) ) * np.array (fftData)
            fftData = [fftData[x] for x in gridIndices]
             
            #spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
            spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
            spectrumChannelData = spectrumChannelData[lowerFrequency:upperFrequency]
            spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)
  
          spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
  
          spectrogramData.append(list(spectrumChannelAvgData))
  
          if windowNum > 10:
  
            break  
  
          #break
  
      print(np.shape(spectrogramData))
      
      spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
      spectrumPSD = np.array(spectrumPSD)/100
      plt.plot(spectrumPSD.transpose())
      
      plt.show()

      break

  except:
    traceback.print_exc(file=sys.stdout)
    return

calcAvgPSD("/Users/smitra/self/andre/MIT-concussion/257802-post_season_20161206_142914.edf")
print (" end 999 ")
