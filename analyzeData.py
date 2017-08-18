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

def analyzeEDFData(filePath):

  try:

    print ( " in analyze data ")
    f = pyedflib.EdfReader(filePath)
    n = f.signals_in_file
  
    signal_labels = f.getSignalLabels()
  
    beginWin = 0
    endWin = 0
    
    samplingFrequency = eegFS = 250
    upperFrequency = 100
    lowerFrequency = 0
    timeBandWidth = 4
    timeWindow = WINDOW_SIZE = 5 # time window in seconds
    STEP_SIZE = 2 # in seconds
    numDataPoints =  timeWindow * samplingFrequency
    stepSize = STEP_SIZE * samplingFrequency
    padding = pad = 1
  
    winLen = WINDOW_SIZE * eegFS  
    paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + pad ) )
  
    numTapers = 2 * timeBandWidth -1
  
    [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )
  
    numTapers = len(tapers)
    print ("numTapers="+ str(numTapers))
    pad = 0
    fpass = [0,100]
    timeBandWidth = 4
    numTapers = 9
  
    gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, eegFS)
  
    #spectrumChannelSumData = [0] * ( upperFrequency - lowerFrequency + 1 )
    spectrumChannelSumData = [0] * len(gridIndices)
    
    for channelIndex in range(n):
  
      spectrogramData = []
  
      channelData = f.readSignal(channelIndex)
  
      numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
  
      print ( " numWindows = " + str(numWindows) )
  
      for windowNum in range ( numWindows ) :
  
          beginWin = windowNum * samplingFrequency * STEP_SIZE
          endWin = beginWin + numDataPoints
  
          #print ( " beginWin = " + str(beginWin) + " endWin = " + str(endWin))
  
          windowData = channelData [ beginWin : endWin]
  
          #print ( " windowData = " + str(windowData) )
  
          if len(windowData) == 0:
  
            break
  
          for taperIndex, taper in enumerate ( tapers ) :
  
            print (" taperIndex = " + str(taperIndex))
            print (" tapers = " + str(tapers))
  
            taperData = [float(a*b) for a,b in zip(windowData,taper)]
            
            fftData = fft(taperData,paddedNumDataPoints)
            
            print ( " fftData before = " + str(fftData) ) 
            
            fftData = [fftData[x] for x in gridIndices]
            
            #print ( " fftData after ^^^^^^^^^ = " + str(fftData) ) 
  
            spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
  
            spectrumChannelData = spectrumChannelData[lowerFrequency:upperFrequency]
  
            #print ( " padded num = " + str(paddedNumDataPoints) + " spectrum len = " + str(spectrumChannelData))
  
            #spectrumChannelData = list(spectrumChannelData[gridIndices])
  
            spectrumChannelData = (1 / float(samplingFrequency) ) * np.array(spectrumChannelData)
            #spectrumChannelData = spectrumChannelData[lowerFrequency:upperFrequency+1]
  
            #print (spectrumChannelData)
  
            #print (" for taper = " + str(taperIndex ) + "  spectrumChannelData = " + str(spectrumChannelData) )
            spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)
  
          spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
  
          spectrogramData.append(list(spectrumChannelAvgData))
  
          #print (spectrogramData)
  
          #print (" for window = " + str(windowNum) + " spectrogramData = " + str(spectrogramData) )
  
          #break
  
      #np.savetxt("outdata/channel_spectrogram_data" + str(channelIndex) + ".txt", spectrogramData )
      #print (spectrogramData)
      print(np.shape(spectrogramData))
      #np.savetxt("outdata/channel_spectrogram_data" + str(channelIndex) + ".txt", spectrogramData )
      #spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
      #spectrumPSD = np.array(spectrumPSD)/100
      #np.savetxt("outdata/channel_spectrum_PSD_data" + str(channelIndex) + ".txt", spectrumPSD )
  
      #plt.plot(spectrumPSD.transpose())
      #plt.savefig("outdata/channel_spectrum_psd_" + str(channelIndex) + ".png" )
  
      #plt.clf()
  
      plt.imshow(spectrogramData)
      plt.savefig("outdata/channel_spectrogram_" + str(channelIndex) + ".png" )
  
      plt.show()
  
      #plt.clf()
  
      #S = imag(S12x);
  
      #WPLI
       #outsum   = nansum(S,1) # compute the sum; this is 1 x size(2:end)
       #outsumW  = nansum(abs(S),1) # normalization of the WPLI
       #outssq   = nansum(S.^2,1)
       #wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq) # do the pairwise thing in a handy way


      break

  except:
    traceback.print_exc(file=sys.stdout)
    return

analyzeEDFData("/Users/smitra/self/andre/MIT-concussion/baseline257802_20161018_194838.edf")
print (" end 999 ")
