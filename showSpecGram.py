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
SELECTED_CHANNELS = ['F7','T5', 'F3','P3', 'F4','P4', 'F8','T6', 'Fp1','O1', 'Fp2','O2']

# get frequency grid
def getGridIndices(lowerFrequency, upperFrequency, paddedNumDataPoints, samplingFrequency):

  try:

      frequencyResolution = float ( samplingFrequency ) / float ( paddedNumDataPoints )
      
      print (" frequencyResolution = " + str(frequencyResolution))
      
      gridValues = np.arange ( 0, samplingFrequency , frequencyResolution )
      
      print(" len gridValues = " + str(len(gridValues)))
      
      gridValues = gridValues[ :paddedNumDataPoints ]
      
      print ( " len grid values 000 = " + str(len(gridValues )))

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

    #data = loadmat("/Users/smitra/projects/webChronux/utils/EEG181.mat")["eegData"]
    filePath = "/Users/smitra/self/andre/MIT-concussion/257802-post_season_20161206_142914.edf"
    f = pyedflib.EdfReader(filePath)
    signal_labels = f.getSignalLabels()
    numChannels = f.signals_in_file

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
    #print ("win len = " + str(winLen))
    paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + padding ) )
  
    numTapers = 9
    [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth)/2, int(numTapers) )

    fpass = [lowerFrequency,upperFrequency]
  
    gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, samplingFrequency)

    print(" len gridValues = " + str(len(gridValues)))
    print(" len gridIndices = " + str(len(gridIndices)))
    
    dataMatrix = []
    
    for channelIndex in range(numChannels):
  
      spectrogramData = []
  
      #channelData = data[channelIndex]
  
      channelData = f.readSignal(channelIndex)
      
      channelLabel = signal_labels[channelIndex]
      
      # only process selected channels
      #if channelLabel not in SELECTED_CHANNELS:
        #continue
        
      dataMatrix.append(channelData)
      if channelIndex != 14:
        continue

      #print ("for channel " + channelLabel)
  
      numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
  
      print ( " numWindows = " + str(numWindows) )

      spectrumChannelSumData = [0] * len(gridIndices) 
  
      spectrogramData = []
  
      #channelData = data[14]
  
      #numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
      numWindows = int ( ( len ( channelData )) / ( numDataPoints  ) )
  
      #numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
  
      #print ( " numWindows = " + str(numWindows) )
      
      #theta (4-7 Hz), alpha (9-13 Hz), beta (15-25 Hz) and gamma (30-50 Hz)       
  
      duration = numWindows * timeWindow
      
      for windowNum in range ( numWindows ) :
    
          print ( " for window = " + str(windowNum) )
  
          beginWin = windowNum * numDataPoints
          endWin = beginWin + numDataPoints
          
          windowData = channelData [ beginWin : endWin]
          if len(windowData) == 0:
            break
          
          #if windowNum > 200:
            #break
            
          spectrumChannelDataList = []
  
  
          print ( " len gridIndices = " + str(len(gridIndices)) )
          for taperIndex, taper in enumerate ( tapers ) :
            #print ("window data = " + str(np.shape(windowData)))
            #print ("taper = " + str(np.shape(taper)))
            taperData = [float(a)*float(b) for a,b in zip(windowData,taper)]
            fftData = scipy.fftpack.fft(taperData,paddedNumDataPoints)
            fftData = [fftData[x] for x in gridIndices]
            fftData = (1.0/float(samplingFrequency) ) * np.array (fftData)
            #print (['({0.real:.4f} + {0.imag:.4f}i)'.format(x) for x in fftData[:20]])
            spectrumChannelData = np.array([abs(x*conj(x)) for x in fftData])
            #break
            spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)
            
            spectrumChannelDataList.append(spectrumChannelData)
  
          spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
          #plt.plot([log10(x) for x in spectrumChannelDataList[5]])
          #plt.show()
          #break
          print (" appending " + str(len(list(spectrumChannelAvgData))))
          spectrogramData.append(list(spectrumChannelAvgData))
          print (" size " + str(np.shape(spectrogramData)))

      #print(np.shape(spectrogramData))
      #spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
      #spectrumPSD = np.array(spectrumPSD)/100

      #plt.plot([log(x) for x in spectrumPSD])
      #plt.show()    
      
      fig = plt.figure(1, figsize = (20,30))
      
      #fig = plt.figure()
      
      axes = fig.add_subplot(211)
      axes.set_title("Spectrogram")
      axes.set_autoscalex_on(False)

      #print (" size spectrogramData = " + str(np.shape(np.array(log(spectrogramData)).transpose())))
      
      #x = np.array(log(spectrogramData)).transpose()
      
      #print(x[:20][:20])
      plt.imshow(np.array(log(spectrogramData)).transpose())
      plt.gca().invert_yaxis()
      
      plt.xlim((0,numWindows))      
      plt.ylim((lowerFrequency,upperFrequency))
      plt.ylabel('Frequency')
      plt.xlabel('Time')
      
      #plt.show()
      break
    savemat("/Users/smitra/projects/webChronux/utils/257802-post_season_20161206_142914.mat", {"eegData":dataMatrix})
    savemat("/Users/smitra/projects/webChronux/utils/spectrogram.mat", {"spectrogramData":spectrogramData})

  except:
    traceback.print_exc(file=sys.stdout)
    return

showFig("/Users/smitra/projects/webChronux/matspectrogram.mat")
print (" end 999 ")
