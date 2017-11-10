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

def analyzeEDFData(filePath):

  try:

    print ( " in analyze data ")
    #data = loadmat("/Users/smitra/projects/webChronux/utils/EEG181.mat")["eegData"]
    filePath = "/Users/smitra/self/andre/MIT-concussion/257802-post_season_20161206_142914.edf"
    f = pyedflib.EdfReader(filePath)
    signal_labels = f.getSignalLabels()
    numChannels = f.signals_in_file
  
    beginWin = 0
    endWin = 0
    
    samplingFrequency = eegFS = 250
    upperFrequency = 100
    lowerFrequency = 0
    timeBandWidth = 4
    timeWindow = 5 # time window in seconds
    STEP_SIZE = 2 # in seconds
    
    numDataPoints =  timeWindow * samplingFrequency
    stepSize = STEP_SIZE * samplingFrequency
    padding = pad = 0
  
    winLen = timeWindow * eegFS  
    paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + pad ) )
  
    numTapers = 2 * timeBandWidth -1
  
    [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )
  
    #numTapers = len(tapers)
    numTapers = 9
    print ("numTapers="+ str(numTapers))

    fpass = [lowerFrequency,upperFrequency]
  
    gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, eegFS)
    
    dataMatrix = []
  
    #spectrumChannelSumData = [0] * ( upperFrequency - lowerFrequency + 1 )
    spectrumChannelSumData = [0] * len(gridIndices) 
    
    for channelIndex in range(numChannels):
  
      spectrogramData = []
  
      #channelData = data[channelIndex]
  
      channelData = f.readSignal(channelIndex)
      
      channelLabel = signal_labels[channelIndex]
      
      # only process selected channels
      if channelLabel not in SELECTED_CHANNELS:
        continue
      #print ("for channel " + channelLabel)
  
      numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
  
      print ( " numWindows = " + str(numWindows) )
      
      #theta (4-7 Hz), alpha (9-13 Hz), beta (15-25 Hz) and gamma (30-50 Hz)       
      thetaPowerSpectra = []
      alphaPowerSpectra = []
      betaPowerSpectra = []
      gammaPowerSpectra = []
  
      for windowNum in range ( numWindows ) :
  
          beginWin = windowNum * samplingFrequency * STEP_SIZE
          endWin = beginWin + numDataPoints
          
          print (" window num = " + str(windowNum) )
  
          #print ( " beginWin = " + str(beginWin) + " endWin = " + str(endWin))
  
          windowData = channelData [ beginWin : endWin]
  
          #print ( " windowData = " + str(windowData) )
  
          if len(windowData) == 0:
  
            break

  
          for taperIndex, taper in enumerate ( tapers ) :
  
            #print (" taperIndex = " + str(taperIndex))
            #print (" tapers = " + str(tapers))
  
            taperData = [float(a)*float(b) for a,b in zip(windowData,taper)]
            
            fftData = scipy.fftpack.fft(taperData,paddedNumDataPoints)
            
            #print ( " fftData before = " + str(fftData) ) 
            
            fftData = [fftData[x] for x in gridIndices]

            fftData = (1.0/float(eegFS) ) * np.array (fftData)
             
            spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
  
            #spectrumChannelData = [ x for x in spectrumChannelData if x > lowerFrequency or x < upperFrequency ]
            #spectrumChannelData = spectrumChannelData[lowerFrequency:upperFrequency]

            #theta (4-7 Hz), alpha (9-13 Hz), beta (15-25 Hz) and gamma (30-50 Hz)       
            #thetaPowerSpectra = spectrumChannelData[4:7]
            #alphaPowerSpectra = []
            #betaPowerSpectra = []
            #gammaPowerSpectra = []

            #thetaPowerSpectra = spectrumChannelData[lowerFrequency:upperFrequency]
            #thetaPowerSpectra = spectrumChannelData[lowerFrequency:upperFrequency]
            #thetaPowerSpectra = spectrumChannelData[lowerFrequency:upperFrequency]

            spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)
            
            #fftData2 = [fftData2[x] for x in gridIndices]
            #fftData2 = (1.0/float(eegFS) ) * np.array (fftData2)
            #crossSpectrumChannelData = np.array ([np.conj(x) * y for x,y in zip ( fftData1, fftData2) ])            
            #crossSpectrumChannelDataList.append(crossSpectrumChannelData)
        
        #crossSpectrumChannelDataAvg = np.mean(crossSpectrumChannelDataList, axis = 0)                
  
          spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
          
          #spectrumChannelAvgData = np.mean( spectrumChannelSumData, axis = 0 ) 
  
          spectrogramData.append(list(spectrumChannelAvgData))
          
          plt.plot(spectrumChannelAvgData)
          plt.show()
  
          #print (spectrogramData)
  
          #print (" for window = " + str(windowNum) + " spectrogramData = " + str(spectrogramData) )
  
          break
  
      #np.savetxt("outdata/channel_spectrogram_data" + str(channelIndex) + ".txt", spectrogramData )
      #print (spectrogramData)
      #print(np.shape(spectrogramData))
      
      #spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
      #spectrumPSD = np.array(spectrumPSD)/100
      #np.savetxt("outdata/channel_spectrum_PSD_data" + str(channelIndex) + ".txt", spectrumPSD )
  
      #plt.plot(spectrumPSD.transpose())
      #plt.savefig("outdata/channel_spectrum_psd_" + str(channelIndex) + ".png" )
      #plt.clf()
      
      #plt.figure(1, figsize = (8.5,11))
      #plt.imshow(np.array(spectrogramData))
      #plt.savefig("outdata/channel_spectrogram_" + str(channelIndex) + ".png" )
      #plt.gca().invert_yaxis()
      #plt.show()
      #plt.clf()

      break

  except:
    traceback.print_exc(file=sys.stdout)
    return

analyzeEDFData("/Users/smitra/self/andre/MIT-concussion/257802-post_season_20161206_142914.edf")
#257807-_post-season_20161129_132254.edf
print (" end 999 ")
