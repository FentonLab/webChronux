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

samplingFrequency = 250
upperFrequency = 100
lowerFrequency = 1
timeBandWidth = 4
numTapers = 9
timeWindow = 3 # time window in seconds
STEP_SIZE = 2 # in seconds
numDataPoints =  timeWindow * samplingFrequency
stepSize = STEP_SIZE * samplingFrequency
padding = 1
paddedNumDataPoints = int ( pow (2, ceil ( math.log( numDataPoints, 2 )) + padding ))

def getGridIndices():

  try:

     if upperFrequency and lowerFrequency:
        paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( numDataPoints ) ) + padding ) )

        frequencyResolution = float ( samplingFrequency ) / float ( paddedNumDataPoints )
        gridValues = np.arange ( 0, paddedNumDataPoints , frequencyResolution )
        gridIndices = np.where ( (gridValues >= lowerFrequency ) & (gridValues <= upperFrequency ) )

        upperFrequencyGrid = 0
        lowerFrequencyGrid = 0
        # print " gridIndices " + str(gridIndices)

        if len( gridIndices ) > 0:

          upperFrequencyGrid = gridIndices[0] [len(gridIndices[0]) -1]

          if int(gridIndices [0][0]) > 0:

              lowerFrequencyGrid = int(gridIndices [0][0])

  except:
         traceback.print_exc(file=sys.stdout)
  print ( " grid indices = " + str(gridIndices))
  return int(lowerFrequencyGrid), int(upperFrequencyGrid) , gridIndices

def analyzeEDFData():

 try:

  print ( " in analyze data ")
  f = pyedflib.EdfReader("data/test.edf")
  n = f.signals_in_file

  signal_labels = f.getSignalLabels()

  #print ( " num signals = " + str(n) )
  #print ( " signal labels = " + str(signal_labels) )

  beginWin = 0
  endWin = 0

  numTapers = 2 * timeBandWidth -1

  [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )

  spectrumChannelSumData = [0] * ( upperFrequency - lowerFrequency + 1 )

  numTapers = len(tapers)

  lowerFrequencyGrid, upperFrequencyGrid, gridIndices = getGridIndices ()

  for channelIndex in range(n):

     spectrogramData = []

     crossSpectrumData = []

     channelData = f.readSignal(channelIndex)

     #numWindows = round ( len(channelData) / numDataPoints ) + 1

     numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )

     print ( " numWindows = " + str(numWindows) )

     for windowNum in range ( numWindows ) :

        beginWin = windowNum * numDataPoints * STEP_SIZE
        endWin = beginWin + numDataPoints

        windowData = channelData [ beginWin : endWin]

        if len(windowData) == 0:

           break

        for channelIndex2 in range ( n ) :

            if channelIndex != channelIndex2:

              channelData2 = f.readSignal(channelIndex2)

              windowData2 = channelData2 [ beginWin : endWin]

        for taperIndex, taper in enumerate ( tapers ) :

           taperData = [float(a*b) for a,b in zip(windowData,taper)]
           taperData2 = [float(a*b) for a,b in zip(windowData2,taper)]

           fftData = fft(taperData,paddedNumDataPoints)

           fftData2 = fft(taperData,paddedNumDataPoints)

           spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
           spectrumChannelData2 = np.array([log(abs(x*conj(x))) for x in fftData2])

           crossSpectrumData = np.array([log(abs(x*conj(y))) for x,y in zip(fftData, fftData2)])

           print ( " padded num = " + str(paddedNumDataPoints) + " spectrum len = " + str(len(spectrumChannelData)))

           spectrumChannelData = list(spectrumChannelData[gridIndices])

           spectrumChannelData = (1 / float(samplingFrequency) ) * np.array(spectrumChannelData)
           spectrumChannelData = spectrumChannelData[lowerFrequency:upperFrequency+1]

           #print (spectrumChannelData)

           #print (" for taper = " + str(taperIndex ) + "  spectrumChannelData = " + str(spectrumChannelData) )
           spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)

        spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers

        spectrogramData.append(spectrumChannelAvgData)

        #print (spectrogramData)

        #print (" for window = " + str(windowNum) + " spectrogramData = " + str(spectrogramData) )

        #break

     np.savetxt("outdata/channel_spectrogram_data" + str(channelIndex) + ".txt", spectrogramData )
     #print (spectrogramData)

     np.savetxt("outdata/channel_spectrogram_data" + str(channelIndex) + ".txt", spectrogramData )
     spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
     spectrumPSD = np.array(spectrumPSD)/100
     np.savetxt("outdata/channel_spectrum_PSD_data" + str(channelIndex) + ".txt", spectrumPSD )

     plt.plot(spectrumPSD.transpose())
     plt.savefig("outdata/channel_spectrum_psd_" + str(channelIndex) + ".png" )

     plt.clf()

     plt.imshow(np.array(spectrogramData))
     plt.savefig("outdata/channel_spectrogram_" + str(channelIndex) + ".png" )

     plt.clf()

     break

 except:
     traceback.print_exc(file=sys.stdout)
 return

analyzeEDFData()
print (" end 999 ")
