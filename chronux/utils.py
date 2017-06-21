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
from django.conf import settings
import pyedflib

# get frequency grid
def getGridIndices(analysisObj):

  try:

    if analysisObj.upperFrequency and analysisObj.lowerFrequency:

        paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( analysisObj.numDataPoints ) ) + analysisObj.padding ) )

        #print ( " paddedNumDataPoints = " + str(paddedNumDataPoints))

        frequencyResolution = float ( analysisObj.samplingFrequency ) / float ( paddedNumDataPoints )
        
        gridValues = np.arange ( 0, paddedNumDataPoints , frequencyResolution )

        #print ( " frequencyResolution = " + str(frequencyResolution))
        
        #print ( " gridValues = " + str(gridValues))
        
        gridIndices = np.where ( (gridValues >= analysisObj.lowerFrequency ) & (gridValues <= analysisObj.upperFrequency ) )

        #print (" gridIndices " + str(gridIndices))

        upperFrequencyGrid = 0
        lowerFrequencyGrid = 0

        if len( gridIndices ) > 0:

          upperFrequencyGrid = gridIndices[0] [len(gridIndices[0]) -1]

          if int(gridIndices [0][0]) > 0:

              lowerFrequencyGrid = int(gridIndices [0][0])

  except:
    traceback.print_exc(file=sys.stdout)
  #print ( " grid indices = " + str(gridIndices))
  return int(lowerFrequencyGrid), int(upperFrequencyGrid) , gridIndices

# 
def analyzeEDFData(analysisObj):
  
  try:

    for datafile in analysisObj.datafiles:
  
      print ( " in analyze data ")
      f = pyedflib.EdfReader(datafile.filePath)
      n = f.signals_in_file
    
      signal_labels = f.getSignalLabels()
    
      print ( " num signals = " + str(n) )
      print ( " signal labels = " + str(signal_labels) )
    
      beginWin = 0
      endWin = 0
    
      numTapers = 2 * analysisObj.bandWidth -1
      
      #print (" numDataPoints =" + str(int(analysisObj.numDataPoints)) + " bandWidth = " + str(float(analysisObj.bandWidth)) + " numtapers= " + str(int(numTapers)))
    
      [tapers, eigenValues] = dpss_windows(int(analysisObj.numDataPoints), float(analysisObj.bandWidth), int(numTapers) )
    
      spectrumChannelSumData = [0] * (analysisObj.upperFrequency - analysisObj.lowerFrequency + 1 )
    
      numTapers = len(tapers)
    
      lowerFrequencyGrid, upperFrequencyGrid, gridIndices = getGridIndices (analysisObj)
      
      paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( analysisObj.numDataPoints ) ) + analysisObj.padding ) )
    
      for channelIndex in range(n):
        
        if signal_labels[channelIndex] not in analysisObj.spectrogramChannels:
          continue
    
        spectrogramData = []
   
        channelData = f.readSignal(channelIndex)
   
        #numWindows = round ( len(channelData) / numDataPoints ) + 1
   
        numWindows = int ( ( len ( channelData ) - analysisObj.numDataPoints + 1) / ( analysisObj.stepSize  ) )
   
        #print ( " numWindows = " + str(numWindows) )
   
        for windowNum in range ( numWindows ) :
   
          beginWin = windowNum * analysisObj.numDataPoints * analysisObj.stepSize
          endWin = beginWin + analysisObj.numDataPoints
  
          print ( " beginWin = " + str(beginWin) + " endWin = " + str(endWin))
  
          windowData = channelData [ int(beginWin) : int(endWin)]
  
          #print ( " windowData = " + str(windowData) )
  
          if len(windowData) == 0:
  
            break
  
          for taperIndex, taper in enumerate ( tapers ) :
   
              taperData = [float(a*b) for a,b in zip(windowData,taper)]
   
              fftData = fft(taperData,paddedNumDataPoints)
   
              spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
   
              print ( " padded num = " + str(paddedNumDataPoints) + " spectrum len = " + str(len(spectrumChannelData)))
   
              spectrumChannelData = list(spectrumChannelData[gridIndices])
   
              spectrumChannelData = (1 / float(analysisObj.samplingFrequency) ) * np.array(spectrumChannelData)
              spectrumChannelData = spectrumChannelData[analysisObj.lowerFrequency : analysisObj.upperFrequency+1]
   
              #print (spectrumChannelData)
   
              #print (" for taper = " + str(taperIndex ) + "  spectrumChannelData = " + str(spectrumChannelData) )
              spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)
    
          spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
    
          spectrogramData.append(spectrumChannelAvgData)
    
            #print (spectrogramData)
    
            #print (" for window = " + str(windowNum) + " spectrogramData = " + str(spectrogramData) )
    
            #break
    
        np.savetxt("outdata/channel_spectrogram_file_" + str(datafile.id) + "_channel_" + str(channelIndex) + ".txt", spectrogramData )
   
        np.savetxt("outdata/channel_spectrogram_file_" + str(datafile.id) + "_channel_" + str(channelIndex) + ".txt", spectrogramData )
        
        spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
        spectrumPSD = np.array(spectrumPSD)/100

        np.savetxt("outdata/channel_spectrum_PSD_file_" + str(datafile.id) + "_channel_" + str(channelIndex) + ".txt", spectrumPSD )
   
        plt.plot(spectrumPSD.transpose())
        plt.savefig("outdata/channel_spectrum_psd_transpose_file_" + str(datafile.id) + "_channel_" + str(channelIndex) + ".png" )
   
        plt.clf()
   
        plt.imshow(np.array(spectrogramData))
        plt.savefig("outdata/channel_spectrogram_file_" + str(datafile.id) + "_channel_" + str(channelIndex) + ".png" )
   
        plt.clf()
  
         #S = imag(S12x);
    
        ##WPLI
         #outsum   = nansum(S,1) # compute the sum; this is 1 x size(2:end)
         #outsumW  = nansum(abs(S),1) # normalization of the WPLI
         #outssq   = nansum(S.^2,1)
         #wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq) # do the pairwise thing in a handy way
    
         #break

  except:
    traceback.print_exc(file=sys.stdout)
  return

