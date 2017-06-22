import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
#from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt

# get frequency grid
def getGridIndices(lowerFrequency, upperFrequency, samplingFrequency, padding):

  try:

    if upperFrequency and lowerFrequency:

        paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( numDataPoints ) ) + padding ) )

        #print ( " paddedNumDataPoints = " + str(paddedNumDataPoints))

        frequencyResolution = float ( samplingFrequency ) / float ( paddedNumDataPoints )
        
        gridValues = np.arange ( 0, paddedNumDataPoints , frequencyResolution )

        #print ( " frequencyResolution = " + str(frequencyResolution))
        
        #print ( " gridValues = " + str(gridValues))
        
        gridIndices = np.where ( (gridValues >= lowerFrequency ) & (gridValues <= upperFrequency ) )

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

def getWPLI():
    
    '''
    compute phase lag variance
    PLV plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    
    '''
    
    WINDOW_SIZE = 3
    
    try:    

        eegData = loadmat('EEG181.mat')
        
        eegFS = 250
        
        eeg1 = eegData[14]
        eeg2 = eegData[15]        

        eeg1 = [ x[50*eegFS:80*eegFS] for x in eeg1] 
        eeg1 = [ x[50*eegFS:80*eegFS] for x in eeg2] 

        eeg1 = eeg1 - np.mean(eeg1, axis = 0)        
        eeg2 = eeg2 - np.mean(eeg2, axis = 0)        
        
        winLen = WINDOW_SIZE * eegFS
        
        tapers = [4,9]
        pad = 1
        Fs = 250
        fpass = [1,100]
        timeBandWidth = 2
        
        paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + winlen ) )

        lowerFrequencyGrid, upperFrequencyGrid, gridIndices = getGridIndices(lowerFrequency, upperFrequency, samplingFrequency, padding)
        
        NF = ( paddedNumDataPoints / 2 ) + 1
 
        [tapers, eigenValues] = dpss_windows(int(winLen), float(bandWidth), int(numTapers) )
        
        Nw = np.floor(len(eeg1)/winLen)
        
        S12x = np.zeros(Nw,NF)
        S1x = np.zeros(Nw,NF)
        S2x = np.zeros(Nw,NF)
        
        for wI in range ( Nw ) :
        
            start = (wI-1) * winLen + 1
            end = wI * winLen
        
            eeg1x = eeg1[start:end]
            eeg2x = eeg2[start:end]
        
            S = performCrossSpectra(eeg1x,eeg2x,tapers,paddedNumDataPoints,samplingFrequency,gridIndices,timeBandWidth, fpass) # FFT

            outsum   = nansum(S,1)      # compute the sum; this is 1 x size(2:end)
            outsumW  = nansum(abs(S),1) # normalization of the WPLI
            outssq   = nansum(S.^2,1)
            wpli     = (outsum**2 - outssq)* np.linalg.inv(outsumW**2 - outssq) #do the pairwise thing in a handy way
            
            plot(fx,wpli)            

@app.task
def performCrossSpectra(eeg1x,eeg2x,tapers,paddedNumDataPoints,samplingFrequency,gridIndices,timeBandWidth, fpass) :
  
  try:
  
      print ( " in analyze data ")
    
      signal_labels = f.getSignalLabels()
    
      print ( " num signals = " + str(n) )
      print ( " signal labels = " + str(signal_labels) )
    
      beginWin = 0
      endWin = 0
    
      numTapers = 2 * timeBandWidth -1
      
      windowData = {}
      
      lowerFrequency = fpass[0]
      upperFrequency = fpass[1]
      
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
   
        #for windowNum in range ( numWindows ) :
   
          #beginWin = windowNum * analysisObj.numDataPoints * analysisObj.stepSize
          #endWin = beginWin + analysisObj.numDataPoints
  
          #print ( " beginWin = " + str(beginWin) + " endWin = " + str(endWin))
  
          #windowData = channelData [ int(beginWin) : int(endWin)]
  
          ##print ( " windowData = " + str(windowData) )
  
          #if len(windowData) == 0:
  
            #break
  
        for taperIndex, taper in enumerate ( tapers ) :
 
            taperData1 = [float(a*b) for a,b in zip(eeg1x,taper)]
            fftData1 = fft(taperData1,paddedNumDataPoints)
            fftData1 = [fftData1[x] for x in gridIndices]
 
            taperData2 = [float(a*b) for a,b in zip(eeg2x,taper)]
            fftData2 = fft(taperData2,paddedNumDataPoints)
            fftData2 = [fftData2[x] for x in gridIndices]
            
            crossSpectrumChannelData = np.array ([x*conj(y) for x,y in zip ( fftData1, fftData2) ])            
        
            crossSpectrumChannelData = [x.imag for x in crossSpectrumChannelData]

        crossSpectrumChannelAvgData = np.mean(crossSpectrumChannelData)
        
        windowData[channelNum] = crossSpectrumChannelAvgData

  except:
    traceback.print_exc(file=sys.stdout)
  return
        

