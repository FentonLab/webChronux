#get artifacts, store results in file
import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from scipy.signal import firwin, lfilter
import matplotlib.pyplot as plt
import pandas as pd
from multitaper import * 

def createSignal():
    
    try:

        x = np.linspace(0,360,1000)

        s1 = np.sin(np.array(x) * np.pi / 180 )
        
        s2 = np.sin(np.array(x) * np.pi / 18 ) # s2 = 10* w1
        
        s2 = np.array([x*y for x,y in zip(s1,s2)]) 
        
        #print ( s2)
        s3 = s1 + s2 

        #print ( len(s3))
        
        #To evaluate the performance of your code.
        
        #You can even add 60 Hz noise to test the hypothesis that the difference in the python and may lab examples is due effects of 60Hz:
        
        s60 = s3 + np.sin(np.array(x) * np.pi / 300 )
        print ( s60)

    except:
        
        traceback.print_exc(file=sys.stdout)    
        
    return s3
        

def runFTest():
    
    try:
        
        createSignal()
        
        #eegFS = 250 # sampling rate
        
        #eegData = loadmat('EEG181.mat')
        
        #eegData = eegData["eegData"]
                
        #eeg = eegData[14]
        
        #eeg = eeg[50 * eegFS - 1 : 80 * eegFS - 1 ]         
    
        #spectrumChannelList = []
        #fftDataList = []
        #taperFFTList = []
        #spectrumList = []
        
        #WINDOW_SIZE = 3
        
        #numDataPoints = WINDOW_SIZE * eegFS
        
        #paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( winLen ) ) + pad ) )

        #gridValues, gridIndices = getGridIndices(fpass[0], fpass[1], paddedNumDataPoints, eegFS)

        #NF = ( paddedNumDataPoints / 2 ) + 1
 
        #timeBandWidth = 4
        #numTapers = 9        
        
        #[tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )        
        
        #numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )        
        
        #for taperIndex, taper in enumerate(tapers):
            ## multiply by taper
            ## taperedData = [ float(a)*float(b) for a,b in zip( channelDataMap [ data ] , tapers [ taperIndex ] ) ]
            ##testData = numpy.ones( (1000), dtype=numpy.float64 )
            #taperData = [float(a*b) for a,b in zip(testData,taper)]
            #specgram(taperData , NFFT=paddedNumDataPoints, Fs=samplingFrequency, Fc=0, detrend=detrend, window = mlab.window_hanning, noverlap=250, cmap=None, xextent=None, pad_to=None, sides='default', scale_by_freq=None)
            #fftData = fft(taperData,paddedLen)
            #fftData = [float(a)/float(samplingFrequency) for a in fftData] 
            #fftDataList.append(fftData)
            #taperFFT = fft(taper,paddedLen)
            #taperFFTList.append(taperFFT)
            #spectrumChannel = fftData [ minFrequency : maxFrequency ] 
            #spectrumChannel = abs ( spectrumChannel )
            #spectrumChannelList.append ( spectrumChannel )
    
        #fValues = getFValues(fftDataList, taperFFTList, numTapers)
        #fittedData = [(a-b) for a,b in zip(testData, fValues)]

    except:
        
        traceback.print_exc(file=sys.stdout)    


def getFValues(fftDataList, taperFFTList, numTapers):

    try:
        evenFFT =  [x[:paddedLen/2] for i,x in enumerate(fftDataList) if i%2 == 0] 
        oddFFT =  [x[:paddedLen/2] for i,x in enumerate(fftDataList) if i%2 == 1] 
        # print " evenFFT " + str(evenFFT[0][1]) + " -- " + str(evenFFT[0][2]) 
        evenTaperSum = map ( sum, zip (  [x for i,x in enumerate(tapers) if i%2 == 0] ) )
        evenTaperSumSqList = []
        evenTaperSumList = []
    
        for taperSum in evenTaperSum:
            evenTaperSumSqListElem = []
            evenTaperSumSqListElem.append(taperSum)
            evenTaperSumSqListElem = evenTaperSumSqListElem*(paddedLen/2)
            evenTaperSumSqListElem = [ a*a for a in evenTaperSumSqListElem ]
            evenTaperSumSqList.append(evenTaperSumSqListElem)
    
            evenTaperSumListElem = []
            evenTaperSumListElem.append(taperSum)
            evenTaperSumListElem = evenTaperSumListElem*(paddedLen/2)
            evenTaperSumList.append(evenTaperSumListElem)
    
        evenTaperSumSq = map(sum,zip(*evenTaperSumSqList))
        taperEvenFFTProductList = [] 
    
        for index, evenTaperSumListElem in enumerate(evenTaperSumList):
            taperEvenFFTProduct = [a*b for a,b in zip(evenFFT[index],  evenTaperSumListElem)]
            taperEvenFFTProductList.append ( taperEvenFFTProduct )
    
        taperEvenFFTProduct = map(sum, zip (*taperEvenFFTProductList))
        amplitude = [a/b for a,b in  zip (taperEvenFFTProduct , evenTaperSumSq)  ]
    
        numEvenTapers = 0
        numOddTapers = 0
        if  numTapers%2 ==0:
            numEvenTapers = int(numTapers/2)
            numOddTapers = int(numTapers/2) -1
        else:
            numEvenTapers = int(ceil(numTapers/2) + 1)
            numOddTapers = int(ceil(numTapers/2))
        amplitudeList = []
        for index in range(numEvenTapers):
            amplitudeListElem = []
            amplitudeList.append ( amplitude )
        fittedFFTValueList = []
        for index, amplitudeValue in enumerate(amplitudeList):
            fittedFFTValue = [a*b for a,b in zip(amplitudeValue,  evenTaperSumList[index] ) ]
            fittedFFTValueList.append ( fittedFFTValue )
        amplitudeSquare = [ abs(a) * abs(a) for a in amplitude ]
        numerator =  [  ( numTapers - 1 )  * a * b for a,b in zip ( amplitudeSquare , evenTaperSumSq ) ] 
        denominatorList = []
        for numTaper in range(numEvenTapers):
            denominatorListElem =  [ abs ( a - b)* abs(a-b) for a,b in zip( evenFFT[numTaper], fittedFFTValueList[numTaper] ) ]  
            denominatorList.append(denominatorListElem)
        oddFFTList = []
        for numTaper in range(numOddTapers):
            oddFFTListElem =  [ abs ( a)* abs(a) for a in oddFFT[numTaper] ]  
            oddFFTList.append(oddFFTListElem)
        denominator1 = map(sum, zip(*denominatorList))
        denominator2 = map(sum,zip(*oddFFTList))
        denominator = [(a+b) for a,b in zip(denominator1, denominator2)]
        fValues =  [a/b for a,b in zip( numerator , denominator ) if  b != 0] 
        significance = robjects.r.qf ( 1-pValue , df1=2 , df2 = 2 * numTapers - 2 ) 
        evenTaperSumSqTaper = [numTapers * a for a in evenTaperSumSq]
        variance  =  [a/b for a,b in zip( denominator, evenTaperSumSqTaper ) if  b != 0] 
        stdDev = [sqrt(a) for a in variance] 
        #fValueObj = FValueObj(fValues, amplitude, significance)
        
    except:
        
        traceback.print_exc(file=sys.stdout)            

    return fValues, amplitude, significance

createSignal()
        