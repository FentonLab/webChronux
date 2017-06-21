import pyedflib
from scipy import *
import numpy as np
import pandas as pd
import sys
import os
import traceback
import shutil
from multitaper import *
import matplotlib.pyplot as plt
from scipy.io import savemat

samplingFrequency = 250
upperFrequency = 70
lowerFrequency = 1
timeBandWidth = 2
numTapers = 3
timeWindow = 5 # time window in seconds
STEP_SIZE = 2 # in seconds
numDataPoints =  timeWindow * samplingFrequency
stepSize = STEP_SIZE * samplingFrequency
padding = 1
paddedNumDataPoints = int ( pow (2, ceil ( math.log( numDataPoints, 2 )) + padding ))

filePath = "/Users/smitra/self/andre/testdata"
ofilePath = "/Users/smitra/self/andre/mat"

ANALYZE_CHANNELS = ["FP2", "O2", "O1", "FP1", "F7", "F8", "T5", "T6", "T4", "A2", "C3","C4","T3","Cz","Fz"]
    
def saveFileData(fileName, opath):

    try:

        print ( " file = " + str(fileName) )
        
        f = pyedflib.EdfReader(fileName)
        n = f.signals_in_file
      
        signal_labels = f.getSignalLabels()
        
        channelDataList = []
      
        for channelIndex in range(n):
            
            if signal_labels[channelIndex] not in ANALYZE_CHANNELS:
                
                continue
            
            #print ( " channel is = " + str(signal_labels[channelIndex]) ) 
        
            channelData = f.readSignal(channelIndex)
        
            channelDataList.append(np.array(channelData))
            
            print (" data len for " + str(signal_labels[channelIndex]) + " is = " + str(len(channelData)))
            
            if signal_labels[channelIndex] == "FP1":
                print ( " !!!!!!!!!!!! in fp1 ")
                plt.plot(channelData)
                plt.show()
                break

        #channelDataList = np.array(channelDataList).transpose()
        
        #print ("min = " + str(min([min(x) for x in channelDataList])))
        #print ("max = " + str(max([max(x) for x in channelDataList])))

        ##print ( " !!!!!! shape is = " + str(np.shape(np.array(channelDataList))) )  

        #savemat(opath + "_outmat.mat", mdict={'arr': channelDataList} )
            
       
    except:
        traceback.print_exc(file=sys.stdout)
        return
            

def saveData():

    try:
        
        files = os.listdir(filePath)
        for filex in files:
            if filex.find(".edf") != -1:
                saveFileData(filePath + "/" + filex, ofilePath   + "/" + filex)

    except:
        traceback.print_exc(file=sys.stdout)
    return

saveData()
print (" end 999 ")
