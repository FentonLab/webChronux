import pyedflib
from scipy import *
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import traceback
from scipy.io import savemat, loadmat

def analyzeEDFData(filePath):

  try:

    print ( " in analyze data ")
    f = pyedflib.EdfReader(filePath)
    n = f.signals_in_file
  
    signal_labels = f.getSignalLabels()

    dataMatrix = []
    
    for channelIndex in range(n):
  
      channelData = f.readSignal(channelIndex)
      dataMatrix.append(channelData)
      
    savemat("/Users/smitra/projects/webChronux/edfsave.mat", {"eegdata":dataMatrix})

  except:
    traceback.print_exc(file=sys.stdout)
    return

analyzeEDFData("/Users/smitra/self/andre/MIT-concussion/baseline257802_20161018_194838.edf")
print (" end 999 ")
