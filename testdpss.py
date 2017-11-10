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
from spectrum import *

def testDpss():

  try:
  
    samplingFrequency = 250
    timeBandWidth = 4
    timeWindow = 3 # time window in seconds
    
    numDataPoints =  timeWindow * samplingFrequency
  
    numTapers = 9
    [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )
    
    print(np.shape(tapers))
    
    plt.figure()
    plotnum = 411
    for index, taper in enumerate(tapers):
      print ( " taper num " + str(index) + " plot num " + str(plotnum))
      plt.subplot(plotnum)
      plt.plot(taper)
      plotnum += 1
      if index >=3:
        break
    plt.plot(eigenValues)
        
    
    plt.show()
    
    #plt.figure()                # the first figure
    #plt.subplot(11)             # the first subplot in the first figure
    #plt.plot([1, 2, 3])
    #plt.subplot(212)             # the second subplot in the first figure
    #plt.plot([4, 5, 6])    

  except:
    traceback.print_exc(file=sys.stdout)
    return

testDpss()
