#%finds pieces of eeg without saturation
#%uses threshold-free code based on diff=0
#%INPUT:
#%maxLengthCrossing - allows brief saturations such as spikes (in samples)
#%OUTPUT:
#%signalOK is boolean - 1s when signal is good
#%last edit Sept 1 2013

import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math
from scipy.signal import kaiser, group_delay, hilbert, angle
import matplotlib.pyplot as plt
