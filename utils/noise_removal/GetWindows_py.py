#get artifacts, store results in file

import numpy as np
import os
import sys
import traceback
from scipy.io import loadmat
import math

#from scipy.signal import kaiser, group_delay, hilbert, angle

import matplotlib.pyplot as plt

def getWindows(k,minLen,winSafe):

    try:

        #%for every island of 0s extend it before and after (winSafe)
        if winSafe != 0:
            k[0] = 1 
            k[end] = 1
            #kDOWN = find(k(1:end-1) == 1 & k(2:end) == 0); 
            #kDOWN = kDOWN + 1;
            #kUP = find(k(1:end-1) == 0 & k(2:end) == 1); 
            for i in range ( len(kUP)):
                st = kDOWN(i)-winSafe
                if st < 1:
                    st=1
                ed = kUP(i)+winSafe
                if ed > length(k):
                    ed = length(k)
                k(st:ed) = 0
        
        #%allow only these which have enough length (minLen)
        if minLen != 0:
            #find islands of 1s
            k[0] = 0
            k[end] = 0
            kUP = find(k(1:end-1) == 0 & k(2:end) == 1)
            kUP = kUP + 1
            kDOWN = find(k(1:end-1) == 1 & k(2:end) == 0)
            d = kDOWN - kUP # %length of overshoots
            kd = d > minLen # %only select large crossings
            kUP = kUP(kd)
            kDOWN = kDOWN(kd)
            
            #redo boolean vector k
            k = zeros(1,length(k))
            for i in range ( len(kUP)):
                k(kUP(i):kDOWN(i)) = 1
        
        #%figure; hold on;
        #%plot(k2,'k*');
        #%plot(k,'r');
        #%return;
                
        windows = cat(2,kUP',kDOWN')
        signalOK = logical(k)

    except:
        
