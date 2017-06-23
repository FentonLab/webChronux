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

def getNoiseSat(eeg,maxLengthCrossing):
    
    '''
    compute phase amplitude coupling, as in Vinnick
    PAC plotted at the end
    Original Matlab code by Dino Dvorak 2012 dino@indus3.net
    Python version by Siddhartha Mitra 2017 mitra.siddhartha@gmail.com
    Input: 
    Output: 
    '''

    try:    
        
        #%{
        #%testing
        #clear all; close all; clc;
        #load('sat1.mat');
        #eeg = x;
        #eegFS = 2000;
        #maxLengthCrossing = 5;
        #%}
        
        #get signal difference
        dif = np.diff(eeg);
        k = dif == 0;
        #add first sample (diff takes one sample out)
        k = [k(1), k]
        
        #skip extremely short saturations (1 value is 0), find 010
        k(1) = 0
        k(end) = 0
        #kUP = find(k(1:end-1) == 0 & k(2:end) == 1))
        #kUP = kUP + 1;
        #kDOWN = find(k(1:end-1) == 1 & k(2:end) == 0))
        d = kDOWN-kUP
        kd = find(d < 2) # allow two samples being same
        for i in range (len(kd)):
            ind = kd(i)
            k(kUP(ind):kDOWN(ind)) = 0
        
        k = ~k
        
        #%{
        #subplot(3,1,1);
        #plot(eeg);
        #axis tight;
        
        #subplot(3,1,2);
        #plot(dif);
        #axis tight;
        
        #subplot(3,1,3);
        #plot(k,'r');
        #axis tight;
        #%}
        
        #%{
        #%old method which required threshold specification
        #if bipolar == 1 %both positive and negative overshoot
            #k = eeg < threshold & eeg > -1*threshold;  %0 = bad
        #else
            #k = eeg < threshold; %only positive overshoot 
        #end
        #%}
                    
        #overshoots allowed
        if maxLengthCrossing != 0:                
            #first find overshoots and allow short ones
            k(1) = 1; k(end) = 1
            kDOWN = find(k(1:end-1) == 1 & k(2:end) == 0)
            kDOWN = kDOWN + 1
            kUP = find(k(1:end-1) == 0 & k(2:end) == 1)
            d = kUP - kDOWN + 1 # %length of overshoots
            kd = d > maxLengthCrossing # %only select large crossings
            kDOWN = kDOWN(kd)
            kUP = kUP(kd)
        
            #create new series and find long enough windows
            k = ones(1,length(k))
            for i in range ( len(kDOWN)):
                st = kDOWN(i)
                ed = kUP(i)
                k(st:ed) = 0
                
        signalOK = logical(k);

    except:
        traceback.print_exc(file=sys.stdout)    
        
    return       
