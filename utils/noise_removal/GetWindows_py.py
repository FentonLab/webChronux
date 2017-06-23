#%generates good signal and list of good windows
#%based on binary signal where 1=good signal
#%INPUT:
#%k - good signal
#%minLen - minimum length of good window of data (in samples)
#%winSafe - adds safe boarder around all saturations (in samples)
#%OUTPUT:
#%signalOK is boolean - 1s when signal is good
#%windows is list of windows starts and ends, which have good signal
#%last edit Sept 1 2013

def getWindows(k,minLen,winSafe)

#%{
#clear all; close all; clc;
#load('sample data/k.mat');
#eegFS = 1000;
#minLen = 0.5*eegFS;
#winSafe = 0.5*eegFS;
#k2 = k;
#%}

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

    
