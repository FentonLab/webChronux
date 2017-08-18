clear; close all; clc;

load('EEG181.mat');

eegFS = 250;

eeg = eegData(15,50*eegFS:80*eegFS);
eeg = eeg-mean(eeg);

%edges for phase
edges = linspace(-pi,pi,21);

%eegFS = 2000; %sampling rate
%winLen = length(eeg); %window length

x = linspace(0,200,5000);

s1 = sin( x * pi / 180 )

s2 = sin( x * pi / 18 ) % s2 = 10* w1

s2 = s1.*s2 

%print ( s2)
s3 = s1 + s2 

%eeg = s3
% figure(1);
% hold;
% plot(s1);
% 
% plot(s2);
% plot(s3);
% hold;

winLen = length(eeg)

%create alpha filter
lcut = 9;
hcut = 13;
attendB = 40;
attenHz = 2;
nyq = eegFS/2;
Fstop1 = (lcut - attenHz)  / nyq;  
Fpass1 = lcut  / nyq; 
Fpass2 = hcut / nyq;
Fstop2 = (hcut + attenHz) / nyq;
Astop1 = attendB;
Apass  = 1;
Astop2 = attendB;
h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                Fpass2, Fstop2, Astop1, Apass, Astop2);
Hd = design(h, 'kaiserwin');
bAlpha = Hd.Numerator;
%group delay
[a,f] = grpdelay(bAlpha,1,nyq,eegFS);
k = f >= lcut & f <= hcut;
gdAlpha = fix(mean(a(k)));
    
   
%create multiple gamma filters
fGamma = 20:5:100;
filtersGamma = cell(length(fGamma),2);

bw = 20; %bandwidth
attendB = 40; %attenuation
attenHz = 4; %transition band
%gamma freqs
for fI = 1:length(fGamma)

    lcut = fGamma(fI)-bw/2;
    hcut = fGamma(fI)+bw/2;
    
    Fstop1 = (lcut - attenHz)  / nyq;  
    Fpass1 = lcut  / nyq; 
    Fpass2 = hcut / nyq;
    Fstop2 = (hcut + attenHz) / nyq;
    Astop1 = attendB;
    Apass  = 1;
    Astop2 = attendB;
    h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                    Fpass2, Fstop2, Astop1, Apass, Astop2);
    Hd = design(h, 'kaiserwin');
    b = Hd.Numerator;
    filtersGamma{fI,1} = b;
    
    %group delay
    [a,f] = grpdelay(b,1,nyq,eegFS);
    k = f >= lcut & f <= hcut;
    gd = fix(mean(a(k)));
    filtersGamma{fI,2} = gd;
end
    

MIs = zeros(1,length(fGamma));

%phase (alpha)
ph = filter(bAlpha,1,eeg);
ph = [ph(gdAlpha+1:end) zeros(1,gdAlpha)];
ph = angle(hilbert(ph));

%amplitude for gamma range + MI
for fI = 1:length(fGamma)
    
    b = filtersGamma{fI,1};
    gd = filtersGamma{fI,2};
    
    %amplitude
    amp = filter(b,1,eeg);
    amp = [amp(gd+1:end) zeros(1,gd)];
    amp = abs(hilbert(amp));
    
    %compute raw modulation index
    MI = getMI(amp,ph,edges);

    MIs(fI) = MI;
end

plot(fGamma,MIs);
