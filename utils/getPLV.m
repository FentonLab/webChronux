clear; close all; clc;

load('EEG181.mat');

eegFS = 250;

eeg1 = eegData(15,50*eegFS:80*eegFS);
eeg2 = eegData(16,50*eegFS:80*eegFS);
eeg1 = eeg1-mean(eeg1);
eeg2 = eeg2-mean(eeg2);

%eegFS = 2000; %sampling rate
    
%create filters
freqs = 5:5:80;
filters = cell(length(freqs),2);

bw = 4; %bandwidth
nyq = eegFS/2;
attendB = 40; %attenuation
attenHz = 2; %transition band
%gamma freqs
for fI = 1:length(freqs)

    lcut = freqs(fI)-bw/2;
    hcut = freqs(fI)+bw/2;
    
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
    filters{fI,1} = b;
    
    %group delay
    [a,f] = grpdelay(b,1,nyq,eegFS);
    k = f >= lcut & f <= hcut;
    gd = fix(mean(a(k)));
    filters{fI,2} = gd;
end

PLVs = zeros(1,length(freqs));
    
%filters
for fI = 1:length(filters)
    b = filters{fI,1};
    gd = filters{fI,2};

    ph1 = filter(b,1,eeg1);
    ph1 = [ph1(gd+1:end) zeros(1,gd)];
    ph1 = angle(hilbert(ph1));

    ph2 = filter(b,1,eeg2);
    ph2 = [ph2(gd+1:end) zeros(1,gd)];
    ph2 = angle(hilbert(ph2));

    %phase diff
    phiD = ph1 - ph2;    

    %PLV
    PLVs(fI) = (1/length(phiD))*(abs(sum(exp(1i*phiD)))); 

end

plot(freqs,PLVs)