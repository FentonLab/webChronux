clear; close all; clc;

load('EEG181.mat');

eegFS = 250;

eeg = eegData(15,50*eegFS:80*eegFS);
eeg = eeg-mean(eeg);

Twindow = 1; %length of subwindow
Nwindow = Twindow * eegFS; %length of window in samples
Noverlap = round ( Nwindow /2); % 50% overlap
NFFT = 2^11; %resolution
[~,f] = pwelch(rand(1,length(eeg)), Nwindow, Noverlap, NFFT, eegFS); %get frequency vector (can be done other ways)

%pwelch
[Pxx,~] = pwelch(eeg, Nwindow, Noverlap, NFFT, eegFS);
Pxx = 20*log10(Pxx); %convert to dB

plot(f,Pxx)
