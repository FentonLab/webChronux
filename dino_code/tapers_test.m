clear; close all; clc;

load('/Users/smitra/projects/webChronux/utils/257802-post_season_20161206_142914.mat');

eeg = eegData(14,:);

%remove mean
eeg = eeg - mean(eeg);

addpath('/Volumes/DATA/matlab/chronux/spectral_analysis/helper/');
addpath('/Volumes/DATA/matlab/chronux/spectral_analysis/continuous/');

eegFS = 250;

winLen = 3*eegFS;
params.fpass = [0,100]
params.Fs = eegFS;
params.tapers = [4,6]; %TW product, number of tapers (less than or equal to 2TW-1). 
params.trialave = 0;

[tapers,pad,Fs,fpass,err,trialave]=getparams(params);

nfft=max(2^(nextpow2(winLen)+pad),winLen);
[fx,findx]=getfgrid(Fs,nfft,fpass);
NF = (nfft/2)+1;
tap=dpsschk(tapers,winLen,Fs); % check tapers

figure; 
hold on;
for tI = 1:6 
    subplot(6,1,tI);
    plot(tap(:,tI));
    %colormap rgb;
end;
hold off;