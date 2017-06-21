clear; close all; clc;

load('EEG181.mat');

eeg1 = eegData(15,10*eegFS:100*eegFS);
eeg2 = eegData(16,10*eegFS:100*eegFS);

%remove mean
eeg1 = eeg1 - mean(eeg1);
eeg2 = eeg2 - mean(eeg2);

%addpath('/Volumes/DATA/matlab/Fieldtrip/');
addpath('/Users/smitra/projects/chronux_2_12/spectral_analysis/helper/');
addpath('/Users/smitra/projects/chronux_2_12/spectral_analysis/continuous/');

eegFS = 250;

winLen = 3*eegFS;

params.Fs = eegFS;
params.tapers = [4,9]; %TW product, number of tapers (less than or equal to 2TW-1). 
params.trialave = 0;

[tapers,pad,Fs,fpass,err,trialave]=getparams(params);

nfft=max(2^(nextpow2(winLen)+pad),winLen);
[fx,findx]=getfgrid(Fs,nfft,fpass);
NF = (nfft/2)+1;
tap=dpsschk(tapers,winLen,Fs); % check tapers

%split eeg1 and eeg2 into windows of length winLen
Nw = floor(length(eeg1)/winLen);

S12x = zeros(Nw,NF);
S1x = zeros(Nw,NF);
S2x = zeros(Nw,NF);

for wI = 1:Nw

    st = (wI-1)*winLen+1;
    ed = wI*winLen;

    eeg1x = eeg1(st:ed)';
    eeg2x = eeg2(st:ed)';

    J1=mtfftc(eeg1x,tap,nfft,Fs); %mt fourier
    J2=mtfftc(eeg2x,tap,nfft,Fs);

    J1=J1(findx,:,:); J2=J2(findx,:,:);
    S12=squeeze(mean(conj(J1).*J2,2)); %cross-spectrum

    S12x(wI,:) = S12';

    S1=squeeze(mean(conj(J1).*J1,2));
    S2=squeeze(mean(conj(J2).*J2,2));

    S1x(wI,:) = S1';
    S2x(wI,:) = S2';
end

S = imag(S12x);

%WPLI
outsum   = nansum(S,1);      % compute the sum; this is 1 x size(2:end)
outsumW  = nansum(abs(S),1); % normalization of the WPLI
outssq   = nansum(S.^2,1);
wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way

plot(fx,wpli)

