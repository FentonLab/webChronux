clear all; close all; clc;

load('EEG181.mat');

eegFS = 250;
wFactor = 8;
bands = 2:60;
wavelets = getWaveletsNorm( bands, wFactor, eegFS );

eeg = eegData(15,50*eegFS:80*eegFS);

pwrs = zeros(1,length(wavelets));

%filter
for filtI = 1:length(wavelets)

    psi = wavelets{filtI};

    %convolution of wavelet and signal
    c = conv(eeg,psi);
    %fix start and end
    N = round((length(psi)-1)/2);
    c = c(N:length(c)-N); 
    if length(c) > size(eeg,2)
        c = c(1:size(eeg,2));
    end
    power = (abs(c)).^2;

	pwrs(filtI) = mean(power);
end

plot(bands,pwrs)