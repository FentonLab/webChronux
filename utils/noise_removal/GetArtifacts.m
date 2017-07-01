%get artifacts, store results in file
clear; close all; clc;

%input directory
%drIn = '/mnt/kant/dinod/ANALYSIS/basma/group1/MAT/';
%output directory
%drOut = '/mnt/kant/dinod/ANALYSIS/basma/group1/ART_border025dif/';

%if exist(drOut,'dir') == 0; mkdir(drOut); end;

doSaturation = 1; %DO saturation test (set to 1)
doSmallSig = 1;   %DO slow signal test
doHF = 1;         %DO high frequency artifact test

load('EEG181.mat');
eegFS = 250;
%eegFS = 2000;       %sampling rate
%eeg = eegData(15,:)
eeg = eegData(15,50*eegFS:80*eegFS);

%saturation 
maxLengthCrossing = 5; %maximum length of crossing in samples (allow spikes to overshoot)

%slow signal
thVar = 1e-4;       %variance threshold - DECREASE for LOWER sensitivity
thVarMin = 1e-2;    %minimum variance - DECREASE for LOWER sensitivity
thVarMax = 2;       %maximum variance
winLenVar = 0.05*eegFS; %minimum variance window under threshold
winLenVarSlide = 0.05*eegFS; %sliding variance window

%high frequency
lenHF = 50; %length of HF filter [samples]
freqHF = 250; %frequency of HF filter [Hz]
winLenVarHF = 0.1*eegFS; %variance window [samples]
winLenSmoothHF = 0.2*eegFS; %smoothing window after HF [samples]
thHF = 1e-5; %threshold for HF signal - INCREASE for LOWER sensitivity

%windows computation
minLen = 1*eegFS;   %minimum required length of good signal in seconds
winSafe = 0.25*eegFS;  %safe window around saturations to exclude

%list files
%files = dir([drIn '*.mat']);

%myPool = parpool('local',10);

%through files
% for fI = 1:length(files)
%     
%     %load MAT file
%     load([drIn files(fI).name]);
%     disp(files(fI).name);
%     
%     Nch = size(eegData,1);
%     
%     signalOK = cell(1,Nch);
%     windowsOK = cell(1,Nch);
%     
%     %through channels
%     %PP
%     parfor chI = 1:Nch

%eeg = eegData(chI,:);

%if isempty(eeg); continue; end;

sigSat = [];
sigSmall = [];
sigHF = [];

%get saturations (non-changing signal)
if doSaturation == 1
    sigSat = GetNoiseSat(eeg,maxLengthCrossing);
end

%get HF noise
if doHF == 1
    sigHF = GetNoiseHF(eeg,eegFS,lenHF,freqHF,winLenVarHF,winLenSmoothHF,thHF);
end

%merge all three signals
k = ones(1,length(eeg));
if doSaturation == 1; k = k & sigSat; end;
if doSmallSig == 1; k = k & sigSmall; end;
if doHF == 1; k = k & sigHF; end;

%compute windows
[sigOKCH, winOKCH] = GetWindows(k,minLen,winSafe);

%{
figure
subplot(5,1,1);
plot(eeg,'k');
title('eeg');
axis tight;

subplot(5,1,2);
plot(sigSat,'r');
title('saturations');
axis tight;

subplot(5,1,3);
plot(sigSmall,'r');
title('small signal');
axis tight;

subplot(5,1,4);
plot(sigHF,'r');
title('hf signal');
axis tight;

subplot(5,1,5); hold on;
plot(k,'k*');
plot(sigOKCH,'r');
title('good signal');
axis tight;
%}

%store data
%signalOK{chI} = sigOKCH;
%windowsOK{chI} = winOKCH;
        
%     end
%     
%     %save data
%     save([drOut files(fI).name], 'signalOK','windowsOK');
%     
% end %files

delete(myPool);
