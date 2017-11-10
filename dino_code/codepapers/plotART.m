%plots events and eeg using EEGPlot
clear all; close all; clc;

addpath('/Volumes/DATA/matlab/eeglab/functions/sigprocfunc/');

dr = 'MAT/';
drART = 'ART/';
%drART = 'ART_sat/';
%drART = 'ART_low/';
%drART = 'ART_hf/';


fI = 1;
files = dir([dr '*.mat']);

load([dr files(fI).name]);
load([drART files(fI).name]);

chI = 1;

sigOK = signalOK{chI};

sum(sigOK)/length(sigOK)

eeg = eegData(chI,:);
eeg = eeg/max(eeg);

%eegplot(eeg,'srate', eegFS, 'winlength', 50, 'spacing', 1 );

k = ~sigOK;

%plot events
k(1) = 0; k(end) = 0;
kUP = find(k(1:end-1) == 0 & k(2:end) == 1); kUP = kUP + 1;
kDOWN = find(k(1:end-1) == 1 & k(2:end) == 0);

%build events - start and duration in SAMPLES !!
tp = cell(1,length(kUP));
dur = cell(1,length(kUP));
start = cell(1,length(kUP));
for i = 1:length(kUP)
   tp(1,i) = {''};
   dur(1,i) = {kDOWN(i)-kUP(i)};
   start(1,i) = {kUP(i)};
end

events = struct('type',tp,'duration',dur,'latency',start) ;
eegplot(eeg,'events', events, 'srate', eegFS, 'winlength', 50, 'spacing', 1, 'ploteventdur','on' );
