% %% load things
clear all; close all; clc;  warning off;
baseFolder = '\\rackstation\backup\Ellie\Int Juxta\032819_dlx32\boutFinder\';
f = 'c1_1';
fullRecord = abfload(['\\rackstation\backup\Ellie\Int Juxta\032819_dlx32\c1_0001.abf']);
load([baseFolder f '\MotifTimes.mat']);
fs = 50000;
spikes = fullRecord(:,1); 
load([baseFolder f '\thresh.mat'])

%% calculate noise
noiseStart = 1; %CHANGE
noiseStop = 2; %CHANGE
noiseWindow = noiseStart*fs:noiseStop*fs-1;
noiseVals = spikes(noiseWindow);
avgNoise = mean(noiseVals);
sdNoise = std(noiseVals);

%% calculate spikes
window = .002 * fs; %spike width in sec
flipped = -1 * spikes;
%flipped = spikes;
[vals, locs] = findpeaks(flipped, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
rasterSpikes = zeros(length(spikes),1);
for i = 1:length(locs)
    rasterSpikes(locs(i)) = 1;
end

spikeAmps =[];
for i=1:length(rasterSpikes)
    if rasterSpikes(i) == 1
       % amp = flipped(i)-(avgNoise*-1);
        amp = flipped(i)-avgNoise;
        spikeAmps = [spikeAmps; amp];
    end
end

avgSignal = mean(spikeAmps);

%% calculate SNR
snr = avgSignal/sdNoise;
disp('signal noise snr')
snrStats = [avgSignal sdNoise snr]
save snrStats.mat snrStats

%% calculate spontaneous FR 
pseudoSpikes = spikes;
starts = Motif.start;
stops = Motif.stop;
for i = 1:length(starts)
    pseudoSpikes(starts(i):stops(i))=nan;
end 
pseudoSpikes = pseudoSpikes(~isnan(pseudoSpikes));
%pseudoFlipped = pseudoSpikes;
pseudoFlipped = -1 * pseudoSpikes;
[pvals, plocs] = findpeaks(pseudoFlipped, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
pseudoRasterSpikes = zeros(length(pseudoSpikes),1);
for i = 1:length(plocs)
    pseudoRasterSpikes(plocs(i)) = 1;
end
avgSpontFR = mean(sum(pseudoRasterSpikes)/(length(pseudoRasterSpikes) / fs))

%% stability test
wing = 10 * fs;
wingSpikes = [];
wingBeg = (starts*fs)-wing;
wingEnd = (starts*fs)-1;
sumWingSpikes = zeros(length(starts), 1);
avgWingFRs = zeros(length(starts), 1);
for i = 1:length(starts)
    wingSpikes(:,i) = spikes(wingBeg(i):wingEnd(i));
    wingSpikes = -1 * wingSpikes;
    [wvals, wlocs] = findpeaks(wingSpikes(:,i), 'MinPeakHeight', thresh, 'MinPeakDistance', window);
    sumWingSpikes(i) = length(wlocs);
    avgWingFRs(i) = length(wlocs)/wing*fs;
end
avgWingFRs

%%