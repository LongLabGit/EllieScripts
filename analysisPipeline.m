%uses motif finder start/stop times to display ephys on playback trials
%% load shit
clear all; close all; clc;  warning off;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\';
boutFolder = 'boutFinder\';
motifFolder = 'motifFinder\';
f = 'c1_1';
fullRecord = abfload([baseFolder boutFolder f '\' f '_filteredspikes.abf']);
fs = 50000;
wing = 5 * fs;
thresh = .2;
load([baseFolder boutFolder f '\' 'MotifTimes.mat']);
Bout = Motif;
clear Motif;
load(strcat(baseFolder, motifFolder, f, '\MotifTimes.mat'));
load([baseFolder boutFolder f '\' 'fixedSpikes.mat']);
goodDataTimes = 1:(1340*fs);
audio = fullRecord(goodDataTimes,3);
spikes = fixedSpikes;
bout = audioread([baseFolder boutFolder 'bout.wav']);
motif = audioread([baseFolder motifFolder 'motif.wav']);
eyeTrig = fullRecord(goodDataTimes,4);
%% align spikes by bouts
lengthRecording = length(spikes)/fs/60;
starts = Bout.start * fs;
stops = Bout.stop * fs;
stopLocs = find(stops < (goodDataTimes(end)-wing));
starts = starts(1:length(stopLocs));
stops = stops(1:length(stopLocs));
numBouts = length(stops);
boutSpikes = {numBouts};
lengthBouts = zeros(numBouts, 1);
boutAudio = {numBouts};
for i = 1:numBouts
    boutSpikes{:,i} = spikes((starts(i)-wing):(stops(i)+wing));
    lengthBouts(i) = length(boutSpikes{:,i});
    boutAudio{:,i} = audio((starts(i)-wing):(stops(i)+wing));
end
for i = 1:numBouts
    if length(boutSpikes{:,i}) > min(lengthBouts)
       vec = boutSpikes{:,i};
       difference = length(vec)-min(lengthBouts);
       fixedVec = vec(1:end-difference);
       boutSpikes{:,i} = fixedVec;
    end
end
%% plot aligned spikes 
window = .002 * fs; %spike width in sec
rasterLocs = [];
numPeaks = [];
valPeaks = [];

for i = 1:numBouts
    flipped = -1 * boutSpikes{:,i};
    [vals, locs] = findpeaks(flipped, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
    rasterLocs = [rasterLocs; locs];
    numPeaks = [numPeaks; length(locs)];
    valPeaks = [valPeaks; vals];
    boutRaster = zeros(length(boutSpikes{:,i}),1);
    for j = 1:length(locs)
        boutRaster(locs(j)) = 1;
    end
    rasterMat{:,i} = boutRaster;
end

figure;
%subplot 1
axis(1) = subplot(2,1,1);
title('spikes aligned by playback bouts')
plot(audio(starts(1)-wing:stops(1)+wing), 'k')
lenAudio = audio((starts(1)-wing):(stops(1)+wing));
pad = length(lenAudio) - length(bout);
if mod(pad, 2) == 1
    pad = pad-1;
end
newSpec = [zeros(pad/2,1); bout; zeros(pad/2,1)];
vigiSpec(newSpec, fs, [], .45);
set(axis,'xtick',[])       
set(axis, 'ytick', [])

%subplot 2
xAxis = linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i}));
for i = 1:numBouts
    axis(2) = subplot(2,1,2);
    if mod(i,5) == 0
        plot(xAxis,boutSpikes{:,i}-1*i, 'r')
    else
        plot(xAxis,boutSpikes{:,i}-1*i, 'k')
    end 
    hold on
end
set(axis,'TickDir','out')    
set(axis, 'ytick', [])
linkaxes(axis, 'x') 

%% plot rasters
%subplot 1
figure;
axis(1) = subplot(3,1,1);
title('raster plots aligned by playback bouts')
vigiSpec(newSpec, fs, [], .5)
set(gca,'xtick',[])       
set(gca, 'ytick', [])

%subplot 2
space = 3;
for j = 1:numBouts
   axis(2) = subplot(3,1,2);
   for k = 1:length(rasterMat{:,j})
       logicals = rasterMat{:,j};
       if logicals(k) == 1
           if mod(j,5) == 0
              plot([k/fs, k/fs], [0, 3]-space*j,'Color', 'r')
           else
              plot([k/fs, k/fs], [0, 3]-space*j,'Color', 'k')
           end
       end
   end
   hold on
end
%set(gca,'xtick', [])       
ylim([-space*numBouts-5, 5])
set(gca,'TickDir','out')  
set(gca, 'ytick', [])

% make psth
bin = .02 * fs;
p = 1;
psthMat = [];
for m = 1:bin:length(boutSpikes{:,1})-bin %not all of it being sampled btw!!!
    sumAll = 0;
    for n = 1:numBouts
        cm = rasterMat{n};
        sumSingle = sum(cm(m:m+bin-1))/bin*fs;
        sumAll = sumAll+sumSingle;
    end
    sumAll = mean(sumAll);
    psthMat(p) = sumAll;
    p=p+1;
end
leftWingAvgFR = psthMat(1:round(wing/bin));

%subplot 3
axis(3) = subplot(3,1,3);
plot(linspace(0, length(cm)/fs, length(psthMat)),psthMat, 'k', 'LineWidth', 2)
xlabel('Time (s)') 
ylabel('Frequency')
set(gca,'TickDir','out')    
linkaxes(axis, 'x') 
%% check eyeball data
rebound = .3*fs;
fsEyes=2;
[valEye,locEyeFrames] = findpeaks(eyeTrig, 'minPeakHeight', 2, 'minPeakDistance', rebound);
numEyeFrames = length(locEyeFrames);
load([baseFolder 'eyeball\' f '\' 'openLocs.mat']);
load([baseFolder 'eyeball\' f '\' 'numImgs.mat']);
goodOpenFrames = find((openLocs/fsEyes*fs) < goodDataTimes(end));
openLocs = openLocs(1:length(goodOpenFrames));

%% eye plots
openEyeFrames = locEyeFrames(openLocs);
eyeRaster = zeros(numImgs,1);
for i = 1:length(openLocs)
  ind = openLocs(i);
  eyeRaster(ind) = 1;
end
xAxEyes = linspace(0, length(eyeRaster)/fsEyes, length(eyeRaster));
figure; 
title('open eye frames')
subplot(2,1,1)
plot(xAxEyes, eyeRaster); shg;
ylim([0 1.5])
xlim([0 length(eyeTrig)/fs])
subplot(2,1,2);
plot(xAxEyes, movmean(eyeRaster, 20)); shg;
xlabel('Time(s)');
xlim([0 length(eyeTrig)/fs])

%% motif alignment

%bout motif analysis and plots
mwing = 0*fs;
newAudio = [zeros(mwing,1); motif; zeros(mwing,1)];
% templateLength = Motif(1).stop(1) - Motif(1).start(1);
templateLength = length(newAudio)/fs;
mRasterLocs = [];
mNumPeaks = [];
mValPeaks = [];
countAll = 0;

figure;
%subplot 1
title('motifs within bouts')
axis(1) = subplot(2,1,1);
plot(linspace(0, length(newAudio)/fs, length(newAudio)), newAudio, 'k')
set(gca,'xtick',[])       
set(gca, 'ytick', [])
%subplot 2
axis(2) = subplot(2,1,2);
for i = 1:numBouts
    for j = 1:length(Motif(i).start)
       currentLength = Motif(i).stop(j) - Motif(i).start(j);
       warp = templateLength/currentLength;
       motifStart = round(Motif(i).start(j)*fs);
       motifStop = round(Motif(i).stop(j)*fs);
       tmp = rasterMat{:,i};
       motifSpikes = tmp(motifStart:motifStop);
       for k = 1:length(motifSpikes)
          if motifSpikes(k) == 1
            plot([k/fs*warp, k/fs*warp], [0, 2]-2*countAll,'Color', 'k')
          end
       end
      hold on 
      countAll = countAll + 1;
      end
end   
set(gca,'TickDir','out')    
set(gca, 'ytick', [])
linkaxes(axis, 'x') 

%% motifs only analysis and plots
load([baseFolder motifFolder f '\motifInterval1.mat']);
load([baseFolder motifFolder f '\motifInterval2.mat']);
spikesm1 = spikes(motifInterval1);
spikesm2 = spikes(motifInterval2);
moSpikes = {spikesm1 spikesm2};
MotifOnly = [Motif(26); Motif(27)];
mrasterLocs = [];
mnumPeaks = [];
mvalPeaks = [];
    
for i = 1:length(MotifOnly)
    flipped = -1 * moSpikes{:,i};
    [mvals, mlocs] = findpeaks(flipped, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
    mrasterLocs = [mrasterLocs; mlocs];
    mnumPeaks = [mnumPeaks; length(mlocs)];
    mvalPeaks = [mvalPeaks; mvals];
    moRaster = zeros(length(moSpikes{:,i}),1);
    for j = 1:length(mlocs)
        moRaster(mlocs(j)) = 1;
    end
    morasterMat{:,i} = moRaster;
end

figure;
%subplot 1
title('isolated motifs')
axis(1) = subplot(2,1,1);
plot(linspace(0, length(newAudio)/fs, length(newAudio)), newAudio, 'k')
set(gca,'xtick',[])       
set(gca, 'ytick', [])

%subplot 2
axis(2) = subplot(2,1,2);
for i = 1:length(moSpikes)
    for j = 1:length(MotifOnly(i).start)
       currentLength = MotifOnly(i).stop(j) - MotifOnly(i).start(j);
       warp = templateLength/currentLength;
       motifStart = round(MotifOnly(i).start(j)*fs);
       motifStop = round(MotifOnly(i).stop(j)*fs);
       tmp = morasterMat{:,i};
       motifSpikes = tmp(motifStart:motifStop);
       for k = 1:length(motifSpikes)
          if motifSpikes(k) == 1
            plot([k/fs*warp, k/fs*warp], [0, 2]-2*countAll,'Color', 'k')
          end
       end
       hold on 
       countAll = countAll + 1;
    end 
end
set(axis,'TickDir','out')    
set(axis, 'ytick', [])
linkaxes(axis, 'x') 
