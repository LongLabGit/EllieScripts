
%uses motif finder start/stop times to display ephys on playback trials

%% load shit
clear all; close all; clc; 
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\boutFinder\';
f = 'c1_1';
load([baseFolder f '\' 'MotifTimes.mat']);
fullRecord = abfload([baseFolder f '\' f '_filteredspikes.abf']);
spikes = fullRecord(:,1); 
audio = fullRecord(:,3);
bout = audioread([baseFolder 'bout.wav']);
load([baseFolder  f  '\thresh.mat']);
%audiowrite([baseFolder 'spikes.wav'], spikes, fs);
fs = 50000;
%% plot aligned spikes to bout trials
lengthRecording = length(spikes)/fs
starts = Motif.start * fs; 
stops = Motif.stop * fs;
numBouts = length(stops);
wing = 100000;

for i = 1:numBouts
    boutSpikes{:,i} = spikes(starts(i)-wing:stops(i)+wing);
end

% figure(20);
% ax(1) = subplot(2,1,1);
% %plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})), audio(starts:stops), 'k')
%  lenAudio = audio(starts(1)-wing:stops(1)+wing);
%  pad = length(lenAudio) - length(bout);
%  newSpec = [zeros(pad/2,1); bout; zeros(pad/2,1)];
%  
%   for i=1:pad/2
%      if newSpec(i) == 0
%          newSpec(i) = .2;
%      end
%      if newSpec(i+length(bout)+pad/2) == 0
%          newSpec(i+length(bout)+pad/2) = .2;
%      end
%   end
%  
% 
%  vigiSpec(newSpec, fs, [], .6);
% 
%   set(gca,'xtick',[])       
%   set(gca, 'ytick', [])
% %end
% 
% 
% %plot(audio(starts(1)-wing:stops(1)+wing), 'k')
% %set(gca,'xtick',[])       
% %set(gca, 'ytick', [])
% hold on
% 
% for i = 1:numBouts
%     ax(2) = subplot(2,1,2);
%     %plot(boutSpikes{:,i}, 'k')
%     plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})),boutSpikes{:,i}-3*i, 'k', 'LineWidth', .5);
%     set(gca,'TickDir','out')       
%     set(gca, 'ytick', [])
%     hold on 
% end
% 
% linkaxes(ax,'x');

%% raster and PSTH maker
window = .002 * fs; %spike width in sec
rasterLocs = [];
numPeak = [];
valPeaks = [];

for i = 1:numBouts
   %flipped = abs(boutSpikes{:,i});
    flipped = -1 * boutSpikes{:,i};
    [vals, locs] = findpeaks(flipped, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
    rasterLocs = [rasterLocs; locs];
    numPeak = [numPeak; length(locs)];
    valPeaks = [valPeaks; vals];
    boutRaster = zeros(length(boutSpikes{:,i}),1);
    for j = 1:length(locs)
        boutRaster(locs(j)) = 1;
    end
    rasterMat{:,i} = boutRaster;
end

figure;
%plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})), audio(starts:stops), 'k')
 lenAudio = audio(starts(1)-wing:stops(1)+wing);
 pad = length(lenAudio) - length(bout);
 newSpec = [zeros(pad/2,1); bout; zeros(pad/2,1)];
 
%   for i=1:pad/2
%      if newSpec(i) == 0
%          newSpec(i) = .2;
%      end
%      if newSpec(i+length(bout)+pad/2) == 0
%          newSpec(i+length(bout)+pad/2) = .2;
%      end
%   end
  axis(1) = subplot(3,1,1);
  vigiSpec(newSpec, fs, [], .6);
  set(gca,'xtick',[])       
  set(gca, 'ytick', []) 

  axis(2) = subplot(3,1,2);
  for k = 1:length(boutSpikes)
  plot(linspace(0,length(boutSpikes{:,k})/fs, length(boutSpikes{:,k})),rasterMat{:,k}-2*k, 'k', 'LineWidth', .5);
  set(gca,'xtick',[])       
  set(gca, 'ytick', [])
  hold on  
end

%make psth
bin = .04 * fs;
p = 1;
clear var psthMat
for m = 1:bin:length(boutSpikes{:,1})-bin %not all of it being sampled btw!!!!!!!!!!!
    sumAll = 0;
    for n = 1:numBouts
        cm = rasterMat{n};
        sumSingle = sum(cm(m:m+bin-1));
        sumAll = sumAll+sumSingle;
    end
    psthMat(p) = sumAll;
    p=p+1;
end

axis(3) = subplot(3,1,3);
plot(linspace(0,length(boutSpikes{:,1})/fs, length(psthMat)),psthMat, 'k', 'LineWidth', 2)
xlabel('Time (s)') 
ylabel('Frequency')
set(gca,'TickDir','out')    
ylim([-2 25])
%set(gca,'xtick',[])       
%set(gca, 'ytick', []) 
linkaxes(axis, 'x')

%% ISI comparisons
spont = spikes(1:starts(1)-1);

for i = 1:numBouts-1
    spont = [spont; spikes(stops(i)+1:starts(i+1)-1)];
end
spont = [spont; spikes(stops(end)+1:end)];
flipSpont = spont * -1;
 [valsSpont, locsSpont] = findpeaks(flipSpont, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
figure(2);
a(1) = subplot(2,1,1)
h1 = histogram(diff(locsSpont/fs),'FaceColor', 'r', 'EdgeColor', 'r');
h1.BinWidth = .020;
set(gca,'TickDir','out');
xlim([0 1])


pb = [];
for j = 1:numBouts
    pb = [pb; spikes(starts(j):stops(j))];
end
flipPb= pb * -1;
[valsPb, locsPb] = findpeaks(flipPb, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
a(2) = subplot(2,1,2)
h2 = histogram(diff(locsPb/fs), 'FaceColor', 'k');
h2.BinWidth = .020;
xlim([0 1])
set(gca,'TickDir','out');
shg

spontFR = length(locsSpont)/(length(spont)/fs)
pbFR = length(locsPb)/(length(pb)/fs)
spontMeanISI = mean(diff(locsSpont/fs))
pbMeanISI = mean(diff(locsPb/fs))
spontStdISI = std(diff(locsSpont/fs))
pbStdISI = std(diff(locsPb/fs))
pbCV = spontStdISI/spontMeanISI
pbCV = pbStdISI/pbMeanISI
