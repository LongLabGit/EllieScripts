         
%uses motif finder start/stop times to display ephys on playback trials

%% fix audio file first!!!!!
clear all; close all; clc;  warning off;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\boutFinder\';
f = 'c1_4';
load([baseFolder f '\' 'MotifTimes.mat']);
% load([baseFolder 'cutoff.mat'])
fullRecord = abfload([baseFolder f '\' f '_filteredspikes.abf']);
%fullRecord = abfload([baseFolder f '\' f '.abf']);
fs = 50000;
wing = 5 * fs;
thresh = .2;
audio = fullRecord(:,3);
%% load shit
spikes = fullRecord(:,1); 
bout = audioread([baseFolder 'bout.wav']);
%motif = audioread([baseFolder 'motif.wav']);
lengthRecording = length(spikes)/fs
starts = Motif.start * fs; %change this later
stops = Motif.stop * fs; %change this later
numBouts = length(stops);
boutSpikes = {numBouts};
lengthBouts = [];
boutAudio = {numBouts};
for i = 1:numBouts
    boutSpikes{:,i} = spikes(starts(i)-wing:stops(i)+wing);
    lengthBouts(:,i) = length(boutSpikes{:,i});
    boutAudio{:,i} = audio(starts(i)-wing:stops(i)+wing);
end

for i = 1:numBouts
        if length(boutSpikes{:,i}) > min(lengthBouts)
           vec = boutSpikes{:,i};
           difference = length(vec)-min(lengthBouts);
           fixedVec = vec(1:end-difference);
           boutSpikes{:,i} = fixedVec;
        end
end
%% plot aligned spikes to bout trials
% figure(20);
% ax(1) = subplot(2,1,1);
% %plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})), audio(starts:stops), 'k')
% plot(audio(starts(1)-wing:stops(1)+wing), 'k')
% set(gca,'xtick',[])       
% set(gca, 'ytick', [])
% hold on
% title(f)

% for i = 1:numBouts
%     ax(2) = subplot(2,1,2);
%    % plot(boutSpikes{:,i}, 'k')
%     plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})),boutSpikes{:,i}-3*i, 'k', 'LineWidth', .5);
%     set(gca,'TickDir','out')       
%   %  set(gca, 'ytick', [])
%     hold on 
% end

window = .002 * fs; %spike width in sec
rasterLocs = [];
numPeak = [];
valPeaks = [];

for i = 1:numBouts
    %flipped = boutSpikes{:,i};
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
 
% % make raster plot
%for k = 1:numBouts
  %  stem(cell2mat(rasterMat{:,k})-10, 'Marker','none');
 axis(1) = subplot(2,1,1);
 %  plot(linspace(0,length(audio(starts(1)-100000:stops(1)+100000))/fs, length(audio(starts(1)-100000:stops(1)+100000)),audio(starts(1)-100000:stops(1)+100000), 'k')
 plot(audio(starts(1)-wing:stops(1)+wing), 'k')
 lenAudio = audio(starts(1)-wing:stops(1)+wing);
 pad = length(lenAudio) - length(bout);
 if mod(pad, 2) == 1
     pad = pad-1;
 end
 newSpec = [zeros(pad/2,1); bout; zeros(pad/2,1)];
 
%  for i=1:length(newSpec)
%      if newSpec(i) == 0
%          newSpec(i) = .2;
%      end
%  end
 
  vigiSpec(newSpec, fs, [], .45);
  set(gca,'xtick',[])       
  set(gca, 'ytick', [])
  
  for i = 1:numBouts
    axis(2) = subplot(2,1,2);
    if mod(i,5) == 0
        plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})),boutSpikes{:,i}-1*i, 'r')
    else
        plot(linspace(0,length(boutSpikes{:,i})/fs, length(boutSpikes{:,i})),boutSpikes{:,i}-1*i, 'k')
    end 
    hold on
  end
set(gca,'TickDir','out')    
set(gca, 'ytick', [])

linkaxes(axis, 'x') 

% % %   axis(2) = subplot(3,1,2);
% % %   plot(linspace(0,length(boutSpikes{:,k})/fs, length(boutSpikes{:,k})),rasterMat{:,k}-2*k, 'k', 'LineWidth', .5);
% % %   set(gca,'xtick',[])       
% % %   set(gca, 'ytick', [])
% % %   hold on  
% % %   end
% % 
% % 
% % 
% %make psth
% bin = .04 * fs;
% p = 1;
% clear var psthMat
% for m = 1:bin:length(boutSpikes{:,1})-bin %not all of it being sampled btw!!!!!!!!!!!
%     sumAll = 0;
%     for n = 1:numBouts
%         cm = rasterMat{n};
%         sumSingle = sum(cm(m:m+bin-1));
%         sumAll = sumAll+sumSingle;
%     end
%     psthMat(p) = sumAll;
%     p=p+1;
% end
% 
% axis(3) = subplot(3,1,3);
% plot(linspace(0,length(boutSpikes{:,1})/fs, length(psthMat)),psthMat, 'k', 'LineWidth', 2)
% xlabel('Time (s)') 
% ylabel('Frequency')
% set(gca,'TickDir','out')    

cd([baseFolder f '\']);
save boutSpikes.mat boutSpikes
save newSpec.mat newSpec
save boutAudio.mat boutAudio
save wing.mat wing
save thresh.mat thresh
data = spikes;
save data.mat data
%%
% fig = figure; 
% %a(1)=subplot(2,1,1)
% %plot(audio)
% %hold on
% %for i=1:length(boutSpikes)
% %    plot(starts(i):stops(i), audio(starts(i):stops(i)), 'r')
% %    hold on
% %end
% %a(2)=subplot(2,1,2)
% plot(linspace(0,lengthRecording, length(spikes)), spikes)
% ylim([-1.5 1.5])
% hold on
% for i=1:length(boutSpikes)
%     %plot((starts(i):stops(i))/fs, (spikes(starts(i):stops(i)))/fs, 'r')
%  	plot((starts(i):stops(i))/fs, zeros(length(starts(i):stops(i)),1), 'r', 'LineWidth', 5)
%     hold on
% end
% title(f)
% set(gca,'TickDir','out')    
% ylabel('mV')
% xlabel('Time (s)')
% savefig(fig, f, 'compact')
% %linkaxes(a, 'x') 
% 
%% eyeball
% cd([baseFolder f '\']);
% eyeTrig = fullRecord(:,4);
% [valEye,locEye] = findpeaks(eyeTrig, 'minPeakHeight', 2);
% numEyeFrames = length(locEye);
% save numEyeFrames.mat numEyeFrames 