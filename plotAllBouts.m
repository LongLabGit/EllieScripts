
%% combine data matrices
clear all; clc; close; 

baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\boutFinder\';%change

cellNum = 'c1_'; %change
files = [1]; %change 
eye = 0; %is there eye data

allBouts = [];
eyeFrames = [];
for i = 1:length(files)
    f = [cellNum + string(files(i))];
    if i == 1
        load([baseFolder + f + '\newSpec.mat'])
        load([baseFolder + f + '\wing.mat'])
        load([baseFolder + f + '\thresh.mat'])
    end
    load(baseFolder + f + '\boutSpikes.mat');
    boutSpikes = cell2mat(boutSpikes);
    allBouts = [allBouts boutSpikes];   
    if eye == 1
        load([baseFolder + f + '\numEyeFrames.mat']);
        eyeFrames = [eyeFrames; numEyeFrames];
    end
end
 cd(baseFolder)
 save allBouts.mat allBouts
 

%% plot spikes and psth
fs = 50000;
window = .002 * fs; %spike width in sec
numBouts = size(allBouts, 2);
rasterLocs = [];
numPeaks = [];
valPeaks = [];

%subplot 1
figure;
axis(1) = subplot(3,1,1);
%plot(linspace(-wing/fs,wing/fs, length(newSpec)), newSpec, 'k')
vigiSpec(newSpec, fs, [], .5)
set(gca,'xtick',[])       
set(gca, 'ytick', [])

for i = 1:numBouts
    flipped = -1 * allBouts(:,i);
   % flipped = allBouts(:,i);
    %flipped = allBouts(:,i);
    [vals, locs] = findpeaks(flipped, 'MinPeakHeight', thresh, 'MinPeakDistance', window);
    rasterLocs = [rasterLocs; locs];
    numPeaks = [numPeaks; length(locs)];
    valPeaks = [valPeaks; vals];
    boutRaster = zeros(length(allBouts(:,i)),1);
    for j = 1:length(locs)
        boutRaster(locs(j)) = 1;
    end
    rasterMat{:,i} = boutRaster;
end

%subplot 2
for j = 1:numBouts
    axis(2) = subplot(3,1,2);
   % plot(linspace(0,length(allBouts)/fs, length(allBouts)),allBouts(:,j)-2*j, 'k', 'LineWidth', 2);
   
   for k = 1:length(rasterMat{:,j})
       logicals = rasterMat{:,j};
       if logicals(k) == 1
           plot([k/fs, k/fs], [0, 3]-3*j,'Color', 'k')
       end
   end
   hold on
end
set(gca,'xtick', [])       
set(gca, 'ytick', [])

% %subplot 3
% bin = .030 * fs;
% p = 1;
% for m = 1:bin:length(allBouts)-bin %not all of it being sampled btw!!!!!!!!!!!
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
% plot(linspace(0, length(allBouts)/fs, length(psthMat)),psthMat, 'k', 'LineWidth', 2)
% xlabel('Time (s)') 
% ylabel('Frequency')
% set(gca,'TickDir','out')    

linkaxes(axis, 'x') 