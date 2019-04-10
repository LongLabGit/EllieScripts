%% load shit
clear all; close all; clc; warning off;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\';
subFolder = 'boutFinder\';
f = 'c1_1';
fs = 50000;
smoothWindow = .002 * fs;
threshY = .7e-4;
windBack = .5 * fs;
windFor = .5 * fs;
spikedata = abfload([baseFolder subFolder f '\' f '_filteredspikes.abf']);
spikedata = spikedata(:,1);
spikes = spikedata((1):(1340*fs));
clear spikedata;    

%% extract mvmt envelopes
mvmtTime = 1:length(spikes)';
mvmtBout = spikes(mvmtTime);
xAxis = (linspace(0,length(mvmtBout)/fs, length(mvmtBout)))';
smoothed = smoothdata(mvmtBout, 'gaussian', smoothWindow);
smoothed = smoothed.^2;
smoothed = normalize(smoothed, 'range');
raster = (smoothed >= threshY);

%calculate mvmt env timepoints
env = movmean(smoothed, [windBack windFor]);
[locs, ~] = find(env > threshY);
startMove = [locs(1)];
stopMove = [];
for i = 1:(length(locs)-1)
    if locs(i+1)-locs(i) ~= 1
        stopMove = [stopMove; locs(i)];
        startMove = [startMove; locs(i+1)];
    end
end
stopMove = [stopMove; locs(end)];

%concatenate close together mvmt bouts
for i = 2:length(startMove)
    if (startMove(i)-stopMove(i-1)) < (.1 * fs)
        startMove(i) = nan;
        stopMove(i-1) = nan;
    end
end
startMove = startMove(~isnan(startMove));
stopMove = stopMove(~isnan(stopMove));

%remove short mvmt timestamps
for i = 1:length(startMove)
    if stopMove(i)-startMove(i) < .1 * fs;
        startMove(i) = nan;
        stopMove(i) = nan;
    end
end
startMove = startMove(~isnan(startMove));
stopMove = stopMove(~isnan(stopMove));

figure; 
ax(1) = subplot(3,1,1);
plot(xAxis, mvmtBout); 
%ylim([-.7 .7]);
hold on;
plot(startMove/fs, threshY, 'ro')
hold on;
plot(stopMove/fs, threshY, 'ko')
xlim([xAxis(1) xAxis(end)]);

ax(2) = subplot(3,1,2);
plot(xAxis, smoothed);
hold on;
plot(startMove/fs, threshY, 'ro')
hold on;
plot(stopMove/fs, threshY, 'ko')
ylim([-.002 .1]);
xlim([xAxis(1) xAxis(end)]);
% 
% ax(3) = subplot(4,1,3);
% plot(xAxis,raster);
% xlim([xAxis(1) xAxis(end)]);

ax(4) = subplot(3,1,3);
plot(xAxis, env);
xlim([xAxis(1) xAxis(end)]);
ylim([0 .1])
hold on;
plot(startMove/fs, threshY, 'ro')
hold on;
plot(stopMove/fs, threshY, 'ko')
hold on;
line([xAxis(1) xAxis(end)],[threshY threshY])
linkaxes(ax, 'x');
%% remove mvmt bouts from spikedata
fixedSpikes = spikes;
for i = 1:length(startMove)
    for j = startMove(i):stopMove(i)
        fixedSpikes(j) = 0;
    end
end
figure; plot(fixedSpikes); shg;
cd([baseFolder subFolder f])
save fixedSpikes.mat fixedSpikes