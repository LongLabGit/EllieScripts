clear all; clc; close all;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\';
f = 'c1_1\';
boutFolder = 'boutFinder\';
motifFolder = 'motifFinder\';
audio = audioread([baseFolder motifFolder 'motif.wav']);
load([baseFolder boutFolder f 'boutSpikes.mat']);
load([baseFolder motifFolder f 'MotifTimes.mat']);
fs = 50000;
wing = 0*fs;

figure;
axis(1) = subplot(2,1,1);
newAudio = [zeros(wing,1); audio; zeros(wing,1)];
plot(linspace(0, length(newAudio)/fs, length(newAudio)), newAudio, 'k')
%  newSpec = [zeros(wing,1); audio; zeros(wing,1)];
%  for i=1:length(newSpec)
%      if newSpec(i) == 0
%          newSpec(i) = .2;
%      end
%  end
%   vigiSpec(newAudio, fs, [], .75);
   set(gca,'xtick',[])       
  set(gca, 'ytick', [])
  %end
axis(2) = subplot(2,1,2);
% templateLength = Motif(1).stop(1) - Motif(1).start(1);
templateLength = length(newAudio)/fs;
window = .005; %spike width in sec
thresh = .1;
rasterLocs = [];
numPeak = [];
valPeaks = [];

countAll = 0;
  for i = 1:length(Motif)
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
    
      for j = 1:length(Motif(i).start)
            currentLength = Motif(i).stop(j) - Motif(i).start(j);
            warp = templateLength/currentLength;
            motifStart = round(Motif(i).start(j)*fs);
            motifStop = round(Motif(i).stop(j)*fs);
%             tmp = boutSpikes{:,i};            
%             motifSpikes = tmp(motifStart:motifStop);
%             plot(linspace(0,length(motifSpikes)/fs*warp, length(motifSpikes)),motifSpikes-2*countAll, 'k');
%             plot(linspace(0, templateLength, length(motifSpikes)),motifSpikes-2*countAll, 'k');
            
           tmp = rasterMat{:,i};
           motifSpikes = tmp(motifStart:motifStop);
            for k = 1:length(motifSpikes)
               if motifSpikes(k) == 1
                   plot([k/fs*warp, k/fs*warp], [0, 3]-2*countAll,'Color', 'k')
               end
            end
            hold on 
            countAll = countAll + 1;
      end
  end   
 set(gca,'TickDir','out')    
 set(gca, 'ytick', [])
 linkaxes(axis, 'x') 




