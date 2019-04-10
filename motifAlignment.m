%% load shit 
clear all; clc; close all; 
baseFolder = 'V:\Ellie\X Juxta\030119_greenxdex_redtag_window\';
cellNum = 'c1_';
files = [1; 2];
fs = 50000;
window = .002 * fs; %spike width in sec
boutFolder = 'boutFinder\';
motifFolder = 'motifFinder\';
%f = [cellNum + string(files))];
load([baseFolder + boutFolder + f + '\thresh.mat']);
audio = audioread([baseFolder + motifFolder + 'motif.wav']);
%% do the things
countAll = 0;
axis(1) = subplot(2,1,1);
wing = 0*fs;
newAudio = [zeros(wing,1); audio; zeros(wing,1)];
plot(linspace(0, length(newAudio)/fs, length(newAudio)), newAudio, 'k')
vigiSpec(newAudio,fs)
set(gca,'xtick',[])       
set(gca, 'ytick', [])

for x = 1:length(files)
    
    f = [cellNum + string(files(x)) + '\'];
    load(strcat(baseFolder, boutFolder, f, 'boutSpikes.mat'));
    load(strcat(baseFolder, motifFolder, f, 'MotifTimes.mat'));

    axis(2) = subplot(2,1,2);
    % templateLength = Motif(1).stop(1) - Motif(1).start(1);
    templateLength = length(newAudio)/fs;
    rasterLocs = [];
    numPeak = [];
    valPeaks = [];

      for i = 1:length(Motif)
    %  flipped = -1 * boutSpikes{:,i};
        flipped = boutSpikes{:,i};
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
                       plot([k/fs*warp, k/fs*warp], [0, 3]-2.5*countAll,'Color', 'k')
                   end
                end
                hold on 
                countAll = countAll + 1; 
          end
      end   
     set(gca,'TickDir','out')    
     set(gca, 'ytick', [])
     xlim([0 1])
     linkaxes(axis, 'x') 
     hold on
end
%% motifs only

