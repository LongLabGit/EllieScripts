clear all; close all; clc;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\boutFinder\';
f = 'c1_3';
audioF = [baseFolder f '\' f '_playback.wav']; 
templateF = 'V:\Ellie\Int Juxta\032819_dlx32\boutFinder\bout.wav';
thresh = .8;
[starts,stops,centers,warps] = findMotifs(audioF, templateF, thresh);

Motif.file = audioF; 
Motif.start = starts;
Motif.stop = stops;
Motif.center = centers;
Motif.warp = warps;
Motif.thresh = thresh; 

cd([baseFolder f '\']);
save MotifTimes.mat Motif