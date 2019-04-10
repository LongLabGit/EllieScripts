clear all; clc; close all;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\';
f = 'c1_1';
boutFolder = 'boutFinder\';
motifFolder = 'motifFinder\';
load([baseFolder boutFolder f '\boutAudio.mat']);
fs = 50000;
mkdir([baseFolder motifFolder f '\'])
cd([baseFolder motifFolder f '\'])
for i = 1:length(boutAudio)
    bout = boutAudio{:,i};
    fname = ['boutAudio_',  num2str(i), '.wav']
   % audiowrite(fname, bout, fs)
end

%% motifs only audio
motifInterval1 = [(60*fs):(150*fs)]; %in seconds
motifInterval2 = [(1000*fs):(1120*fs)];
fullPB=audioread([baseFolder boutFolder f '\' f '_playback.wav']);
motifPB1 = fullPB(motifInterval1);
motifPB2 = fullPB(motifInterval2);
audiowrite('motifsOnly1.wav', motifPB1, fs)
audiowrite('motifsOnly2.wav', motifPB1, fs)
save motifInterval1.mat motifInterval1
save motifInterval2.mat motifInterval2