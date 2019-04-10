%% open song file
clear all; close all; clc;
baseFolder='V:\Robert\INT_connectivity\SiProbe\BOS\surgery_040119\040119\040119\';
boutFile = audioread([baseFolder 'bout.wav']);
motifFile = audioread([baseFolder 'motif.wav']);
fs = 40000;
%% make playback file

figure;
vigiSpec(boutFile,fs)

a = 1 * fs; 
b = 5 * fs;
numMotifs = 20;
numSongs = 10;
itiJitter = randi([a b], 1, numMotifs + numSongs);
itiBoutStandard = zeros(30*fs, 1);
itiMotifStandard = zeros(10*fs, 1);
playbackFile = zeros(15*fs, 1);

for i = 1:numMotifs
    playbackFile = [playbackFile; motifFile; itiMotifStandard; zeros(itiJitter(i), 1)]; 
end 

playbackFile = [playbackFile; zeros(20*fs, 1)];

for i = 1:numSongs
    playbackFile = [playbackFile; boutFile; itiBoutStandard; zeros(itiJitter(i), 1)];
end
lengthFile = length(playbackFile)/fs/60

wavName = [baseFolder  'playbackFile.wav'];
audiowrite(wavName, playbackFile, fs);

figure;
plot(linspace(0,lengthFile, length(playbackFile)), playbackFile);
shg
