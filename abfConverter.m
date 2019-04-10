%% 
clear all; close all; clc;

fs = 50000;
baseFolder = '\\rackstation\backup\Ellie\Int Juxta\032819_dlx32\';
cellName = 'c1_0004'
fileMat = abfload([baseFolder cellName '.abf']);
cellName(regexp(cellName, '0')) = [];
if cellName(end) == '_'
    cellName = [cellName '0'];
end

rawAudio = fileMat(:,3);
% rawAudio2 = fileMat(72700001:end,3);
cd([baseFolder 'boutFinder\']);
mkdir(cellName)
audio = [baseFolder 'boutFinder\' cellName '\' cellName '_playback.wav'];
audiowrite(audio, rawAudio, fs);
% audio = [baseFolder 'boutFinder\' cellName '\' cellName '_playback2.wav'];
% audiowrite(audio, rawAudio2, fs);

%% filteredspikes