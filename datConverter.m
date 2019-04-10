clear all; close all; clc;
baseFolder='V:\Robert\INT_connectivity\SiProbe\BOS\surgery_040119\040119\040119\';
fs=40000;
fileNames = dir([baseFolder '*.dat']);
for i = 1:length(fileNames)
%for i = 4
    fname = [baseFolder fileNames(i).name];
    [audioFileRaw, samplingRate, dateandtime, label, props] = LoadEGUI_daq(fname, 1);
    audioFileRaw = audioFileRaw./max(abs(audioFileRaw));
    filename=fileNames(i).name;
    audiowrite([baseFolder filename(1:end-4) '.wav'],audioFileRaw,fs);
    figure;
    vigiSpec(audioFileRaw,samplingRate);
    title(i)
    shg
 end

% %%%to display good ones
% good = [8;10;36]';
% goodFileNames = dir([baseFolder '*.wav']);
% for j = 1:length(good)    
%     fname = [baseFolder fileNames(j).name];
%     goodfname = [baseFolder baseFolder dayRecord '_d000' num2str(good(j)) '*.wav' .name];
%     [audioFileRaw, samplingRate, dateandtime, label, props] = LoadEGUI_daq(goodfname, 1);
%     audioFileRaw = audioFileRaw/max(abs(audioFileRaw));
%     figure
%     title('Good Specs')
%     subplot(length(good), 1, j) 
%     vigiSpec(audioFileRaw,samplingRate);
%     title(good(j))
% end
% save wavFile.mat wav