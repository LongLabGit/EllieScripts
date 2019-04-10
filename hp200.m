clear all; close all; clc;
fs = 50000;
hp = 200;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\';
cellName = 'c3_0000';
fileMat = abfload([baseFolder cellName '.abf']);
newFolder = 'hp200\';
cd([baseFolder]);
mkdir(newFolder);

cellName(regexp(cellName, '0')) = [];
if cellName(end) == '_'
    cellName = [cellName '0'];
end
rawSpikes = fileMat(:,1);
hpSpikes = highpass(rawSpikes, hp, fs);
cd([baseFolder newFolder]);
hpName = [cellName '_hp200.mat'];
save(hpName, 'hpSpikes');

% figure;
% ax(1)=subplot(2,1,1);
% plot(rawSpikes);
% ax(2)=subplot(2,1,2);
% plot(hpSpikes);
% linkaxes(ax, 'x');