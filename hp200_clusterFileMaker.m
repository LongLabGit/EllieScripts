clear all; close all; clc;
fs = 50000;
hp = 200;
baseFolder = 'V:\Ellie\Int Juxta\031419_dlx31\';
cellName = 'c2_0006';
fileMat = abfload([baseFolder cellName '.abf']);
newFolder = 'hp200\';
cd([baseFolder]);
mkdir(newFolder);
cellName(regexp(cellName, '0')) = [];
if cellName(end) == '_'
    cellName = [cellName '0'];
end

cd([baseFolder newFolder])
mkdir(cellName)
rawSpikes = fileMat(:,1);
%data = highpass(rawSpikes, hp, fs);
data = rawSpikes;
cd([baseFolder newFolder cellName]);
hpName = [cellName '_hp200.mat'];
save(hpName, 'data');

% figure;
% ax(1)=subplot(2,1,1);
% plot(rawSpikes);
% ax(2)=subplot(2,1,2);
% plot(hpSpikes);
% linkaxes(ax, 'x');