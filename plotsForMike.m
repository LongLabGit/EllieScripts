clear all; close all; clc;
%file = abfload('V:\Ellie\Int Juxta\031419_dlx31\boutFinder\c2_6\c2_6_filteredspikes.abf'];
clear all; close all; clc;
fs = 50000;
hp = 200;
baseFolder = 'V:\Ellie\Int Juxta\032819_dlx32\';
cellName = 'c3_0';
files = load([baseFolder 'hp200\' cellName '_hp200.mat']);
file = files.hpSpikes;
%%
fs = 50000;
spikes = file(100*fs:115*fs,1); 
%axis(1) = subplot(2,1,1)
plot(linspace(0, length(spikes)/fs, length(spikes)),spikes)
shg;
%axis(2) = subplot(2,1,2)
%plot(linspace(0, length(spikes)/fs, length(spikes)), file(:,3))
%linkaxes(axis,'x');

%%
% figure(1)
% s1 = spikes(50*fs:60*fs-1);
% plot(s1)
% shg

% s1 = spikes(50*fs:60*fs-1);
% s2 = spikes(70*fs:80*fs-1);
% s3 = spikes(135*fs:145*fs-1);
% s4 = spikes(182*fs:192*fs-1);
% s5 = spikes(208*fs:218*fs-1);
