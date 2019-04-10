%set things up
clear all; close all; clc;
baseFolder = '\\rackstation\backup\Ellie\Int Juxta\032819_dlx32\eyeball\';
cellName = 'cell1_1';
thresh = -.08;
files = dir([baseFolder cellName '*.tif']);
fs = 2;
numImgs = length(files);

%get eye data
loaded = importdata([baseFolder 'Results_' cellName '.txt']);
roiLum = loaded.data(:,3);
%minLum = loaded.data(:,6);
maxLum = loaded.data(:,7);
roiVals = roiLum./(maxLum);
%roiVals = roiLum./(maxLum./minLum);
roiVals = roiVals - mean(roiVals);
deltaRoi = roiVals;

%threshold images
figure(1)
axis(1) = subplot(2,1,1);
plot(linspace(0, length(deltaRoi)/fs/60, length(deltaRoi)), deltaRoi);
hold on
line([0 numImgs/fs/60], [thresh thresh], 'Color', 'r', 'LineWidth', 3);
imgs = {files.name};
[locs, ~] = find(deltaRoi < thresh); 

%make sliding window raster plot of open frames
raster = zeros(numImgs,1);
for i = 1:length(locs)
  ind = locs(i);
  raster(ind) = 1;
end
axis(1) = subplot(2,1,2);
plot(linspace(0, length(raster)/fs/60, length(raster)), movmean(raster, 20)); shg;
xlabel('Time(min)');
linkaxes(axis,'x');

%make montage
openImgFiles = {};
for i = 1:length(locs)
    ind = locs(i);
    openImgFiles(i) = {[baseFolder char(imgs{ind})]};
end 
figure(2);
montage(openImgFiles);

%save variables
cd(baseFolder);
save([cellName '_eyelocs.mat'], locs)
