clear all; close all; clc;
baseFolder='\\rackstation\backup\Ellie\gcamp\dlx26_trainingrig\20190116\20190116_behavior\18start\';
fs=50000;
fileNames = dir([baseFolder '*.wav']);

for i = 1:length(fileNames)
    file = audioread([baseFolder fileNames(i).name]);
    figure;
    plot(file)
    title(i)
    shg
 end
