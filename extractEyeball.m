clear all; close all; clc;
data.name = 'cell2_6'; 

F='V:\Ellie\Int Juxta\031419_dlx31\eyeball\';

data.imgSampRate=2;%filming sampling rate
plotit=1;%optional parameter to look at it
files = [];
imgs=dir([F data.name '*tif']);
tIcreated=[imgs.datenum];%this will be useful later when it isnt always buffering 
imgN={imgs.name}';

%sort imgs based on 
number=nan(size(imgN));%get index from name
for ii=1:length(imgN)
    dashes=strfind(imgN{ii},'-');
    number(ii)=str2double(imgN{ii} (dashes(end)+1:end-4));
end
[~,inds]=sort(number);
imgN=imgN(inds);
tIcreated=tIcreated(inds);%check here that tI created is linear

for k = 1:length(imgN)
    files{k} = strcat(imgN{k});
end

%%
data.eyeball=nan(1,length(files));%initialize score

%open(v);
f = 1186;
    %load image
    G = [F,files{f}];
    img = imread(G);
    mg = double(img(:,:,1)); 
    
    %draw rois
    figure(1); 
    imagesc(img,[30,100])
    colormap gray
    axis equal
    axis tight
    disp(['Draw ROI']);
    roi = roipoly; %draw your roi for right eye

%% Run Parallel Pools for ROI extraction and crop
parfor f=1:length(imgN) %for each file
    %load image
    G = [F,files{f}];
    img = imread(G);
    img = double(img(:,:,1)); 
    imgRoi=img.*double(roi);% extract pixel intensity for roi by multiplying mask by image
    eyeball(f) = mean(imgRoi(:));% extract mean luminance value for entire roi 
end

data.eyeball = eyeball;

%%  eyeball score
final = data.eyeball-mean(data.eyeball);
data.final = final/max(final);
     
figure(100); 
plot(data.final);
title(['Eye score for ' (data.name)]); box off; set(gca,'TickDir','out'); 
ylabel('Score (a.u.)'); xlabel('Frames'); 

%%
figure(100); 
[x,~]=ginput(1);
indI=round(x);
figure; 
clickedImg=imread([F,files{indI}]);
imagesc(clickedImg, [0 255]);
colormap gray; 
numMontage = 5;
imgRange = floor(numMontage/2);
series = {};
j = 1;
for i = (indI-imgRange):(indI+imgRange)
    series{j} = [F, files{i}];
    j = j+1;
end
figure;
montage(series);
shg;

%% threshold eye score
thresh = -0.5; % threshold to separate open from closed conditions
In = data.final;

imgEpochs = [1 length(In)];

data.imgEpochs = imgEpochs;

data.closed = In > thresh; % value greater than threshold = eye closed
data.open = In < thresh; % value less than threshold = eye open
data.sumOpen = sum(data.open);
data.perClosed = data.sumOpen/length(data.open);

% for i = 1:length(data.open)
%     if data.open(i) == 1
%         
% end