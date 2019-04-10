%reboot;
addpath Fcns
F='V:\Vigi\Datasets\SiliconProbe\SingingandListening\11_9_16\Singing\';
load([F,'batches\Quant.mat']);
load([F,'\behavior\Behavior.mat']);
stability=cell(1,length(spk_t));
stability(ismember(clust,stable))={'Stable'};
stability(~ismember(clust,stable))={'Unstable'};
stability(strcmp(clust_type,'mua'))={'MU'};
set(0,'DefaultFigureWindowStyle','docked')
%% Set up some things
bw=.001;sm=20;%bw is resolution, sm is smoothign temporal window
%singing
SsongEdges=[-.28,.21];
StraceEdges=[-.5,.5];
tSing=[B(1).song{:}];
[audio,fs]=audioread([F,'behavior\audio.wav']);
sInds=round(((SsongEdges(1):1/fs:SsongEdges(2))+tSing(5))*fs);
singEx=audio(sInds);
Sbins=StraceEdges(1):bw:StraceEdges(2);
StHist=Sbins(1:end-1)+bw/2;
%listening
[playback,fs_s]=audioread([F,'\behavior\stim.wav']);
tListen=B(3).song;

% LsongEdges=[0,length(playback)/fs_s];
% LtraceEdges=[-.5,LsongEdges(2)+.5];
LsongEdges=SsongEdges;
LtraceEdges=StraceEdges;
Lbins=LtraceEdges(1):bw:LtraceEdges(2);
LtHist=Lbins(1:end-1)+bw/2;
lInds=round(((LsongEdges(1):1/fs:LsongEdges(2))+tListen(3))*fs);
listenEx=audio(lInds);
clear audio;
%%
clear h hh
% cells=find(ismember(clust,stable)&strcmp(clust_type,'good'))';
% cells=find(ismember(clust,[131 143 169 204 232 236 355 358]));
% cells=find(ismember(clust,[358 204 232 236 131 169 143 355]));
%cells=1:length(spk_t);
left = [1:5];
right = [1:5];

cells = find(ismember(clust,[266]));
for c=1:length(cells)
    figure(1);clf;   
    tSpike=spk_t{cells(c)};
    
    %First do Singing
    h(1)=subplot(5,2,1);
    plot(linspace(SsongEdges(1),SsongEdges(2),length(singEx)),singEx, 'r')
    %vigiSpec(singEx, fs)
    set(gca,'xtick',[])       
    set(gca, 'ytick', [])

    %Spike Raster
    h(2)=subplot(5,2,3:2:8);
        title('Singing')
        [times,trial]=timeLock(tSpike/2e4,tSing(left),StraceEdges);%wrap it around
    plotRaster(times,trial,length(tSing(left)),SsongEdges,'	[0, 0.4470, 0.7410]') %plot it
    ylabel('Trials')
   % xlabel('time (s)')
    set(gca,'xtick',[])       
    set(gca, 'ytick', [])
  %  linkaxes(h,'x')
%     h(3)=subplot(529);
%     NS = histcounts(times,Sbins);
%     if sm>1
%         NS=conv(NS,gausswin(sm)/sum(gausswin(sm)),'same');
%     end
%     NS=NS/bw/length(tSing);
%     plot(StHist,NS,'b');
%     line(SsongEdges(1)*[1,1],ylim,'color','k','linestyle','-','linewidth',2)
%     line(SsongEdges(2)*[1,1],ylim,'color','k','linestyle','-','linewidth',2)
%     xlabel('time (s)')
%     ylabel('Hz')


    
    
    %Then do Listening
    h(3)=subplot(5,2,2);
    plot(linspace(LsongEdges(1),LsongEdges(2),length(listenEx)),listenEx,'r')
    set(gca,'xtick',[])       
    set(gca, 'ytick', [])

    %Spike Raster
    h(4)=subplot(5,2,4:2:8);
    title('Listening')
vigi    [times,trial]=timeLock(tSpike/2e4,tListen(right),LtraceEdges);%wrap it around
    plotRaster(times,trial,length(tListen(right)),LsongEdges,'[0, 0.4470, 0.7410]')%plot it
    ylabel('Trials')
   % xlabel('time (s)')
    set(gca,'xtick',[])       
    set(gca, 'ytick', [])
    
%     hh(3)=subplot(5,2,10);hold on;
%     NL = histcounts(times,Lbins);
%     if sm>1
%         NL=conv(NL,gausswin(sm)/sum(gausswin(sm)),'same');
%     end
%     NL=NL/bw/length(tListen);
%     plot(LtHist,NL,'r');
%     plot(StHist,NS/max(NS)*max(NL),'b');
%     line(LsongEdges(1)*[1,1],ylim,'color','k','linestyle','-','linewidth',2)
%     line(LsongEdges(2)*[1,1],ylim,'color','k','linestyle','-','linewidth',2)
%     xlabel('time (s)')
%     ylabel('Hz')
     linkaxes(h,'x')
%     xlim(LtraceEdges)
%     axes(h(3));hold on;
%     plot(LtHist,NL/max(NL)*max(NS),'r');
    %Title it
%     suptitle(['Shank #' num2str(shank(cells(c))),', Cluster #' num2str(clust(cells(c))),...
%         ', Avg Spike Rate: ' num2str(round(length(tSpike)/range(tSpike)*2e4)) 'Hz, ' stability{cells(c)}])
    fr=num2str(round(length(tSpike)/range(tSpike)*2e4));
%     inds=LtHist>LsongEdges(1)&LtHist<LsongEdges(2);
%     r=corr(NS(inds)',NL(inds)');
%    suptitle(['Cluster #' num2str(clust(cells(c))), ', FR: ' fr 'Hz' ' Type: ' clust_type{c} ', ' stability{c} ])
%     drawnow;
    orient landscape
 %   export_fig(['SingingClusters\' num2str(c) '_Ellie.pdf'])
end

