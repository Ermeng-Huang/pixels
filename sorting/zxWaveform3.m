clc
clear
addpath('D:\code\spikes-master\analysis') %code package downloaded from github cortex-lab/spikes
ddpath('D:\code\npy-matlab\npy-matlab\')
cd(lwd);
rootpath=lwd;

%% Extract waveform
if ~isfile([rootpath '\waveform.mat'])
cluster_ids=FindGoodCluster(rootpath);
waveform=cell(length(cluster_ids),4);
spikeTimes= readNPY('spike_times.npy');
spikeClusters =readNPY('spike_clusters.npy');

for i=1:length(cluster_ids)    
    
    WaveForms=getWaveForms(getParams(rootpath,cluster_ids(i),spikeTimes,spikeClusters));
    waveform{i,1}=rootpath;
    waveform{i,2}=cluster_ids(i);    
    
    BestChan=FindBestChannel(squeeze(WaveForms.waveFormsMean));
    waveform{i,3}=squeeze(WaveForms.waveForms(:,:,BestChan(1:4),:));
    waveform{i,4}=mean(squeeze(WaveForms.waveForms(:,:,BestChan(1),:)),1);
end
save([rootpath '\waveform.mat'],'waveform');
else
    load([rootpath '\waveform.mat'])
end

%% Indentify good waveform and plot good waveform
goodWaveform=oneFolder(rootpath);
w=waveform(goodWaveform(:,2)==1,:);
d=ceil(sqrt(size(w,1)));
fh=figure('Color','w','Position',[100,100,800,800]);
for sidx=1:size(w,1)
    subplot(d,d,sidx)
    onetrace=squeeze(mean(w{sidx,3},1));    
    plot(mean(onetrace,1));
    set(gca(),'XTick',[],'YTick',[])
    xlim([0,92])
    box off
    title(w{sidx,2})
end

saveas(fh,fullfile(rootpath,'waveform.png'),'png')
close(fh);
h5create(fullfile(rootpath,'goodWaveform.hdf5'),'/goodWaveform',size(goodWaveform),'Datatype','double')
h5write(fullfile(rootpath,'goodWaveform.hdf5'),'/goodWaveform',goodWaveform)

%% Function
function GoodCluster=FindGoodCluster(rootpath)
s1s=30000;
FR_Th=1;
metaf=ls(fullfile(rootpath,'*.ap.meta'));
fh=fopen(fullfile(rootpath,metaf));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/s1s/2*FR_Th;
clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
waveformGood=strcmp(clusterInfo{:,4},'good');
freqGood=clusterInfo{:,10}>spkNThresh;
GoodCluster = table2array(clusterInfo(waveformGood & freqGood,1));
end

function gwfparams=getParams(rootpath,cluster_ids,spikeTimes,spikeClusters)
gwfparams.dataDir = rootpath;    % KiloSort/Phy output folder
gwfparams.fileName = ls(fullfile(rootpath,'*.ap.bin'));         % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-30 60];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 500;                    % Number of waveforms per unit to pull out

n=spikeClusters==cluster_ids;
gwfparams.spikeClusters=spikeClusters(n);
gwfparams.spikeTimes=spikeTimes(n); 
end
function BestChannel=FindBestChannel(waveform_temp)
for nchannel=1:size(waveform_temp,1)
    amplitudes(nchannel)=max(waveform_temp(nchannel,:))-min(waveform_temp(nchannel,:));
end
   [~,BestChannel]=sort(amplitudes,'descend');
end

function out=oneFolder(folder)
stats=[];% type, trough-peak, fwhm
showcaseCount=0;
if ~isfile(fullfile(folder,'waveform.mat'))
    disp('Missing waveform file');
    disp(folder);
    out=[];
    keyboard();
    return
end

fstr=load(fullfile(folder,'waveform.mat'));
wf_all=fstr.waveform;


% for one_wf=wf_all(:,4)'
for i=1:size(wf_all,1)
    wf=wf_all{i,4};
    %criteria 1
    [mmin,mi]=min(wf);
    [mmax,xi]=max(wf((mi+1):end));
    if mmax>-mmin
%         disp('type 1 bad wf')
        stats(end+1,:)=[-1,0,0];
        continue
    end
    %criteria 2
    [lc_pk,lc_pk_ts]=findpeaks(wf,'MinPeakProminence',-0.05*mmin);
    if numel(lc_pk)>=6
%         disp('type 2 bad wf')
        stats(end+1,:)=[-2,0,0];
        continue
    end
    %criteria 3
    if ~isempty(xi) && xi>3
        [lc_pk,lc_pk_ts]=findpeaks(wf((mi+1):(mi+xi-1)),'MinPeakProminence',-0.05*mmin);
    end
    if isempty(xi) || (xi<=3) || (~isempty(lc_pk))
%         disp('type 3 bad wf')
        stats(end+1,:)=[-3,0,0];
        continue
    end
        
    %criteria 4
    %trough_peak dist
    wf=spline(1:length(wf),wf,1:0.03:length(wf));
    scale=max(abs(wf));
    wf=wf./scale;
    [trough,troughTS]=min(wf);
    [lp,deltaTS]=max(wf((troughTS+1):end));%to late peak
    
    %fwhm
    lcross=find(wf<-0.5,1);
    rcross=find(wf(troughTS:end)>-0.5,1)+troughTS;
    %         figure()
    if numel(lcross)==1 && numel(rcross)==1
        fwhm=rcross-lcross;
%         stats(end+1,:)=[0,deltaTS,fwhm];
        showcaseCount=showcaseCount+1;
        
        if showcaseCount>200 && showcaseCount<=250
            hold on
            plot(wf)
            plot(lcross,-0.5,'ro')
            plot(rcross,-0.5,'ro')
            plot(troughTS,trough,'bo')
            plot(troughTS+deltaTS,lp,'bo')
            hold off
        end
    else
%         disp('type 4 bad wf')
        %             plot(wf)
        stats(end+1,:)=[-4,0,0];
        pause(0.2)
        continue
    end
    stats(end+1,:)=[1,0,0];

end
out=[cell2mat(wf_all(:,2)),stats(:,1)];

end
