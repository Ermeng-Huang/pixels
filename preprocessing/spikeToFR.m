%%% Calculate firing rate from spike trains
clear
clc
cd('D:\pixel-optogenetic')
addpath('D:\code\npy-matlab\npy-matlab')
addpath('D:\code\fieldtrip-20200320') 
ft_defaults
fl=dir('D:\pixel-optogenetic\DataSum\*\*\cluster_info.tsv');
binSize=0.5; % sec
for f=1:length(fl)
    rootpath=fl(f,1).folder;
    
    if exist(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'file')        
        continue
    end
    [SU_id,trials,FR_All]=plotOneTrack(rootpath,binSize);

    h5create(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'/FR_All',size(FR_All),'Datatype','double')
    h5write(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'/FR_All',FR_All)
    h5create(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'/Trials',size(trials),'Datatype','double')
    h5write(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'/Trials',trials)
    h5create(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'/SU_id',size(SU_id),'Datatype','double')
    h5write(fullfile(rootpath,sprintf('FR_All_%dms.hdf5',binSize*1000)),'/SU_id',SU_id)
end

function [SU_id,trials,FR_All]=plotOneTrack(rootpath,binsize)
SU_id=plotOneDir(rootpath);

[FT_SPIKE,trials]=plot_su(rootpath,SU_id);
FR_All=plotOne(FT_SPIKE,binsize);
%     FR_All(i,:,:)=FR;

end

function [cluster_ids]=plotOneDir(rootpath)
sps=30000;
FR_Th=1.0;
metaf=ls(fullfile(rootpath,'*.meta'));
fh=fopen(fullfile(rootpath,metaf));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/sps/2*FR_Th;
clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
waveformGood=strcmp(clusterInfo{:,4},'good');
freqGood=clusterInfo{:,10}>spkNThresh;
cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
end

function [FT_SPIKE,trials]=plot_su(rootpath,SU_id)
sps=30000;
timeLim=[-3,11];
folder=regexp(rootpath,'(?<=DataSum\\)(.*)','match','once');
date=str2double(folder(5:12));
trials=markLPerf(h5read(fullfile(rootpath,'events.hdf5'),'/trials')');

% 20200730 laser-negative trials
if date==20200810
    trials(320:end,:)=[];
elseif date==20200812
    trials(241:end,:)=[];
elseif date>20200729     
    if ~isempty(find(diff(trials(:,3))>300, 1))        
        trials(find(diff(trials(:,3))>300)+1:end,:)=[];
    end
end

spkTS=h5read(fullfile(rootpath,'spike_times.hdf5'),'/spkTS');
spkId=readNPY(fullfile(rootpath,'spike_clusters.npy'));
FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(SU_id)));
FT_SPIKE.timestamp=cell(1,numel(SU_id));
for i=1:numel(SU_id)
    FT_SPIKE.timestamp{i}=spkTS(spkId==SU_id(i))';
end

cfg=struct();
cfg.trl=[trials(:,1)+timeLim(1)*sps,trials(:,1)+timeLim(2)*sps,zeros(size(trials,1),1)+timeLim(1)*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;
    
FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

end
function pfHist=plotOne(FT_SPIKE,binSize)

for j=1:size(FT_SPIKE.label)
    for i=1:size(FT_SPIKE.trialtime,1)        
        pfHist(i,:,j)=histcounts(FT_SPIKE.time{1,j}(1,FT_SPIKE.trial{1,j}==i),FT_SPIKE.trialtime(1,1):binSize:FT_SPIKE.trialtime(1,2))./binSize; 
    end
end
end
function [out]=markLPerf(facSeq)
% facSeq(:,9)=0;
% i=40;
% while i<=length(facSeq)
%     goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
%     if goodOff>=30 %.75 correct rate
%         facSeq(i-39:i,9)=1;
%     end
%     i=i+1;
% end
out=[facSeq,xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0)];
end
