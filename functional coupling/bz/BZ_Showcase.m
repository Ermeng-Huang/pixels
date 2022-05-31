%%%% 
clear
clc

if ispc
    CodePath='D:\code-hem';
    HomePath='F:\pixel-dualtask';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end
addpath(fullfile('D:\code','npy-matlab','npy-matlab'))
addpath(fullfile('D:\code','fieldtrip-20200320'))
addpath(CodePath)
ft_defaults
%% quick and dirty
id=unique(floor(h5read(fullfile(HomePath,'Selectivity_1108.hdf5'),'/cluster_id')/100000));
id=12
set(0,'DefaultFigureVisible', 'on') 
for fidx=id' 
    disp(fidx)
    BZ_SC_session(HomePath,fidx)
end

%% single pair
cell1=4000409;
cell2=4020494;
BZ_SC(HomePath,cell1,cell2)


%% function
function BZ_SC_session(homedir,fidx)
reg=regexp(h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/reg'),'(\w|\\|-)*','match','once');
cluster_id=h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/cluster_id');
path=regexp(h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/path'),'(\w|\\|-)*','match','once');
load(fullfile(homedir,'xcorr','bz0902',sprintf('BZ_XCORR_duo_f%d.mat',fidx)),'folder','mono')
reg_OBM1=isOBM1Map(reg(:,7));
conn=arrayfun(@(x)mono.completeIndex(mono.completeIndex(:,3)==x,2),mono.sig_con);

[avail,FT_SPIKE,trials]=pre_process(fullfile(homedir,'DataSum',folder),fidx,unique(conn));    
ts_id_tagged=cell(0);

for i=1:size(mono.sig_con)
    conn_reg=reg(ismember(cluster_id,conn(i,:)),7);
    if ~all(ismember(conn_reg,reg_OBM1(:,1)))
        continue
    end
    for trlIdx=1:size(trials,1)
        ts1=FT_SPIKE.time{unique(conn)==conn(i,1)}(FT_SPIKE.trial{unique(conn)==conn(i,1)}==trlIdx)';
        ts2=FT_SPIKE.time{unique(conn)==conn(i,2)}(FT_SPIKE.trial{unique(conn)==conn(i,2)}==trlIdx)';
        ts1(:,2)=1;
        ts2(:,2)=2;
        ts_id=[ts1;ts2];
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged{trlIdx}=fc_tag(ts_id); %%???find spike pairs (>2ms, <10ms)
        clear ts1 ts2
    end
    trial_good=find(cellfun(@(x)any(x(:,3)==1),ts_id_tagged));
    if nnz(trial_good)<3
        continue
    end    
    CCH=squeeze(mono.ccgR(:,mono.sig_con(i,1),mono.sig_con(i,2)));   
    [M,I]=findpeaks(CCH,'MinPeakHeight',median(CCH)*1.5,'SortStr','descend','NPeaks',2);
    if ~isempty(I)&&I(1)<(251+25)&&I(1)>(251-25)
        if length(M)==1||(length(M)>1&&M(1)>M(2)*1.5)
            hold on
            plot(1:501,CCH,'r-');
            arrayfun(@(x)plot(repmat(x,1,2),[min(CCH),M(1)],'k--'),[251-25,251,251+25])
            title(sprintf('%d-%d,%s-%s',conn(i,1),conn(i,2),conn_reg{1},conn_reg{2}))
            saveas(gcf,fullfile(homedir,'xcorr','showcase',sprintf('%d-%d_%s-%s.png',conn(i,1),conn(i,2),conn_reg{1},conn_reg{2})))
            close
        end
    end
    
end
end

function BZ_SC(homedir,cell1,cell2)
reg=regexp(h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/reg'),'(\w|\\|-)*','match','once');
cluster_id=h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/cluster_id');
path=regexp(h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/path'),'(\w|\\|-)*','match','once');

load(fullfile(homedir,'xcorr','bz0902',sprintf('BZ_XCORR_duo_f%d.mat',floor(cell1/100000))),'folder','mono')

fh=figure('Color','w','Position',[50,50,250,200]);
sgtitle(sprintf('%d-%d,%s-%s',cell1,cell2,reg{cluster_id==cell1,7},reg{cluster_id==cell2,7}))

%% spike event raster
[avail,FT_SPIKE,trials]=pre_process(fullfile(homedir,'DataSum',folder),floor(cell1/100000),[cell1,cell2]);    
ts_id_tagged=cell(0);
for trlIdx=1:size(trials,1)
    ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
    ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';
    ts1(:,2)=1;
    ts2(:,2)=2;
    ts_id=[ts1;ts2];    
    [~,s]=sort(ts_id(:,1));
    ts_id=ts_id(s,:);
    ts_id_tagged{trlIdx}=fc_tag(ts_id); %%???find spike pairs (>2ms, <10ms)    
    clear ts1 ts2
end
trial_good=find(cellfun(@(x)any(x(:,3)==1),ts_id_tagged));  

subplot(2,2,[1,2])
num=5; % trial plot num
b=125;
hold on
% cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'EdgeColor','none'),{[0,1]})
% cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'EdgeColor','none'),{[3,4],[6,7]})
alpha(0.2)
for itr =1:num % trial
    spike1=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==1,1)';
    spike2=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==2,1)';
    spike1_fc=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==1&ts_id_tagged{itr+b}(:,3)==1,1)';
    spike2_fc=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==2&ts_id_tagged{itr+b}(:,3)==1,1)';
    if isempty(spike1)
        spike1=-3;
    end
    if isempty(spike2)
        spike2=-3;
    end
    if isempty(spike1_fc)
        spike1_fc=-3;
    end
    if isempty(spike2_fc)
        spike2_fc=-3;
    end
    plot([spike1;spike1], [itr*2-1-0.25 itr*2-1+0.25],'color',[1,0.6,0.6]);
    plot([spike2;spike2], [itr*2-0.25 itr*2+0.25],'color',[0.6,0.6,1]);
    alpha(0.2)
    plot([spike1_fc;spike1_fc], [itr*2-1-0.25 itr*2-1+0.25],'r');
    plot([spike2_fc;spike2_fc], [itr*2-0.25 itr*2+0.25],'b');
end
%% CCH
CCH=squeeze(mono.ccgR(:,mono.completeIndex(mono.completeIndex(:,2)==cell1,3),mono.completeIndex(mono.completeIndex(:,2)==cell2,3)));
subplot(2,2,3)
hold on
plot(1:501,smooth(CCH),'r-');
arrayfun(@(x)plot(repmat(x,1,2),[600,1500],'k--'),[251-25,251,251+25])
set()


load(fullfile(replace(HomePath,'D','F'),'xcorr','bz0313',sprintf('BZ_XCORR_duo_f%d.mat',fidx)))

end

function [avail,out,trials]=pre_process(path,sessionIdx,cluster_ids)
sps=30000;
time=[1,2];

folderIdx=unique(floor((cluster_ids-100000*sessionIdx)/10000));
spkTS=[];
spkId=[];

for f=1:size(folderIdx,1)  
    folder=dir(fullfile(path,sprintf('*imec%d*',folderIdx(f)),'spike_clusters.npy'));
    spkId_temp=[];
    spkId_temp=readNPY(fullfile(folder.folder,'spike_clusters.npy'))+sessionIdx*100000;
    spkId=[spkId;spkId_temp+10000*str2num(regexp(folder.folder,'(?<=imec)(\d)','match','once'))];
    
    spkTS_temp=[];
    spkTS_temp=h5read(fullfile(folder.folder,'spike_times.hdf5'),'/spkTS');
    spkTS=[spkTS;spkTS_temp];
end


FT_SPIKE=struct();
FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
end
    
%  continuous format F T struct file
trials=markLPerf(h5read(fullfile(folder.folder,'FR_All_250ms.hdf5'),'/Trials'));
if isempty(trials)
    avail=false;
    out=[];
    return
end

cfg=struct();
cfg.trl=[trials(:,1)+time(1)*sps,trials(:,1)+time(2)*sps,zeros(size(trials,1),1)+time(1)*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end

function [out]=markLPerf(facSeq)
trial_min=120;
lick_criteria=0.7;
Perf_criteria=0.65;
i=trial_min;
facSeq_WT=zeros(length(facSeq),1);
while i<=length(facSeq)
    goodOff=nnz(xor(facSeq(i-trial_min+1:i,5)==facSeq(i-trial_min+1:i,6) , facSeq(i-trial_min+1:i,7)>0));
    lickOff=nnz(facSeq(i-trial_min+1:i,5)~=facSeq(i-trial_min+1:i,6)&facSeq(i-trial_min+1:i,7)>0);
    if goodOff>=Perf_criteria*trial_min && lickOff>=0.5*lick_criteria*trial_min %.60 correct rate
        facSeq_WT(i-trial_min+1:i,1)=1;
    end
    i=i+1;
end
out=[facSeq,facSeq_WT];

end


function out=isOBM1Map(reg)
load('D:/code-hem/OBM1Map.mat','OBM1map')

out=[];
for i=1:length(reg)
    try
        out{end+1,2}=OBM1map(reg{i});
        out{end,1}=reg{i};
    catch     
    end
end
end

function out=fc_tag(in)

out=in;
out(:,3)=0;

for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.01 && in(j,1)>in(i,1)+0.002
                out(i,3)=1;
                out(j,3)=1;
            end
            j=j+1;
        end
    end
end
end