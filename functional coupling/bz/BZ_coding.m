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
id=unique(floor(h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/cluster_id')/100000));
for fidx=1:length(id)
    disp(id(fidx))
    BZcoding(HomePath,id(fidx))
end
function BZcoding(HomePath,fidx)

% ts_sep=0.2;
ts_sep=0;
pre_thresh=0;
trl_thresh=0;
load(fullfile(HomePath,'session_list.mat'),'session')
% load(fullfile(HomePath,'xcorr','bzdata',sprintf('BZ_XCORR_duo_f%d.mat',fidx)),'folder','mono')
mono=load(fullfile(HomePath,'xcorr','bzdata',sprintf('BZ_XCORR_inh_f%d.mat',fidx)));
folder=session{fidx};

% load(fullfile(replace(HomePath,'D','F'),'xcorr','bz0313',sprintf('BZ_XCORR_duo_f%d.mat',fidx)))
[avail,FT_SPIKE,trials]=pre_process(fullfile(HomePath,'DataSum',folder),fidx,'dualtask');
if ~avail
    return
end
%%
sums0=cell(0);
sums1=cell(0);
for idx=1:size(mono.sig_con,1)    
    %dimord={ts1,ts2, fcevents, anti-causal fc, 50ms-lag fc} x trl(~240) x bin(-3:11)    
    stats=zeros(6,size(trials,1),17);
    for trlIdx=1:size(trials,1)
        ts1=FT_SPIKE.time{mono.sig_con(idx,1)}(FT_SPIKE.trial{mono.sig_con(idx,1)}==trlIdx)';
        ts2=FT_SPIKE.time{mono.sig_con(idx,2)}(FT_SPIKE.trial{mono.sig_con(idx,2)}==trlIdx)';
%         ts1=ts1([false;diff(ts1)>ts_sep]);
        if isempty(ts1)
            stats(1,trlIdx,:)=0;
            stats(2,trlIdx,:)=0;
%             stats(3,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)==1 & ts_id_tagged(:,2)==1,1),-3:11); 
%             stats(4,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,4)==1 & ts_id_tagged(:,2)==1,1),-3:11);
%             stats(5,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,5)==1 & ts_id_tagged(:,2)==1,1),-3:11); 
%             stats(6,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)~=1 & ts_id_tagged(:,2)==1,1),-3:11);
            continue
        end
        ts1(:,2)=1;
        ts2(:,2)=2;
        ts_id=[ts1;ts2];
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged=fc.fc_tag(ts_id,false); %%???find spike pairs (>2ms, <10ms)
        stats(1,trlIdx,:)=histcounts(ts1(:,1),-3:14); %cell 1 spikes/s
        stats(2,trlIdx,:)=histcounts(ts2(:,1),-3:14); %cell 2 spikes/s
        stats(3,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)==1 & ts_id_tagged(:,2)==1,1),-3:14);% cell 1 spikes  followed by cell 2 spikes (2ms~10ms)
        stats(4,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,4)==1 & ts_id_tagged(:,2)==1,1),-3:14);
        stats(5,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,5)==1 & ts_id_tagged(:,2)==1,1),-3:14);
        stats(6,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)~=1 & ts_id_tagged(:,2)==1,1),-3:14);  
        clear ts1 ts2
    end
    
    %%
    trial_type=["distractorNo-correct","distractorNoGo-correct","distractorGo-correct"...
        ,"distractorNo-error","distractorNoGo-error","distractorGo-error"];
    onepair=nan(17,26,length(trial_type)); %-3~11 21???
    fc_spk=cell(0);
    for i=1:size(trial_type,2)
        [sel_S1,sel_S2]=util.ExtractTrial(trials,'task','dualtask','trials',trial_type(i));
        fc_spk{end+1}=stats([3,6,1,2],sel_S1,:);
        fc_spk{end+1}=stats([3,6,1,2],sel_S2,:);
        for sbin=1:14
            pre_sel=stats(1,:,sbin)'>pre_thresh; %choose trials with spiking of cell 1 (pre-cell)
            
            s1_fc=stats(3,intersect(sel_S1,find(pre_sel)),sbin);
            s1_pre=stats(1,intersect(sel_S1,find(pre_sel)),sbin);
            s1_post=stats(2,intersect(sel_S1,find(pre_sel)),sbin);
            s1_nofc=stats(6,intersect(sel_S1,find(pre_sel)),sbin);
            
            s2_fc=stats(3,intersect(sel_S2,find(pre_sel)),sbin);
            s2_pre=stats(1,intersect(sel_S2,find(pre_sel)),sbin);
            s2_post=stats(2,intersect(sel_S2,find(pre_sel)),sbin);
            s2_nofc=stats(6,intersect(sel_S2,find(pre_sel)),sbin);
            
            if nnz([s1_pre,s2_pre])==0 || (isempty(s1_pre)||isempty(s2_pre))
                p_pre=1;
            else
                p_pre=ranksum(s1_pre,s2_pre);
            end
            
            if nnz([s1_post,s2_post])==0 || (isempty(s1_post)||isempty(s2_post))
                p_post=1;
            else
                p_post=ranksum(s1_post,s2_post);
            end
            
            if nnz([s1_fc,s2_fc])==0 || (isempty(s1_fc)||isempty(s2_fc))
                p_fc=1;
            else
                p_fc=ranksum(s1_fc,s2_fc);
            end
            
            if nnz([s1_nofc,s2_nofc])==0 || (isempty(s1_nofc)||isempty(s2_nofc))
                p_fc=1;
            else
                p_fc=ranksum(s1_nofc,s2_nofc);
            end
            
            onepair(sbin,:,i)=[mean(s1_fc),std(s1_fc),numel(s1_fc),...
                mean(s2_fc),std(s2_fc),numel(s2_fc),...
                mean(s1_pre),std(s1_pre),numel(s1_pre),...
                mean(s2_pre),std(s2_pre),numel(s2_pre),...
                mean(s1_post),std(s1_post),numel(s1_post),...
                mean(s2_post),std(s2_post),numel(s2_post),...
                p_fc,p_pre,p_post,mean(s1_nofc),mean(s2_nofc),std(s1_nofc),std(s2_nofc),p_fc];
            
        end
    end
    sums0(end+1,:)={idx,mono.inh_conn(idx,:),fc_spk};
    sums1(end+1,:)={idx,mono.sig_con(idx,:),mono.inh_conn(idx,:),onepair};
    
%     sums0(end+1,:)={idx,[mono.completeIndex(mono.completeIndex(:,3)==mono.sig_con(idx,1),2),...,
%         mono.completeIndex(mono.completeIndex(:,3)==mono.sig_con(idx,2),2)],fc_spk};
%     sums1(end+1,:)={idx,mono.sig_con(idx,:),[mono.completeIndex(mono.completeIndex(:,3)==mono.sig_con(idx,1),2),...,
%         mono.completeIndex(mono.completeIndex(:,3)==mono.sig_con(idx,2),2)],onepair};
end
save(fullfile(HomePath,'xcorr','coding_inh',sprintf('fc_decoding_f%d.mat',fidx)),'sums0','trials','folder')
save(fullfile(HomePath,'xcorr','coding_inh',sprintf('fc_coding_f%d.mat',fidx)),'sums1','folder','trials')
if isunix
    quit(0)
else
    return
end
end

function [avail,out,trials]=pre_process(path,sessionIdx,model)
sps=30000;
folder=dir(fullfile(path,'*','FR_All_250ms.hdf5'));

spkTS=[];
spkId=[];
cluster_ids=[];
for f=1:size(folder,1)
    cluster_ids_temp=[];
    cluster_ids_temp=h5read(fullfile(folder(f,1).folder,folder(1,1).name),'/SU_id')+sessionIdx*100000;    
    cluster_ids=[cluster_ids;cluster_ids_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];

    spkId_temp=[];
    spkId_temp=readNPY(fullfile(folder(f,1).folder,'spike_clusters.npy'))+sessionIdx*100000;
    spkId=[spkId;spkId_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
    
    spkTS_temp=[];
    spkTS_temp=h5read(fullfile(folder(f,1).folder,'spike_times.hdf5'),'/spkTS');
    spkTS=[spkTS;spkTS_temp];
end

FT_SPIKE=struct();
FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));

FT_SPIKE.timestamp=arrayfun(@(x)spkTS(spkId==x)',cluster_ids,'UniformOutput',false);   
%  continuous format F T struct file
trials=util.markLPerf(h5read(fullfile(folder(1,1).folder,folder(1,1).name),'/Trials'),0.7,0.8,120,'dualtask');
if isempty(trials)
    avail=false;
    out=[];
    return
end

cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+14*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end

function [out]=markLPerf(facSeq)
i=40;
facSeq_WT=zeros(length(facSeq),1);
while i<=length(facSeq)
    goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
    if goodOff>=24 %.60 correct rate
        facSeq_WT(i-39:i,1)=1;
    end
    i=i+1;
end
out=[facSeq,facSeq_WT];

end
