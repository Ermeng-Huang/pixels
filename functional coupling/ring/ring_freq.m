clc;
clear;
homedir='D:\pixel-optogenetic\';
%% ring freq
load(fullfile(homedir,'xcorr','rings_bz.mat'))  
load(fullfile(homedir,'session_list.mat'))
for sess=1:size(session,1)
    [spkTS,spkID]=spk(session{sess},homedir,sess);
%     load(fullfile(homedir,'DataSum',session{sess},'FT_SPIKE.mat'))    
    for rsize=3        
        sums=cell(0);
        for ring_id=1:size(rings{sess,rsize-2},1)
            ts_id=[];ts_id_downsample=[];
            disp(ring_id);
            %     disp('RRRRR');
            cids=rings{sess,rsize-2}(ring_id,:);
            per_cid_spk_cnt=cids;
            for in_ring_pos=1:rsize               
                per_cid_spk_cnt(in_ring_pos)=nnz(spkID==cids(in_ring_pos));
                ts_id=cat(1,ts_id,[spkTS(spkID==cids(in_ring_pos))...
                    ,ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
%                 per_cid_spk_cnt(in_ring_pos)=nnz(FT_SPIKE.timestamp{str2num(cell2mat(FT_SPIKE.label))==cids(in_ring_pos)});
%                 ts_id=cat(1,ts_id,[FT_SPIKE.timestamp{str2num(cell2mat(FT_SPIKE.label))==cids(in_ring_pos)}'...
%                     ,ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
            end
            ts_id=sortrows(ts_id,1); 
            ts_id_downsample=ts_id(cell2mat(arrayfun(@(x)1000*(x-1)+1:1000*(x-1)+100,1:floor(size(ts_id,1)/1000)...
                ,'UniformOutput',false)),:);            
            ring_stats_downsample=relax_tag(ts_id_downsample,rsize);
            coact_count_downsample=sum(ring_stats_downsample.spk_cnt); 
         
            if coact_count_downsample<20
               continue 
            end
            ring_stats=relax_tag(ts_id,rsize);
            coact_count=sum(ring_stats.spk_cnt);   
            sums={sess,ring_id,cids,per_cid_spk_cnt,ring_stats,coact_count,ts_id(end,1)*0.1/30000};            
        end
        if ~isempty(sums)
            save(fullfile(homedir,'xcorr','ring',sprintf('ring_stats_%d_%d.mat',rsize,sess)),'sums','-v7.3');
        end    
    end
end
return

%%
onefile=dir(fullfile(homedir,'xcorr','ring','ring_stats_*_*.mat'));
sums=[];
for f=onefile'
fstr=load(fullfile(f.folder,f.name));
sums=cat(1,sums,fstr.sums);
[spkTS,spkID]=spk(session{sess},homedir,sess);
    load(fullfile(homedir,'DataSum',session{sess},'FT_SPIKE.mat'))  
end

cell2mat(sums(1:13,6))./cellfun(@sum,sums(1:13,4))
%% function
function [spkTS,spkID]=spk(Sessionfolder,HomePath,sessionIdx)
addpath('D:/code/npy-matlab/npy-matlab')
folder=dir(fullfile(HomePath,'DataSum',Sessionfolder,'*','FR_All.hdf5'));
cluster_ids=[];
for f=1:size(folder,1)
    cluster_ids_temp=[];
    cluster_ids_temp=h5read(fullfile(folder(f,1).folder,folder(1).name),'/SU_id')+sessionIdx*100000;
    cluster_ids=[cluster_ids;cluster_ids_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];   
end

spkTS=[];
spkID=[];
for f=1:size(folder,1)
    spkId_temp=[];
    spkId_temp=double(readNPY(fullfile(folder(f,1).folder,'spike_clusters.npy'))+sessionIdx*100000);
    spkID=[spkID;spkId_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
    
    spkTS_temp=[];
    spkTS_temp=double(h5read(fullfile(folder(f,1).folder,'spike_times.hdf5'),'/spkTS'));
    n=folder(f,1).folder;
    spkTS=[spkTS;spkTS_temp];
end
spk_bad=~ismember(spkID,cluster_ids);
spkTS(spk_bad)=[];
spkID(spk_bad)=[];
end
function out=relax_tag(in,rsize)
% called from \jpsth\+bz\+rings\rings_freq.m

arguments
    in (:,2) double
    rsize (1,1) double {mustBeMember(rsize,3:5)}
end

tags=zeros(size(in,1),1);
ring_idx=1;
curr_pre_ptr=1;
% curr_post_ptr=-1;
curr_ring=[];
starts=[];
ends=[];
spk_cnt=[];
durs=[];
% fprintf('000000');
while curr_pre_ptr<size(in,1)
%     if rem(curr_pre_ptr,100)==0, fprintf('%06d.',curr_pre_ptr);end

    %matching time window, assuming 30kHz
    syn_win_ubound=find(in((curr_pre_ptr+1):end,1)>in(curr_pre_ptr,1)+300,1); %%%%% 0.01s
    if isempty(syn_win_ubound), break;end %TODO use max available instead
    syn_win_lbound=find(in((curr_pre_ptr+1):end,1)>in(curr_pre_ptr,1)+24,1); 
    if isempty(syn_win_lbound), break;end %matching time window, assuming 30kHz
    
    %matching post spike
    cyc_post_pos=rem(in(curr_pre_ptr,2)+1,rsize);
    if cyc_post_pos==0, cyc_post_pos=rsize;end
    diff_post_ptr=find(in(curr_pre_ptr+1:end,2)==cyc_post_pos,1); %post unit
    
    if ~isempty(diff_post_ptr) && diff_post_ptr>=syn_win_lbound && diff_post_ptr<syn_win_ubound
        %TODO temp list ring spk
        if isempty(curr_ring), curr_ring=curr_pre_ptr;end
        curr_pre_ptr=curr_pre_ptr+diff_post_ptr;
        curr_ring=vertcat(curr_ring,curr_pre_ptr);
    else
        if numel(curr_ring)>rsize
            tags(curr_ring)=ring_idx;
            ring_idx=ring_idx+1;
            starts(end+1)=in(curr_ring(1),2);
            ends(end+1)=in(curr_ring(end),2);
            spk_cnt(end+1)=numel(curr_ring);
            durs(end+1)=in(curr_ring(end),1)-in(curr_ring(1),1);
        end
        curr_ring=[];
        curr_pre_ptr=curr_pre_ptr+1;
%         continue
    end
end

out=struct();
out.tags=tags;
out.starts=starts;
out.ends=ends;
out.spk_cnt=spk_cnt;
out.durs=durs;
end