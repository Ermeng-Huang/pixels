%  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
%  chan2. Thus, a negative lag represents a spike in the second dimension of
%  stat.xcorr before the channel in the third dimension of stat.stat.

if ispc
    CodePath='D:\code';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
addpath(fullfile(CodePath,'npy-matlab','npy-matlab'))
addpath(fullfile(CodePath,'fieldtrip-20200320'))
ft_defaults
savepath=fullfile(HomePath,'xcorr','laseroff');
currmodel='full-all'; %neuron-trials
prefix='0325';
delay=6;
bin_range=[-3 14];

load(fullfile(HomePath,'session_list.mat'))
range=[69,69]


for i=range(1):range(2)
    tic    
    if isfile([savepath sprintf('%s_%s_XCORR_duo_f%d_%1.1f_%1.1f_2msbin.mat',prefix,currmodel,i,bin_range(1),bin_range(2))])
        continue
    end
    folder=session{i}; 
    disp(folder)
    [avail,spktrial]=pre_process(HomePath,i,currmodel,folder); % posix 
   if contains(currmodel, 'resample')
       for r=1:100
          [xc_s1(:,:,:,r),xcshuf_s1(:,:,:,r),xc_s2(:,:,:,r),xcshuf_x2(:,:,:,r)]=plotxcorr(spktrial,delay,bin_range);
       end
       sums={i,folder,mean(xc_s1,4),mean(xcshuf_s1,4),mean(xc_s2,4),mean(xcshuf_x2,4)}; %per folder save
   else
       if avail
           [xc_s1,xcshuf_s1,xc_s2,xcshuf_x2]=plotxcorr(spktrial,delay,bin_range);
       end
       sums={i,folder,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save
   end
    
%    sums={i,folder,sustIds,transIds,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save    
   
   save(fullfile(savepath,sprintf('%s_%s_XCORR_duo_f%d_%d_%d_2msbin.mat',prefix,currmodel,i,bin_range(1),bin_range(2))),'sums','-v7.3')
    toc
end


return


function [avail,out]=pre_process(homedir,sessionIdx,model,Sessionfolder)
sps=30000;
folder=dir(fullfile(homedir,'DataSum',Sessionfolder,'*','FR_All.hdf5'));

cluster_ids=[];
for f=1:size(folder,1)
    cluster_ids_temp=[];
    cluster_ids_temp=h5read(fullfile(folder(f,1).folder,'FR_All.hdf5'),'/SU_id')+sessionIdx*100000;    
    cluster_ids=[cluster_ids;cluster_ids_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
end

if exist(fullfile(homedir,'DataSum',Sessionfolder,'FT_SPIKE.mat'),'file')
    load(fullfile(homedir,'DataSum',Sessionfolder,'FT_SPIKE.mat'))
else
    spkTS=[];
    spkId=[];
    for f=1:size(folder,1)
        spkId_temp=[];
        spkId_temp=readNPY(fullfile(folder(f,1).folder,'spike_clusters.npy'))+sessionIdx*100000;
        spkId=[spkId;spkId_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
        
        spkTS_temp=[];
        spkTS_temp=h5read(fullfile(folder(f,1).folder,'spike_times.hdf5'),'/spkTS');
        n=folder(f,1).folder;
        spkTS=[spkTS;spkTS_temp];
    end

    FT_SPIKE=struct();
    FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
    FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
    for i=1:numel(cluster_ids)
        FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
    end
    % save(fullfile(homedir,'DataSum',Sessionfolder,'FT_SPIKE.mat'),'FT_SPIKE')
    
    %  continuous format F T struct file
    
    Trials=h5read(fullfile(folder(1,1).folder,'FR_All.hdf5'),'/Trials');
    
    trials=clearBadPerf(Trials,model);
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
	save(fullfile(homedir,'DataSum',Sessionfolder,'FT_SPIKE.mat'),'FT_SPIKE')
end
out=FT_SPIKE;
avail=true;


end



function out=clearBadPerf(facSeq, model)
if strcmp(model, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
elseif contains(model, 'laseron')
    if length(facSeq)>=40
        out=facSeq(facSeq(:,9)==2,:);
    else
        out=[];
    end
elseif contains(model, 'laseroff')
    if length(facSeq)>=40
        out=facSeq(facSeq(:,9)==-1,:);
    else
        out=[];
    end
elseif contains(model, 'laseroff+resample')
    if length(facSeq)>=40
        out0=facSeq(facSeq(:,9)==2,:);
        out1=facSeq(facSeq(:,9)==-1,:);
        
        out=datasample(out1,length(out0),'Replace',false);
        [~,index]=sort(out(:,1));
        out=out(index,:);
    else
        out=[];
    end
elseif strcmp(model, 'correct_resampled')    
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        errornum=size(facSeq(errorsel,:),1);
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
            if goodOff>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out0=facSeq(facSeq(:,9)==1,:);
        try
            out=out0(randsample(size(out0,1),errornum),:);
        catch
            out=datasample(out0,errornum);
        end
    else
        out=[];
    end
elseif contains(model, 'all')
    out=facSeq;
elseif strcmp(model, 'learning_correct')
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
            if goodOff>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(facSeq(:,9)==1,:);
    else
        out=[];
    end
else    
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
            if goodOff>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(facSeq(:,9)==1,:);
    else
        out=[];
    end
end
end


function [Xc_S1,Xshuff_S1,Xc_S2,Xshuff_S2]=plotxcorr(spikeTrials,delay,bin_range)
% https://www.nature.com/articles/nn799
% A role for inhibition in shaping the temporal flow of information in prefrontal cortex
% Christos Constantinidis, Graham V. Williams & Patricia S. Goldman-Rakic
% Nature Neuroscience volume 5, pages175-180(2002)
%
% Neuron, Volume 76
% Functional Microcircuit Recruited during Retrieval of Object Association Memory in Monkey Perirhinal Cortex
% Toshiyuki Hirabayashi, Daigo Takeuchi, Keita Tamura, and Yasushi Miyashita


cfg             = [];
cfg.maxlag      = 0.1; % maximum 100 ms
cfg.binsize     = 0.002; % bins of 2 ms
cfg.outputunit  = 'raw'; % make unit area
cfg.latency     = bin_range; % time bin based on sample onset
cfg.vartriallen = 'no'; % allow variable trial lengths
cfg.debias      = 'no';

cfg.trials      = find(spikeTrials.trialinfo(:,5)==4 & spikeTrials.trialinfo(:,8)==delay);
if numel(cfg.trials)<2
    Xc_S1=[];
    Xshuff_S1=[];
else
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc_S1 = ft_spike_xcorr(cfg,spikeTrials);
    cfg.method      = 'shiftpredictor'; % compute the shift predictor
    Xshuff_S1 = ft_spike_xcorr(cfg,spikeTrials);
end

cfg.trials      = find(spikeTrials.trialinfo(:,5)==8 & spikeTrials.trialinfo(:,8)==delay);
if numel(cfg.trials)<2
    Xc_S2=[];
    Xshuff_S2=[];
else
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc_S2 = ft_spike_xcorr(cfg,spikeTrials);
    cfg.method      = 'shiftpredictor'; % compute the shift predictor
    Xshuff_S2 = ft_spike_xcorr(cfg,spikeTrials);
end
% compute the shuffled correlogram
end








