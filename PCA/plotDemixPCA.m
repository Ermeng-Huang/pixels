% calculated plane for AIopto
close all
clc;
clear;
% homedir='D:\pixel-optogenetic\';
homedir='F:\pixel-dualtask\';
addpath('D:\code-hem\dPCA')
%% Dataset
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','sel','sel');
% load(fullfile(homedir,'difftype_ByPEV.mat'),'lost','gain')
% NeuronSet=structfun(@(x)x(ismember(NeuronSet.id,lost.md.id)),NeuronSet,'UniformOutput',false);

[FR,FR_mean,trialNum]=NeuronSetForPCA_dualtask(homedir,NeuronSet,'delay');
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T

%%
demixPCA(FR,FR_mean,trialNum,'sel_delay',homedir)
%% Main Function
function [out1,out2,trialNum]=NeuronSetForPCA_dualtask(homedir,NeuronSet,phase)
idx=1;
trialNum=[];avail=[];
trial_type=["distractorNo-correct","distractoro-correct"];
% trial_type=["distractorNo-correct","distractorGo-correct","distractorNoGo-correct"...
%     ,"distractorNo-error","distractorGo-error","distractorNoGo-error"];
switch(phase)
    case 'delay'
        time=17:48;
    case 'MD'
        time=21:36;
    case 'LD'
        time=41:56;        
end
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    switch(length(trial_type))
        case 2
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                trialNum(idx:idx+ncell-1,1,i)=nnz(sel_S1);
                trialNum(idx:idx+ncell-1,2,i)=nnz(sel_S2);
            end
        case 3
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                trialNum(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1))=nnz(sel_S1);
                trialNum(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1))=nnz(sel_S2);
            end
        case 4
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                trialNum(idx:idx+ncell-1,1,i-2*(ceil(i/2)-1))=nnz(sel_S1);
                trialNum(idx:idx+ncell-1,2,i-2*(ceil(i/2)-1))=nnz(sel_S2);
            end
        case 6
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                trialNum(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1),ceil(i/3))=nnz(sel_S1);
                trialNum(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1),ceil(i/3))=nnz(sel_S2);
            end
        case 12
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                trialNum(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1),ceil(i/6),ceil(i/3))=nnz(sel_S1);
                trialNum(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1),ceil(i/6),ceil(i/3))=nnz(sel_S2);
            end
    end
    if any(reshape(trialNum(idx:idx+ncell-1,:,:,:),1,[])<10)
        avail=[avail;zeros(ncell,1)];
    else
        avail=[avail;ones(ncell,1)];
    end
    idx=idx+ncell;    
end
trialNum(~avail,:,:,:)=[];
trialNumMax=max(trialNum(:));
switch(length(trial_type))
    case 2
        out1=nan(length(NeuronSet.id),2,2,length(time),trialNumMax);
        out2=nan(length(NeuronSet.id),2,2,length(time));
    case 3
        out1=nan(length(NeuronSet.id),2,3,length(time),trialNumMax);
        out2=nan(length(NeuronSet.id),2,3,length(time));
    case 6
        out1=nan(length(NeuronSet.id),2,3,2,length(time),trialNumMax);
        out2=nan(length(NeuronSet.id),2,3,2,length(time));
    case 12
        out1=nan(length(NeuronSet.id),2,3,length(time),trialNumMax);
        out2=nan(length(NeuronSet.id),2,3,length(time));
end

idx=1;
while idx <= length(NeuronSet.id)
    
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));  
    
    temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
  
    temp.FR_Normailized=(temp.FR(:,time,:)-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
 
%     temp.FR_Normailized=(temp.FR(:,time,:)-mean(temp.FR(:,time,:),3))./std(temp.FR(:,time,:),0,3);
    switch(length(trial_type))
        case 2
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                out1(idx:idx+ncell-1,1,i,:,1:nnz(sel_S1))=permute(temp.FR_Normailized(:,:,sel_S1)...
                    ,[1,4,5,2,3]);
                out1(idx:idx+ncell-1,2,i,:,1:nnz(sel_S2))=permute(temp.FR_Normailized(:,:,sel_S2)...
                    ,[1,4,5,2,3]);
                out2(idx:idx+ncell-1,1,i,:)=permute(mean(temp.FR_Normailized(:,:,sel_S1),3)...
                    ,[1,3,4,2]);
                out2(idx:idx+ncell-1,2,i,:)=permute(mean(temp.FR_Normailized(:,:,sel_S2),3)...
                    ,[1,3,4,2]);
            end
        case 3
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                out1(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1),:,1:nnz(sel_S1))=permute(temp.FR_Normailized(:,:,sel_S1)...
                    ,[1,4,5,2,3]);
                out1(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1),:,1:nnz(sel_S2))=permute(temp.FR_Normailized(:,:,sel_S2)...
                    ,[1,4,5,2,3]);
                out2(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1),:)=permute(mean(temp.FR_Normailized(:,:,sel_S1),3)...
                    ,[1,3,4,2]);
                out2(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1),:)=permute(mean(temp.FR_Normailized(:,:,sel_S2),3)...
                    ,[1,3,4,2]);
            end
        case 6
            for i=1:size(trial_type,2)
                [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));
                out1(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1),ceil(i/3),:,1:nnz(sel_S1))=permute(temp.FR_Normailized(:,:,sel_S1)...
                    ,[1,4,5,6,2,3]);
                out1(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1),ceil(i/3),:,1:nnz(sel_S2))=permute(temp.FR_Normailized(:,:,sel_S2)...
                    ,[1,4,5,6,2,3]);
                out2(idx:idx+ncell-1,1,i-3*(ceil(i/3)-1),ceil(i/3),:)=permute(mean(temp.FR_Normailized(:,:,sel_S1),3)...
                    ,[1,3,4,5,2]);
                out2(idx:idx+ncell-1,2,i-3*(ceil(i/3)-1),ceil(i/3),:)=permute(mean(temp.FR_Normailized(:,:,sel_S2),3)...
                    ,[1,3,4,5,2]);
            end
    end


    idx=idx+ncell;
    clear temp
end
% out1(~avail,:,:,:,:,:)=[];
% out2(~avail,:,:,:,:)=[];
out1(isnan(out1))=0;
out2(isnan(out2))=0;
end
