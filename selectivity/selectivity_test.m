% Only calculated selectivity of each bin
% Following code is selectivity.py

clear
clc
homedir='F:\pixel-dualtask';


fl=dir(fullfile(homedir,'DataSum','*','*','FR_All_250ms.hdf5'));
% fl=dir('F:\neupix\DataSum\*\*\FR_All.hdf5');

addpath('D:\code\jpsth\npy-matlab-master\npy-matlab')
addpath('D:\code\fieldtrip-20200320')
ubegin=1;
uend=206;
sust_trans=[];path=[];
for i=1:uend
    disp(i)
    onefile=fl(i);
    [sust_trans_temp,path_temp]=FRwithSample(onefile.folder,onefile.name,homedir);
    sust_trans=[sust_trans;sust_trans_temp];
    path=[path;path_temp];
end
save(fullfile(homedir,sprintf('sust_trans_%d_%d_test.mat',ubegin,uend)),'sust_trans','ubegin','uend','i','path')

%% main function

function [out,path]=FRwithSample(rootpath,filename,homedir)
path=cell(0);
if isfile(fullfile(homedir,'session_list.mat'))
    load(fullfile(homedir,'session_list.mat'));
else
    session=struct2cell(dir(fullfile(homedir,'DataSum','M*')))';  
    session(:,2:6)=[];
end
FR_All=h5read(fullfile(rootpath,filename),'/FR_All');
SU_id=h5read(fullfile(rootpath,filename),'/SU_id')+100000*find(ismember(session,regexp(rootpath,'(?<=DataSum\\)(\w*)','match','once')))...
    +10000*str2double(regexp(rootpath,'(?<=imec)(\d)','match','once'));
trials=util.markLPerf(h5read(fullfile(rootpath,filename),'/Trials'),0.7,0.8,120,'dualtask');
%     trials=h5read(fullfile(rootpath,filename),'/Trials'); %all trials for calculating

if isempty(trials)
    out=[];
    return
end
% dualtask
if nnz(trials(:,5)==4&trials(:,14)==1&trials(:,15)==1&trials(:,9)==-1)<3 ...
        && nnz(trials(:,5)==8&trials(:,14)==1&trials(:,15)==1&trials(:,9)==-1)<3
    out=[];
    return
end

%%% criterion for different trial type
trials1=trials(:,6)==4&all(trials(:,end-1:end)==1,2); % test1
trials2=trials(:,6)==8&all(trials(:,end-1:end)==1,2); % test2
% trials1=trials(:,15)==1&trials(:,9)==2&trials(:,13)==1;
% trials2=trials(:,15)==1&trials(:,9)==16&trials(:,13)==-1;
if nnz(trials1)<5 && nnz(trials2)<5
    out=[];
    return
end
delay=round((trials(1,2)-trials(1,1))/30000)-1;

for i=1:size(SU_id,1)
    for t=1:size(squeeze(FR_All(:,:,i)),2)/4
        FR_S1(:,t)=mean(squeeze(FR_All(trials1,4*t-3:4*t,i)),2); 
        FR_S2(:,t)=mean(squeeze(FR_All(trials2,4*t-3:4*t,i)),2);
    end
    [out(i,1)]=TestSecondSelectivity(FR_S1(:,13),FR_S2(:,13),'sample',trials(trials1,:),trials(trials2,:),'anova1');
    [out(i,2)]=TestSecondSelectivity(FR_S1(:,13),FR_S2(:,13),'sample',trials(trials1,:),trials(trials2,:),'ranksum');
    [out(i,3)]=TestSecondSelectivity(FR_S1(:,13),FR_S2(:,13),'sample',trials(trials1,:),trials(trials2,:),'anova2');
    clear FR_S1 FR_S2
    path=[path;rootpath(regexp(rootpath,'(?<=DataSum\\)(\w*)'):end)];
end
out=[out,SU_id];
  
end

function [out]=TestSecondSelectivity(FR_S1,FR_S2,time,trial1,trial2,method)

if strcmp(time,'delay')      
    for i=1:size(FR_S1,2)
        a=FR_S1(:,i);
        b=FR_S2(:,i);        
        if strcmp(method,'ranksum')
            if sum(a(:))+sum(b(:))>0               
                p=ranksum(a,b,'method','approximate','tail','both');
                if mean(a(:))>mean(b(:))
                    out(:,i)=p;
                elseif mean(a(:))<mean(b(:))
                    out(:,i)=-p;
                else
                    out(:,i)=1;
                end
            else
                out(:,i)=1;
            end
        else
            if sum(a(:))+sum(b(:))>0
                g1=[ones(size(a));2*ones(size(b))]; %sample
                g2=ones(size(a,1)+size(b,1),1); %distractor
                g2([trial1(:,9)==2;trial2(:,9)==2])=2;
                g2([trial1(:,9)==16;trial2(:,9)==16])=3;
                p=anovan([a;b],{g1,g2},'model','interaction','display','off');
                if mean(a(:))>mean(b(:))
                    out(:,i)=p(1);
                    out(:,i+8)=p(3);
                elseif mean(a(:))<mean(b(:))
                    out(:,i)=-p(1);
                    out(:,i+8)=-p(3);
                else
                    out(:,i)=-p(1);
                    out(:,i+8)=-p(3);
                    %                 out(:,i)=1;
                    %                 out(:,i+8)=1;
                end
            else
                out(:,i)=1;
                out(:,i+8)=1;
            end
        end
    end     
elseif strcmp(time,'sample') || strcmp(time,'test') || strcmp(time,'baseline')
    a=mean(FR_S1,2);
    b=mean(FR_S2,2);
    if strcmp(method,'ranksum')
        if sum(a(:))+sum(b(:))>0
            p=ranksum(a,b,'method','approximate','tail','both');
            if mean(a(:))>mean(b(:))
                out=p;
            elseif mean(a(:))<mean(b(:))
                out=-p;
            else
                out=1;
            end
        else
            out=1;
        end
    elseif strcmp(method,'anova1')        
        if sum(a(:))+sum(b(:))>0
            g1=[ones(size(a));2*ones(size(b))]; %sample            
            p=anovan([a;b],g1,'display','off');
            if mean(a(:))>mean(b(:))
                out(:,1)=p(1);
            elseif mean(a(:))<mean(b(:))
                out(:,1)=-p(1);
            else
                out(:,1)=-p(1);
                %                 out(:,i)=1;
                %                 out(:,i+8)=1;
            end
        else
            out(:,1)=1;
        end
    elseif strcmp(method,'anova2')        
        if sum(a(:))+sum(b(:))>0
            g1=[ones(size(a));2*ones(size(b))]; %sample
            g2=ones(size(a,1)+size(b,1),1); %distractor
            g2([trial1(:,9)==2;trial2(:,9)==2])=2;
            g2([trial1(:,9)==16;trial2(:,9)==16])=3;
            p=anovan([a;b],{g1,g2},'model','interaction','display','off');
            if mean(a(:))>mean(b(:))
                out(:,1)=p(1);
            elseif mean(a(:))<mean(b(:))
                out(:,1)=-p(1);
            else
                out(:,1)=-p(1);
                %                 out(:,i)=1;
                %                 out(:,i+8)=1;
            end
        else
            out(:,1)=1;
        end
        
    end
%     if sum(FR_S1(:))+sum(FR_S2(:))>0
%         p=ranksum(mean(FR_S1,2),mean(FR_S2,2));
%         if mean(FR_S1(:))>mean(FR_S2(:))
%             out=p;
%         elseif mean(FR_S1(:))<mean(FR_S2(:))
%             out=-p;
%         else
%             out=1;
%         end
%     else
%         out=1;
%     end
end
end


function [ZStatistics,ShuffledZStatistics,IsTransientDelayFR1,P,IsSignificant1,IsSignificantLess1]=TestSecondSelectivityChange_lite(Samp1TrialsFR...
    ,Samp2TrialsFR,ShuffleTimes,Permuted_BinID,DelayBinID)
addpath('D:\code\CQDecoding')

Samp1TrialsFR = Samp1TrialsFR(:,Permuted_BinID);
Samp2TrialsFR = Samp2TrialsFR(:,Permuted_BinID);
%% construct difference score matrix for each trial
% Stable and dynamic coding for working memory in primate prefrontal cortex,
% Eelke Spaak,Kei Watanabe,Shintaro Funahashi,and Mark G. Stokes
% umconstruct the FR difference score matrix for each trial and each bin during delay period: nTrialNum X DelayBinNum X DelayBinN
BinNum=size(Samp1TrialsFR,2);
Sam1DiffScore=zeros(size(Samp1TrialsFR,1),BinNum,BinNum);
for i=1:size(Samp1TrialsFR,1)%go through each trial
    tempTrialBinnedFR=Samp1TrialsFR(i,:);
    for iBin=1:BinNum
        Sam1DiffScore(i,iBin,:)=tempTrialBinnedFR-tempTrialBinnedFR(iBin);
    end
end
Sam2DiffScore=zeros(size(Samp2TrialsFR,1),BinNum,BinNum);
for i=1:size(Samp2TrialsFR,1)%go through each trial
    tempTrialBinnedFR=Samp2TrialsFR(i,:);
    for iBin=1:BinNum
        Sam2DiffScore(i,iBin,:)=tempTrialBinnedFR-tempTrialBinnedFR(iBin);
    end
end
%% construct the Z-stastics (based on ranksum test) for each nTimeBins X nTimeBins matrix
ZStatistics=zeros(BinNum,BinNum);
for i=1:BinNum
    for j=1:BinNum
        if j~=i
            [~,~,stats] =ranksum(Sam1DiffScore(:,i,j),Sam2DiffScore(:,i,j),'method','approximate');
            if isfield(stats,'zval') && ~isnan(stats.zval)
                ZStatistics(i,j)=stats.zval;
            else
                ZStatistics(i,j)=0; 
            end
        end
    end
end
%% construct the shuffled difference score and Z-statistics by permuting the time labels for ShuffleTimes times

ShuffledZStatistics=zeros(ShuffleTimes,BinNum,BinNum);
for iShuffleTimes=1:ShuffleTimes
    tempShuffledZStatistics=ComputeFRDiffScore(Samp1TrialsFR,Samp2TrialsFR,BinNum);
    ShuffledZStatistics(iShuffleTimes,:,:)= tempShuffledZStatistics;
end


DelayBinIndex=find(Permuted_BinID>=min(DelayBinID)&Permuted_BinID<=max(DelayBinID));
%% perform permutation test
[IsSignificant1,IsSignificantLess1,P]=PermutationTest(ZStatistics,ShuffledZStatistics,1);
DelayIsSignificant=IsSignificant1(DelayBinIndex,DelayBinIndex);
DelayIsSignificantLess=IsSignificantLess1(DelayBinIndex,DelayBinIndex);
if sum(sum(DelayIsSignificant+DelayIsSignificantLess))-sum(diag(DelayIsSignificant))-sum(diag(DelayIsSignificantLess))>0
    IsTransientDelayFR1=1;
else
    IsTransientDelayFR1=0;
end
end

