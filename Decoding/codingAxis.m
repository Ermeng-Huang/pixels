clc
clearvars -Except decdata
close all

[CodePath,HomePath]=util.ComputerPath();
Path=util.Path_default;
sus_trans=h5read(Path.selectivity,'/sus_trans_noPermutaion')';
cluster_id=h5read(Path.selectivity,'/cluster_id');
path=regexp(h5read(Path.selectivity,'/path'),'(\w|\\|-)*','match','once');

sample=find(sus_trans(:,6)~=0);
delay=find(sus_trans(:,9)~=0); %delay 3s

% NeuronSet_train.id=cluster_id(sample);
% NeuronSet_train.path=path(sample);
% [decdata.train1,~]=NeuronSetForDecoding(NeuronSet_train,'trial','laseroff','rps',100 ...
%     ,'trialnum',60,'neuronnum',100);
% NeuronSet_train.id=cluster_id(delay);
% NeuronSet_train.path=path(delay);
% [decdata.train2,~]=NeuronSetForDecoding(NeuronSet_train,'trial','laseroff','rps',100 ...
%     ,'trialnum',60,'neuronnum',100);
% decdata.train.s1=cat(1,decdata.train1.s1,decdata.train2.s1);
% decdata.train.s2=cat(1,decdata.train1.s2,decdata.train2.s2);
% decdata.test=decdata.train;

%% Decoding
cvcorr=Decoding(decdata,'decoder','SVM');
axis1=mean(cvcorr.beta(:,1:2:2000),2);
axis2=mean(cvcorr.beta(:,2:2:2000),2);
s=[axis1(1:100,:),axis2(1:100,:)];
d=[axis1(101:200,:),axis2(101:200,:)];
[s_B,s_BINT,s_R]=regress(axis1(1:100,:),axis2(1:100,:));
[d_B,d_BINT,d_R]=regress(axis1(101:200,:),axis2(101:200,:));
s_mdl = fitlm(axis1(1:100,:),axis2(1:100,:));
d_mdl = fitlm(axis1(101:200,:),axis2(101:200,:));

figure
hold on
scatter(axis1(1:100,:),axis2(1:100,:),'r')
scatter(axis1(101:200,:),axis2(101:200,:),'b')
x=-0.1:0.01:0.1;
plot(x,x*s_mdl.Coefficients{2,1}+x*s_mdl.Coefficients{2,1},'r-')
plot(x,x*d_mdl.Coefficients{2,1}+x*d_mdl.Coefficients{2,1},'b-')

save(fullfile(HomePath,'decoding',sprintf('%s_decoding.mat',prefix)),'opt','cvcorr','-v7.3')
acos(dot(cvcorr.beta{1},cvcorr.beta{2})/(norm(cvcorr.beta{1})*norm(cvcorr.beta{2})))
figure
hold on
plot(cell2mat(cellfun(@mean,cvcorr.cvcorr,'UniformOutput',false)),'k')
plot(cell2mat(cellfun(@mean,cvcorr.shufcorr,'UniformOutput',false)),'k--')
arrayfun(@(x)plot(repmat(x,1,2),[0.5,68.5],'k--'),[12.5,16.5,40.5,44.5])
saveas(gcf,fullfile(HomePath,'decoding',sprintf('%s_decoding',prefix)),'png')


plot(cvcorr.beta{1})
%% Function
function out=Decoding(decdata,opt)
%%%
% Input:
%     decdata is structure, including two structures (train and test, each including two 4-D matrix s1 and s2)
%     decdata.train.s1 = (cellID * trlN * bins * rpts) matrix
%%%
arguments
    decdata (1,1) struct
    opt.decoder (1,:) char {mustBeMember(opt.decoder,{'SVM','LDA','NB'})} ='SVM';
end

bins=size(decdata.train.s1,3);
trlN=size(decdata.train.s1,2);
rpts=size(decdata.train.s1,4);

cv=cvpartition(trlN,'KFold',10);
y=[zeros(trlN,1);ones(trlN,1)];

out=struct();

out.beta=[];
out.bias=[];
out.angle=[];
cv_trainingID=cell2mat(arrayfun(@(x)training(cv,x),[1:10],'UniformOutput',false));
cv_testID=cell2mat(arrayfun(@(x)test(cv,x),[1:10],'UniformOutput',false));

for rpt=1:rpts
    fprintf('rpts %d of %d\n',rpt,rpts);
    for kf=1:cv.NumTestSets
        for bin=[13,25]
            s1kf=mean(decdata.train.s1(:,cv_trainingID(:,kf),bin:bin+3,rpt),3);
            s2kf=mean(decdata.train.s2(:,cv_trainingID(:,kf),bin:bin+3,rpt),3);
            varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
            s1kf=s1kf(varsel,:);
            s2kf=s2kf(varsel,:);
            
            Xkf=cat(2,s1kf,s2kf)';
            ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
            
            s1Tkf=mean(decdata.test.s1(varsel,cv_testID(:,kf),bin:bin+3,rpt),3);
            s2Tkf=mean(decdata.test.s2(varsel,cv_testID(:,kf),bin:bin+3,rpt),3);
            XTkf1=cat(2,s1Tkf,s2Tkf)';
            yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);
            yshufTkf=yTkf(randperm(numel(yTkf)));            
           
            if strcmp(opt.decoder,'SVM')
                SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
                modelPredict=SVMM.predict(XTkf1);
            elseif strcmp(opt.decoder,'LDA')
                LDAM=fitcdiscr(Xkf,ykf);
                modelPredict=LDAM.predict(XTkf1);
            elseif strcmp(opt.decoder,'NB')
                NBM=fitcnb(Xkf,ykf);
                modelPredict=NBM.predict(XTkf1);
            end
            out.beta(:,end+1)=SVMM.Beta;
            out.bias(:,end+1)=SVMM.Bias; 
        end        
    end
end
end


function [out,NeuronSet]=NeuronSetForDecoding(NeuronSet,opt)
arguments
    NeuronSet (1,1) struct
    opt.trial (1,:) char
    opt.rps (1,:) double
    opt.NeuronNumConstrcit (1,:) logical = true
    opt.trialnum (1,1) double
    opt.neuronnum (1,1) double
end

Path=util.Path_default;
out.s1=[];
out.s2=[];
idx=1;
if length(NeuronSet.id) > opt.neuronnum && opt.NeuronNumConstrcit
    id=sort(randperm(length(NeuronSet.id),length(NeuronSet.id)-opt.neuronnum),'ascend');
    NeuronSet.id(id)=[];
    NeuronSet.path(id)=[];
end

while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    temp.FR_All=h5read(fullfile(Path.home,'DataSum',NeuronSet.path{idx},'FR_All.hdf5'),'/FR_All');
    temp.trials=h5read(fullfile(Path.home,'DataSum',NeuronSet.path{idx},'FR_All.hdf5'),'/Trials');
    temp.SU_id=h5read(fullfile(Path.home,'DataSum',NeuronSet.path{idx},'FR_All.hdf5'),'/SU_id');
    temp.FR_Normailized=(temp.FR_All-mean(temp.FR_All,1))./repmat(std(temp.FR_All),size(temp.FR_All,1),1,1);
    if strcmp(opt.trial,'laseroff-correct')
        sel_S1=find(temp.trials(:,5)==4 & temp.trials(:,9)== -1 & temp.trials(:,11)== 1);
        sel_S2=find(temp.trials(:,5)==8 & temp.trials(:,9)== -1 & temp.trials(:,11)== 1);
    elseif strcmp(opt.trial,'laseroff-error')
        sel_S1=find(temp.trials(:,5)==4 & temp.trials(:,9)== -1 & temp.trials(:,11)== 0);
        sel_S2=find(temp.trials(:,5)==8 & temp.trials(:,9)== -1 & temp.trials(:,11)== 0);
    elseif strcmp(opt.trial,'laseron-correct')
        sel_S1=find(temp.trials(:,5)==4 & temp.trials(:,9)== 2 & temp.trials(:,11)== 1);
        sel_S2=find(temp.trials(:,5)==8 & temp.trials(:,9)== 2 & temp.trials(:,11)== 1);
    elseif strcmp(opt.trial,'laseron-error')
        sel_S1=find(temp.trials(:,5)==4 & temp.trials(:,9)== 2 & temp.trials(:,11)== 0);
        sel_S2=find(temp.trials(:,5)==8 & temp.trials(:,9)== 2 & temp.trials(:,11)== 0);
    else
        sel_S1=find(temp.trials(:,5)==4);
        sel_S2=find(temp.trials(:,5)==8);
    end
    
    FR=temp.FR_Normailized(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    
    temp.s1=reshape(arrayfun(@(x)FR(x,:,:),cell2mat(arrayfun(@(x)sel_S1(x),randi([1,length(sel_S1)],opt.trialnum,opt.rps)...
        ,'UniformOutput',false)),'UniformOutput',false),[],1);
    temp.s2=reshape(arrayfun(@(x)FR(x,:,:),cell2mat(arrayfun(@(x)sel_S2(x),randi([1,length(sel_S2)],opt.trialnum,opt.rps)...
        ,'UniformOutput',false)),'UniformOutput',false),[],1);
    s1=permute(reshape(cat(1,temp.s1{:}),opt.trialnum,opt.rps,68,[]),[4,1,3,2]);
    s2=permute(reshape(cat(1,temp.s2{:}),opt.trialnum,opt.rps,68,[]),[4,1,3,2]);
    
    out.s1=cat(1,out.s1,s1);
    out.s2=cat(1,out.s2,s2);
    idx=idx+nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    clear sel_S1 sel_S2 FR temp s1 s2
end
end



