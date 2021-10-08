% 1 traning for n test (same neurons for different condition/result)
clc
clear
close all

[CodePath,HomePath]=util.ComputerPath();

% addpath(fullfile(CodePath,'+util'))

%% Decoding parameter
prefix='Tran100_CorrectTrialsDecodeCorrectandErrorTrial_SampleDecoding_100';
opt.rps=100;
opt.IsSameNeuronSet=true; %traning and test for same neuron set
opt.decoder='SVM';
opt.parallel=true;
opt.learning='L'; % L(learning)/ W(well-trained)/ all phase(default)
opt.region='all'; % all(default)
opt.sel='transient'; %sustained/transient/sel/nonsel
opt.opto='none';
opt.trial='laseroff-correct';

opt.trialnum=20;
opt.neuronnum=100;
    
%% Decoding neuron set 
disp('NeuronSet')
tic
if isfile(fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.mat',prefix)))
    load(fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.mat',prefix)))
else    
    %%% Train    
    
    NeuronSet_train=util.ChooseNeuronSet('learning',opt.learning,'region',opt.region,'regionLevel',5,'sel',opt.sel,'opto',opt.opto);
    %NeuronSet(ID,Path,Region) for decoding
    [decdata.train,decdata.NeuronSet_train]=NeuronSetForDecoding(NeuronSet_train,'trial',opt.trial,'rps',opt.rps...
        ,'trialnum',opt.trialnum,'neuronnum',opt.neuronnum);
    
    %%% Testing        
          
    opt.test.trial='laseroff-correct';
    [decdata.test1,~]=NeuronSetForDecoding(decdata.NeuronSet_train,'trial',opt.test.trial,'rps',opt.rps...
        ,'trialnum',opt.trialnum,'neuronnum',opt.neuronnum);
    opt.test.trial='laseroff-error';
    [decdata.test2,~]=NeuronSetForDecoding(decdata.NeuronSet_train,'trial',opt.test.trial,'rps',opt.rps...
        ,'trialnum',opt.trialnum,'neuronnum',opt.neuronnum);
    
    save(fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.mat',prefix)),'decdata','opt','-v7.3')
    writetable(struct2table(opt),fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.csv',prefix)))
end
toc
%% Decoding
cvcorr=Decoding(decdata,'decoder','SVM');
save(fullfile(HomePath,'decoding',sprintf('%s_decoding.mat',prefix)),'opt','cvcorr','-v7.3')

figure
hold on
plot(cell2mat(cellfun(@mean,cvcorr.cvcorr,'UniformOutput',false)),'k')
plot(cell2mat(cellfun(@mean,cvcorr.cvcorr2,'UniformOutput',false)),'b')
plot(cell2mat(cellfun(@mean,cvcorr.shufcorr,'UniformOutput',false)),'k--')

arrayfun(@(x)plot(repmat(x,1,2),[0.5,1],'k--'),[12.5,16.5,40.5,44.5])

saveas(gcf,fullfile(HomePath,'decoding',sprintf('%s_decoding',prefix)),'png')

%% CrossTimeDecoding
cvcorr=CrossTimeDecoding(decdata,'decoder','SVM');
save(fullfile(HomePath,'decoding',sprintf('%s_CTD.mat',prefix)),'opt','cvcorr','-v7.3')

fh=figure('Color','w','Position',[100,100,250,200]);
colormap('jet');
imagesc(cell2mat(cellfun(@mean,cvcorr(1).cvcorr,'UniformOutput',false)));
set(gca,'YDir','normal'); 
xlim([0.5,68.5])
ylim([0.5,68.5])
colorbar
caxis([0 1])
hold on
arrayfun(@(x)plot([0.5,68.5],repmat(x,1,2),'k--'),[12.5,16.5,40.5,44.5])
arrayfun(@(x)plot(repmat(x,1,2),[0.5,68.5],'k--'),[12.5,16.5,40.5,44.5])
contour(1:68,1:68,cvcorr.p,[1 1],'-r','linewidth',1)
saveas(gcf,fullfile(HomePath,'decoding',sprintf('%s_CTD',prefix)),'png')


%% Function
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
    NeuronSet.reg(id)=[];    
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
%     out.s1=cat(1,out.s1,(s1-mean(s1,2))./repmat(permute(std(permute(s1,[2,1,3,4])),[2,1,3,4]),1,size(s1,2),1,1));
%     out.s2=cat(1,out.s2,(s1-mean(s2,2))./repmat(permute(std(permute(s2,[2,1,3,4])),[2,1,3,4]),1,size(s2,2),1,1));
    idx=idx+nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    clear sel_S1 sel_S2 FR temp s1 s2
end
  
end




