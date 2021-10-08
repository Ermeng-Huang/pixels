clc
clear
close all

[CodePath,HomePath]=util.ComputerPath('dualtask');

% addpath(fullfile(CodePath,'+util'))

%% Decoding parameter
opt.rps=100;
opt.IsSameNeuronSet=true; %traning and test for same neuron set
opt.decoder='SVM';
opt.parallel=false;
% opt.learning='L'; % L(learning)/ W(well-trained)/ all phase(default)
opt.region='PL'; % all(default)
opt.sel='all'; %%% sustained/transient/sel/nonsel or d1/d2/d1d2
% opt.opto='none';
opt.trial='correct-distratorNo'; 
%%% laseroff-correct/laseroff-error/laseron-correct/laseron-error or correct-distratorGo
opt.trialnum=20;
opt.neuronnum=250;
prefix=sprintf('%s_%d_%s_%d',opt.region,opt.neuronnum,opt.trial,opt.rps);   
%% Decoding neuron set 
disp('NeuronSet')
tic
if isfile(fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.mat',prefix)))
    load(fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.mat',prefix)))
else    
    %%% Train   
    % opt.train.opto='none'; % laser-N,N means time bin following laser, 0-laser time bin
    % opt.train.opto0=true; % including neuron in laser-0 ??
    
%     NeuronSet_train=util.ChooseNeuronSet_dualtask('region',opt.region,'regionLevel',5,'sel',opt.sel);
    %%%
    load(fullfile(HomePath,'decoding','NeuronSet','PL_250_correct-distratorNoGo_100_decdata.mat'),'decdata')
    NeuronSet_train=decdata.NeuronSet_train;
    %%%
%     NeuronSet_train=util.ChooseNeuronSet('learning',opt.learning,'region',opt.region,'regionLevel',5,'sel',opt.sel,'opto',opt.opto);
    %%% NeuronSet(ID,Path,Region) for decoding
    [decdata.train,decdata.NeuronSet_train]=NeuronSetForDecoding(HomePath,NeuronSet_train,'trial',opt.trial,'rps',opt.rps...
        ,'trialnum',opt.trialnum,'neuronnum',opt.neuronnum);
    
    %%% Testing        
    if opt.IsSameNeuronSet
        decdata.test=decdata.train;
        %     opt.test=opt.train;
    else       
        opt.test.trial='laseroff-error';
        [decdata.test,~]=NeuronSetForDecoding(decdata.NeuronSet_train,'trial',opt.test.trial,'rps',opt.rps);
    
    end
    save(fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.mat',prefix)),'decdata','opt','-v7.3')
    writetable(struct2table(opt),fullfile(HomePath,'decoding','NeuronSet',sprintf('%s_decdata.csv',prefix)))
end
toc
%% Decoding
% cvcorr=Decoding(decdata,'decoder','SVM');
% save(fullfile(HomePath,'decoding',sprintf('%s_decoding.mat',prefix)),'opt','cvcorr','-v7.3')
% 
% figure
% hold on
% plot(cell2mat(cellfun(@mean,cvcorr.cvcorr,'UniformOutput',false)),'k')
% plot(cell2mat(cellfun(@mean,cvcorr.shufcorr,'UniformOutput',false)),'k--')
% arrayfun(@(x)plot(repmat(x,1,2),[0.4,1],'k--'),[12.5,16.5,40.5,44.5])
% saveas(gcf,fullfile(HomePath,'decoding',sprintf('%s_decoding',prefix)),'png')

%% CrossTimeDecoding
cvcorr=CrossTimeDecoding(decdata,'decoder','SVM');
save(fullfile(HomePath,'decoding',sprintf('%s_CTD.mat',prefix)),'opt','cvcorr','-v7.3')

fh=figure('Color','w','Position',[100,100,250,200]);
colormap('jet');
imagesc(cell2mat(cellfun(@mean,cvcorr(1).cvcorr,'UniformOutput',false)));
set(gca,'YDir','normal'); 
xlim([0.5,size(cvcorr.cvcorr,1)+0.5])
ylim([0.5,size(cvcorr.cvcorr,1)+0.5])
colorbar
caxis([0 1])
hold on
arrayfun(@(x)plot([0.5,size(cvcorr.cvcorr,1)+.5],repmat(x,1,2),'k--'),[12.5,16.5,24.5,28.5,48.5,52.5])
arrayfun(@(x)plot(repmat(x,1,2),[0.5,size(cvcorr.cvcorr,1)+0.5],'k--'),[12.5,16.5,24.5,28.5,48.5,52.5])
contour(1:size(cvcorr.cvcorr,1),1:size(cvcorr.cvcorr,1),cvcorr.P<0.001/56,[1 1],'-w','linewidth',1)
saveas(gcf,fullfile(HomePath,'decoding',sprintf('%s_CTD',prefix)),'png')


%% Function
function [out,NeuronSet1]=NeuronSetForDecoding(HomePath,NeuronSet,opt)
arguments  
    HomePath (1,:) char
    NeuronSet (1,1) struct
    opt.trial (1,:) char 
    opt.rps (1,:) double 
    opt.NeuronNumConstrcit (1,:) logical = true
    opt.trialnum (1,1) double
    opt.neuronnum (1,1) double    
end

out.s1=[];
out.s2=[];
idx=1;
while idx <= length(NeuronSet.id)    
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
   
    temp.trials=markLPerf(h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'));
    temp.SU_id=h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');  
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    [sel_S1,sel_S2]=util.ExtractTrial(temp,'task','dualtask','trials',opt.trial); 
    if isempty(sel_S1)||isempty(sel_S2)||nnz(sel_S1)<20||nnz(sel_S2)<20
        trial_good(idx:idx+ncell-1)=0;
        idx=idx+ncell;
        clear sel_S1 sel_S2 FR temp 
        continue
    end 
    trial_good(idx:idx+ncell-1)=1;
    idx=idx+ncell;
    
end

if length(NeuronSet.id) > opt.neuronnum && opt.NeuronNumConstrcit
    id_good=find(trial_good);
    id=id_good(randperm(length(id_good),opt.neuronnum)); 
    NeuronSet1.id=NeuronSet.id(id(1:opt.neuronnum));
    NeuronSet1.path=NeuronSet.path(id(1:opt.neuronnum));   
else
    NeuronSet1=NeuronSet;
end
idx=1;
while idx <= length(NeuronSet1.id)    
    if ~ispc
        NeuronSet1.path{idx}=replace(NeuronSet1.path{idx},'\','/');
    end
   
    temp.FR_All=h5read(fullfile(HomePath,'DataSum',NeuronSet1.path{idx},'FR_All_250ms.hdf5'),'/FR_All');    
    temp.trials=markLPerf(h5read(fullfile(HomePath,'DataSum',NeuronSet1.path{idx},'FR_All_250ms.hdf5'),'/Trials'));
    temp.SU_id=h5read(fullfile(HomePath,'DataSum',NeuronSet1.path{idx},'FR_All_250ms.hdf5'),'/SU_id');  
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet1.id(strcmp(NeuronSet1.path,NeuronSet1.path{idx})),10000)));
    [sel_S1,sel_S2]=util.ExtractTrial(temp,'task','dualtask','trials',opt.trial); 
    if isempty(sel_S1)||isempty(sel_S2)||nnz(sel_S1)<20||nnz(sel_S2)<20
        idx=idx+ncell;
        clear sel_S1 sel_S2 FR temp s1 s2
        continue
    end
    temp.FR=temp.FR_All(:,:,ismember(temp.SU_id,rem(NeuronSet1.id(strcmp(NeuronSet1.path,NeuronSet1.path{idx})),10000)));
%     FR=(temp.FR-mean(temp.FR,1))./repmat(std(temp.FR),size(temp.FR,1),1,1);
    temp.FR_base=reshape(temp.FR([sel_S1;sel_S2],1:10,:),10*(nnz(sel_S1)+nnz(sel_S2)),ncell);
    FR=(temp.FR-reshape(mean(temp.FR_base,1),1,1,ncell))./reshape(std(temp.FR_base,0,1),1,1,ncell);
    
    temp.s1=reshape(arrayfun(@(x)FR(x,:,:),cell2mat(arrayfun(@(x)sel_S1(x),randi([1,length(sel_S1)],opt.trialnum,opt.rps)...
        ,'UniformOutput',false)),'UniformOutput',false),[],1);
    temp.s2=reshape(arrayfun(@(x)FR(x,:,:),cell2mat(arrayfun(@(x)sel_S2(x),randi([1,length(sel_S2)],opt.trialnum,opt.rps)...
        ,'UniformOutput',false)),'UniformOutput',false),[],1);
    s1=permute(reshape(cat(1,temp.s1{:}),opt.trialnum,opt.rps,size(FR,2),[]),[4,1,3,2]);
    s2=permute(reshape(cat(1,temp.s2{:}),opt.trialnum,opt.rps,size(FR,2),[]),[4,1,3,2]);
    out.s1=cat(1,out.s1,s1);
    out.s2=cat(1,out.s2,s2);
    
    
%     FR=temp.FR_All(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
%     
% %     FR=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:2,:),ncell,[]),0,2);
% %     FR=temp.FR_Normailized(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
%     
%     temp.s1=reshape(arrayfun(@(x)FR(x,:,:),cell2mat(arrayfun(@(x)sel_S1(x),randi([1,length(sel_S1)],opt.trialnum,opt.rps)...
%         ,'UniformOutput',false)),'UniformOutput',false),[],1);
%     temp.s2=reshape(arrayfun(@(x)FR(x,:,:),cell2mat(arrayfun(@(x)sel_S2(x),randi([1,length(sel_S2)],opt.trialnum,opt.rps)...
%         ,'UniformOutput',false)),'UniformOutput',false),[],1);
%     s1=permute(reshape(cat(1,temp.s1{:}),opt.trialnum,opt.rps,size(FR,2),[]),[1,3,2,4]);
%     s2=permute(reshape(cat(1,temp.s2{:}),opt.trialnum,opt.rps,size(FR,2),[]),[1,3,2,4]);
%     [s1,s2]=norm(s1,s2);
%     out.s1=cat(1,out.s1,permute(s1,[4,1,2,3]));
%     out.s2=cat(1,out.s2,permute(s2,[4,1,2,3]));
    idx=idx+ncell;
    clear sel_S1 sel_S2 FR temp s1 s2
end
% if length(NeuronSet.id) > opt.neuronnum && opt.NeuronNumConstrcit
%     id=sort(randperm(length(NeuronSet.id),length(NeuronSet.id)-opt.neuronnum),'ascend');
%     NeuronSet.id(id)=[];
%     NeuronSet.path(id)=[];
%     NeuronSet.reg(id)=[];    
% end
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

function [s1,s2]=norm(s1,s2)
reshape(cat(1,s1(:,1:10,:,:),s2(:,1:10,:,:)),[2*size(s1,1)*10,size(s1,3),size(s1,4)]);
mm=reshape(mean(reshape(cat(1,s1(:,1:10,:,:),s2(:,1:10,:,:)),[2*size(s1,1)*10,size(s1,3),size(s1,4)]),1),[1,1,size(s1,3),size(s1,4)]);
ss=reshape(std(reshape(cat(1,s1(:,1:10,:,:),s2(:,1:10,:,:)),[2*size(s1,1)*10,size(s1,3),size(s1,4)]),0,1),[1,1,size(s1,3),size(s1,4)]);
s1=(s1-mm)./ss;
s2=(s2-mm)./ss;
end
