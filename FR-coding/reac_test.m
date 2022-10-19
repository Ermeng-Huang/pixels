% coding re-actived in test phase
clear
clc

homedir='F:\pixel-dualtask\';

%% calculate pev to find reactived neurons
rpts=100;
trial_type=["distractorNo-correct","distractor-correct"];
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','all');
% NeuronSet=trial_limit(homedir,NeuronSet,trial_type);
for t=1:length(trial_type)
    for rpt=1:rpts
        FR=NeuronSetForPEV(homedir,NeuronSet,trial_type(t),true,false);
        for i=1:size(FR,1)
            [~,tbl]=anova1(reshape(FR(i,:,:),1,[]),repmat([0,1],1,size(FR,3)),'off');
            pev{t}(i,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            if t==2
                [p,tbl]=anova1(reshape(FR(i,:,:),1,[]),rem(randperm(size(FR,3)*2,size(FR,3)*2),2),'off');
                pev{3}(i,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            end
        end
    end
end

for i=1:22685%size(FR,1)
    p(i,1)=statistics.permutationTest(pev{2}(i,:)',pev{1}(i,:)',500,'sidedness','larger'); 
    p(i,2)=statistics.permutationTest(pev{2}(i,:)',pev{3}(i,:)',500,'sidedness','larger');    
    p(i,3)=statistics.permutationTest(pev{1}(i,:)',pev{3}(i,:)',500,'sidedness','larger');    
end
NeuronSet=structfun(@(x)x(p(:,2)<0.05),NeuronSet,'UniformOutput',false);
save(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet','trial_type','pev','p')


%% plot pev
load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet')
load('F:\pixel-dualtask\pev_testphase.mat','pev1')
% nnz(p(:,2)<0.05&p(:,3)<0.05)
NeuronSet_all=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','all');
NeuronSet_sel=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');
NeuronSet_nonsel=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','nonsel');
pev1=cellfun(@(x)x(p(:,3)<0.05,:,:),pev1,'UniformOutput',false);

fh=figure('Color','w','Position',[50,50,450,215]);
subplot(2,3,1)
pie([nnz(ismember(NeuronSet2.id,NeuronSet_sel.id)),length(NeuronSet_sel.id)-nnz(ismember(NeuronSet2.id,NeuronSet_sel.id))]);
colormap([102,102,102;255,255,255]/255)
title('reac-sel')

subplot(2,3,2)
pie([nnz(ismember(NeuronSet.id,NeuronSet_nonsel.id)),length(NeuronSet_nonsel.id)-nnz(ismember(NeuronSet.id,NeuronSet_sel.id))]);
colormap([102,102,102;255,255,255]/255)
title('reac-nonsel')

subplot(2,3,3)
pie([length(NeuronSet.id),length(NeuronSet_all.id)-length(NeuronSet.id)]);
colormap([102,102,102;255,255,255]/255)
title('reac')

subplot(2,3,4)
hold on
plotLine(pev1{2}(ismember(NeuronSet.id,NeuronSet_sel.id),:,:),'r','-')
plotLine(pev1{1}(ismember(NeuronSet.id,NeuronSet_sel.id),:,:),'k','-')
% plotLine(pev1{3}(p(:,2)<0.05&ismember(NeuronSet.id,NeuronSet_sel.id),:,:)*100,'r',':')
% title(sprintf('memory neuron (n=%d)',nnz(p(:,2)<0.05&ismember(NeuronSet.id,NeuronSet_sel.id))))
arrayfun(@(x)fill([x,x+4,x+4,x],[-0.002,-0.002,0.02,0.02],'w','FaceColor',[0.97,0.72,0.18],'FaceAlpha',0.5,'EdgeColor','none'),[8.5])

subplot(2,3,5)
hold on
plotLine(pev1{2}(ismember(NeuronSet.id,NeuronSet_nonsel.id),:,:),'r','-')
plotLine(pev1{1}(ismember(NeuronSet.id,NeuronSet_nonsel.id),:,:),'k','-')
arrayfun(@(x)fill([x,x+4,x+4,x],[-0.002,-0.002,0.02,0.02],'w','FaceColor',[0.97,0.72,0.18],'FaceAlpha',0.5,'EdgeColor','none'),[8.5])

subplot(2,3,6)
hold on
plotLine(pev1{2},'r','-')
plotLine(pev1{1},'k','-')
arrayfun(@(x)fill([x,x+4,x+4,x],[-0.002,-0.002,0.02,0.02],'w','FaceColor',[0.97,0.72,0.18],'FaceAlpha',0.5,'EdgeColor','none'),[8.5])


exportgraphics(fh,fullfile(homedir,'pev_reactived-neuron.pdf'))

%% frac of reactivated-neurons cross region

%% with perfomance
NeuronSet_all=util.ChooseNeuronSet_dualtask('sel','all');
load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet')
trial_type=["distractor"];
NeuronSet1=util.ChooseNeuronSet_dualtask('sel','sel');
NeuronSet=structfun(@(x)x(ismember(NeuronSet.id,NeuronSet1.id)),NeuronSet,'UniformOutput',false);
%%% different trials
for i=1:length(trial_type)    
    [FR{i},avails(:,i)]=NeuronSetForPCA(homedir,NeuronSet,'trial',trial_type(i),'phase','test','pref',true);     %col1-prefer+correct col2-nonprefer_correct
end   

for i=1:length(trial_type)
    [aucc(i),auce(i)]=auc(FR{i},true);
end

fh=figure('Color','w','Position',[100,100,450,150]);
for i=1:2
    subplot(1,3,i)
    hold on
    histogram(FR{1}(:,2*i-1),-2:0.1:2,'FaceColor','r');
    histogram(FR{1}(:,2*i),-2:0.1:2,'FaceColor','k');
    xlabel('Normalized FR');
    ylabel('Probability')
end

subplot(1,3,3)
hold on;
hc=plot(aucc.xc,aucc.yc,'-r','LineWidth',1);
he=plot(auce.xe,auce.ye,'-k','LineWidth',1);
legend([hc,he],{sprintf('Correct trials AUC=%0.3f',aucc.auc),...
    sprintf('Error trials AUC=%0.3f',auce.auc)},'Location','southeast');
xlabel('False positive rate (fpr)');
ylabel('True positive rate (tpr)');

exportgraphics(fh,fullfile(homedir,'reactived-FR_Performance.pdf'))
%%
load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet');    

rpts=500;
trial_type=["distractor-correct","distractor-error"];

for t=1:length(trial_type)
    for rpt=1:rpts
        FR=NeuronSetForPEV(homedir,NeuronSet,trial_type(t),true,true);        
            for i=1:size(FR,1)
                [~,tbl]=anova1(reshape(FR(i,:,:),1,[]),repmat([0,1],1,size(FR,3)),'off');
                pev{t}(i,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            end           
    end
end
NeuronSet_Sel=util.ChooseNeuronSet_dualtask('sel','sel');
NeuronSet_nonsel=util.ChooseNeuronSet_dualtask('sel','nonsel');
m(1:2)=cellfun(@(x)nanmean(nanmean(x(ismember(NeuronSet.id,NeuronSet_Sel.id),:),2),1),pev);
m(3:4)=cellfun(@(x)nanmean(nanmean(x(ismember(NeuronSet.id,NeuronSet_nonsel.id),:),2),1),pev);
s(1:2)=cellfun(@(x)nanstd(nanmean(x(ismember(NeuronSet.id,NeuronSet_Sel.id),:),2),1)/nnz(ismember(NeuronSet.id,NeuronSet_Sel.id)),pev);
s(3:4)=cellfun(@(x)nanstd(nanmean(x(ismember(NeuronSet.id,NeuronSet_nonsel.id),:),2),1)/nnz(ismember(NeuronSet.id,NeuronSet_nonsel.id)),pev);
fh=figure('Color','w','Position',[100,100,214,214]);
hold on
bar(1:4,m*100,'w')
errorbar(1:4,m*100,s*100,'k','LineStyle','none')
set(gca,'XTick',1.5:2:4,'XTickLabel',{'memory','non-memory'},'YLim',[2,4])
ylabel('PEV')

p(1)=anova1(cell2mat(cellfun(@(x)nanmean(x(ismember(NeuronSet.id,NeuronSet_Sel.id),:),2),pev,'UniformOutput',false)')....
    ,kron([1,2],ones(1,nnz(ismember(NeuronSet.id,NeuronSet_Sel.id))))','display','off');
p(2)=anova1(cell2mat(cellfun(@(x)nanmean(x(ismember(NeuronSet.id,NeuronSet_nonsel.id),:),2),pev,'UniformOutput',false)')....
    ,kron([1,2],ones(1,nnz(ismember(NeuronSet.id,NeuronSet_nonsel.id))))','display','off');

exportgraphics(fh,fullfile(homedir,'reactive-correct-error.pdf'))

%%

m=cellfun(@(x)nanmean(x,2),pev,'UniformOutput',false);

anova1([m{1};m{2}],[zeros(length(m{1}),1);ones(length(m{2}),1)])
fh=figure('Color','w','Position',[100,100,214,214]);
hold on
bar(1:4,m*100,'w')
errorbar(1:4,m*100,s*100,'k','LineStyle','none')
set(gca,'XTick',1.5:2:4,'XTickLabel',{'memory','non-memory'},'YLim',[2,4])
ylabel('PEV')

p(1)=anova1(cell2mat(cellfun(@(x)nanmean(x(ismember(NeuronSet.id,NeuronSet_Sel.id),:),2),pev,'UniformOutput',false)')....
    ,kron([1,2],ones(1,nnz(ismember(NeuronSet.id,NeuronSet_Sel.id))))','display','off');
p(2)=anova1(cell2mat(cellfun(@(x)nanmean(x(ismember(NeuronSet.id,NeuronSet_nonsel.id),:),2),pev,'UniformOutput',false)')....
    ,kron([1,2],ones(1,nnz(ismember(NeuronSet.id,NeuronSet_nonsel.id))))','display','off');

exportgraphics(fh,fullfile(homedir,'reactive-correct-error.pdf'))



%% FR histogram
fh=figure('Color','w','Position',[100,100,428,428]);
sgtitle(sprintf('averaged cross-trial'));

for i=1:length(trial_type) 
    subplot(2,3,i)
    hold on
    histogram(FR{i}(:,1),-2:0.1:2,'FaceColor','r');
    histogram(FR{i}(:,2),-2:0.1:2,'FaceColor','k');
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')
    
    subplot(2,3,i+3)
    hold on
    histogram(FR{i}(:,3),-2:0.1:2,'FaceColor','r');
    histogram(FR{i}(:,4),-2:0.1:2,'FaceColor','k');
    
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')
end
exportgraphics(fh,fullfile(homedir,'FR_Performance','LD-InnerTask.pdf'))
%% Main function
function [out,avails]=NeuronSetForPCA(homedir,NeuronSet,opt)
arguments
    homedir (1,:) char 
    NeuronSet (1,1) struct
    opt.trial (1,:) char  = 'all'
    opt.phase (1,:) char    
    opt.pref (1,:) logical = false
end
if strcmp(opt.phase,'d1')
    time=5:6;
elseif strcmp(opt.phase,'d2')
    time=8:9;
elseif strcmp(opt.phase,'d3')
    time=12;
elseif strcmp(opt.phase,'test')
    time=13;

else
    time=5:12;
end
[temp1.FR1,temp1.FR2,temp1.FR3,temp1.FR4]=deal([]); avails=[];
idx=1;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    
    temp.FR=mean(permute(temp.FR_all(:,time,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]),2);
    
    [temp.trial1,temp.trial2]=ExtractTrial(temp.trials,sprintf('%s-correct',opt.trial));
    [temp.trial3,temp.trial4]=ExtractTrial(temp.trials,sprintf('%s-error',opt.trial));
    if isempty(temp.trial1)|| isempty(temp.trial2)||isempty(temp.trial3)|| isempty(temp.trial4)...
            || nnz(temp.trial4)<5 || nnz(temp.trial4)<5 || nnz(temp.trial4)<5 || nnz(temp.trial4)<5
        
        idx=idx+ncell;
        avails=[avails;zeros(ncell,1)];
        [temp1.FR1(end+1:end+ncell,1),temp1.FR2(end+1:end+ncell,1),temp1.FR3(end+1:end+ncell,1),temp1.FR4(end+1:end+ncell,1)]=deal(zeros(ncell,1));
        clear temp
        continue
    end
    [temp.trial5,temp.trial6]=util.ExtractTrial(temp.trials,'task','dualtask','trials','distractor-correct');
    delaymm=mean(cat(3,temp.FR(:,:,temp.trial5),temp.FR(:,:,temp.trial6)),3);
    delaystd=std(cat(3,temp.FR(:,:,temp.trial5),temp.FR(:,:,temp.trial6)),0,3);
%     delaymm=mean(cat(3,temp.FR(:,1:2,temp.trial1),temp.FR(:,1:2,temp.trial2)),3);
%     delaystd=std(cat(3,temp.FR(:,1:2,temp.trial1),temp.FR(:,1:2,temp.trial2)),0,3);
    

    temp1.FR1(end+1:end+ncell,1)=(mean(temp.FR(:,:,temp.trial1),3)-delaymm)./delaystd;
    temp1.FR2(end+1:end+ncell,1)=(mean(temp.FR(:,:,temp.trial2),3)-delaymm)./delaystd;
    temp1.FR3(end+1:end+ncell,1)=(mean(temp.FR(:,:,temp.trial3),3)-delaymm)./delaystd;
    temp1.FR4(end+1:end+ncell,1)=(mean(temp.FR(:,:,temp.trial4),3)-delaymm)./delaystd;
    
    idx=idx+ncell;
    avails=[avails;delaystd~=0];
    clear temp    
end
temp1.FR=[temp1.FR1,temp1.FR2,temp1.FR3,temp1.FR4]; % S1-c,S2-c,S1-e,S1-e
if opt.pref
    if strcmp(opt.phase,'test')
        for i=1:size(temp1.FR,1)
            if ~avails(i)
                continue
            end
            if temp1.FR(i,1)>=temp1.FR(i,2)
                out(i,:)=temp1.FR(i,:);
            else
                out(i,:)=temp1.FR(i,[2,1,4,3]);
            end
        end    
        
    else     
    per_sample=NeuronSet.pref;
    try
        for i=1:size(per_sample,1)
            if ~avails(i)
                continue
            end
            out(i,1)=temp1.FR(i,per_sample(i,:));
            out(i,2)=temp1.FR(i,2/per_sample(i,:));
            out(i,3)=temp1.FR(i,2+per_sample(i,:));
            out(i,4)=temp1.FR(i,2+2/per_sample(i,:));
        end
    catch
        d
    end
    end
else
    out=temp1.FR;

end
% id=NeuronSet.id(logical(avails),1);
end

function NeuronSet=trial_limit(homedir,NeuronSet,trial_type)
idx=1;
avails=[];
while idx <= length(NeuronSet.id)
    trial_min=[];
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    for t=1:length(trial_type)
        [sel_S1,sel_S2]=ExtractTrial(temp.trials,trial_type(t));
        trial_min=min([trial_min,length(sel_S1),length(sel_S2)]);
    end
    if trial_min<5
        avails=[avails;zeros(ncell,1)];
    else
        avails=[avails;ones(ncell,1)];
    end
    idx=idx+ncell;   
    clear temp 
end
NeuronSet=structfun(@(x)x(logical(avails)),NeuronSet,'UniformOutput',false);
end

function [out,NeuronSet]=NeuronSetForPEV(homedir,NeuronSet,trial_type,single_neuron,Isdiff_phase)
out=[];
idx=1;
trial_min=5;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    out0=[];
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));  
    
    if Isdiff_phase
%         diff_phase=[17,24;29,36;45,48];       
        temp.FR_all=mean(temp.FR_all(:,49:52,:),2); 
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.FR_Normailized=temp.FR;
    else
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.FR_Normailized=temp.FR(:,41:56,:);%[-1,11]
    end   
    
%     temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
    
    [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type);
    try
    if single_neuron
        out=cat(1,out,[temp.FR_Normailized(:,:,sel_S1(randperm(numel(sel_S1),trial_min))),temp.FR_Normailized(:,:,sel_S2(randperm(numel(sel_S2),trial_min)))]);
    else
        out=cat(1,out,[mean(temp.FR_Normailized(:,:,sel_S1),3),mean(temp.FR_Normailized(:,:,sel_S2),3)]);
    end
    catch
        disp('aiza')
    end
    idx=idx+ncell;
    clear temp
end

end

function [sel_S1,sel_S2]=ExtractTrial(in,trial_type)
if strcmp(trial_type,'Inner-correct-correct')
    sel_S1=find(in(:,5)== 4 & all(in(:,end-1:end),2) & in(:,9)~=-1 & xor(in(:,9)==in(:,11),in(:,13)==1));
    sel_S2=find(in(:,5)== 8 & all(in(:,end-1:end),2) & in(:,9)~=-1 &  xor(in(:,9)==in(:,11),in(:,13)==1));
elseif strcmp(trial_type,'Inner-correct-error')
    sel_S1=find(in(:,5)== 4 & in(:,14)==0 & in(:,9)~=-1 &  xor(in(:,9)==in(:,11),in(:,13)==1));
    sel_S2=find(in(:,5)== 8 & in(:,14)==0 & in(:,9)~=-1 &  xor(in(:,9)==in(:,11),in(:,13)==1));
elseif strcmp(trial_type,'Inner-error-correct')
    sel_S1=find(in(:,5)== 4 & all(in(:,end-1:end),2) & in(:,9)~=-1 &  ~xor(in(:,9)==in(:,11),in(:,13)==1));
    sel_S2=find(in(:,5)== 8 & all(in(:,end-1:end),2) & in(:,9)~=-1 &  ~xor(in(:,9)==in(:,11),in(:,13)==1));
elseif strcmp(trial_type,'Inner-error-error')
    sel_S1=find(in(:,5)== 4 & in(:,14)==0 & in(:,9)~=-1 &  ~xor(in(:,9)==in(:,11),in(:,13)==1));
    sel_S2=find(in(:,5)== 8 & in(:,14)==0 & in(:,9)~=-1 &  ~xor(in(:,9)==in(:,11),in(:,13)==1));
elseif strcmp(trial_type,'distractorNo-correct')
    sel_S1=find(in(:,5)==4 & in(:,9)== -1 & all(in(:,end-1:end),2));
    sel_S2=find(in(:,5)==8 & in(:,9)== -1 & all(in(:,end-1:end),2));
elseif strcmp(trial_type,'distractorNo-error')
    sel_S1=find(in(:,5)==4 & in(:,9)== -1 & in(:,14)== 0);
    sel_S2=find(in(:,5)==8 & in(:,9)== -1 & in(:,14)== 0);   
elseif strcmp(trial_type,'distractor-correct')
    sel_S1=find(in(:,5)==4 & in(:,9)~= -1 & in(:,14)== 1 & in(:,15)== 1);
    sel_S2=find(in(:,5)==8 & in(:,9)~= -1 & in(:,14)== 1 & in(:,15)== 1);
elseif strcmp(trial_type,'distractor-error')
    sel_S1=find(in(:,5)==4 & in(:,9)~= -1 & in(:,14)== 0);
    sel_S2=find(in(:,5)==8 & in(:,9)~= -1 & in(:,14)== 0);
end
end
            
function [out1,out2]=auc(in,IsError)
if iscell(in)
    [out1.xc,out1.yc,~,out1.auc]=perfcurve([zeros(size(in{:,1})),ones(size(in{:,2}))],[in{:,1},in{:,2}],0);
elseif ismatrix(in)
    onecolumn=size(in,1);
    [out1.xc,out1.yc,~,out1.auc]=perfcurve((1:2*onecolumn)>onecolumn,[in(:,1);in(:,2)],0);
end
if IsError
    if iscell(in)
        [out2.xe,out2.ye,~,out2.auc]=perfcurve([zeros(size(in{:,3})),ones(size(in{:,4}))],[in{:,3},in{:,4}],0);
    elseif ismatrix(in)
        onecolumn=size(in,1);
        [out2.xe,out2.ye,~,out2.auc]=perfcurve((1:2*onecolumn)>onecolumn,[in(:,3);in(:,4)],0);
    end
else
    out2=[];
end

end

function out=judgePerf(in) 
if nnz(in(:,end)&in(:,9)==-1)<10 || nnz(in(:,end)&in(:,9)~=-1&xor(in(:,9)==in(:,11),in(:,13)==1))<10 || nnz(in(:,end)&in(:,9)~=-1&~xor(in(:,9)==in(:,11),in(:,13)==1))<10
    out=nan(1,3);
    return
end
out(1)=nnz(all(in(:,end-1:end),2) & in(:,9)==-1)/nnz(in(:,end) & in(:,9)==-1);
out(2)=nnz(all(in(:,end-1:end),2) & in(:,9)~=-1 & xor(in(:,9)==in(:,11),in(:,13)==1))...
    /nnz(in(:,end) & in(:,9)~=-1 & xor(in(:,9)==in(:,11),in(:,13)==1));
out(3)=nnz(all(in(:,end-1:end),2) & in(:,9)~=-1 & ~xor(in(:,9)==in(:,11),in(:,13)==1))...
    /nnz(in(:,end) & in(:,9)~=-1 & ~xor(in(:,9)==in(:,11),in(:,13)==1));


end

function plotLine(in,color,LineStyle)
for i=1:size(in,1)
    for bin=1:size(in,2)
        out(i,bin)=nansum(in(i,bin,:),3)/nnz(~isnan(in(i,bin,:)));
    end
    
end
m=nanmean(out,1);
% s=nanstd(out,0,1)/sqrt(sum(~isnan(out),1));
plot(1:length(m),m,'color',color,'LineStyle',LineStyle)
ci=bootci(1000,{@(x) nanmean(x,1),out},'type','normal');
fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],color,'EdgeColor','none','FaceAlpha',0.2);
set(gca,'XTick',4.5:4:12.5,'XTickLabel',-1:1:1)
xlabel('time from test onset(s)')
ylabel('pev')
end