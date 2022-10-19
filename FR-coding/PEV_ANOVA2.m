% calculated plane for AIopto

clc;
clear;
rpts=100;
homedir='F:\pixel-dualtask\';

%% PEV of dualtask distractor sample
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');
trial_type=["distractor-correct","InnerTask-correct"];
for i=1:size(NeuronSet.id,1)    
    [FR,trial]=NeuronSetForPCA_dualtask(homedir,structfun(@(x)x(i),NeuronSet,'UniformOutput',false),"distractor-correct",true,false);
    g1=double(trial(:,5)==4);
    g2=double(trial(:,9)==2);
    
    for bin=1:size(FR,2)       
        [p,tbl]=anovan(squeeze(FR(1,bin,:)),{g1,g2},'model','interaction','display','off');
        pev_sample(i,bin)=(tbl{2,2}-tbl{2,3}*tbl{5,5})/(tbl{6,2}+tbl{5,5});
        pev_distractor(i,bin)=(tbl{3,2}-tbl{3,3}*tbl{5,5})/(tbl{6,2}+tbl{5,5});
        pev_inter(i,bin)=(tbl{4,2}-tbl{4,3}*tbl{5,5})/(tbl{6,2}+tbl{5,5});
        %             [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),rem(randperm(size(FR,3)*2,size(FR,3)*2),2),'off');
        %             pev_shuffle(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});        
    end
end
return


%%

fh=figure('Color','w','Position',[50,50,450,215]);
hold on

plotLine(pev_sample,'k','-')
plotLine(pev_distractor,'r','-')
plotLine(pev_inter,'b','-')
arrayfun(@(x)fill([x,x+4,x+4,x],[0,0,0.1,0.1],'r','FaceAlpha',0.5,'EdgeColor','none'),[4.5,16.5,40.5])
%% Main Function
function [out,trials]=NeuronSetForPCA_dualtask(homedir,NeuronSet,trial_type,single_neuron,Isdiff_phase)
out=[];
idx=1;
trial_min=10;
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
        temp.FR_all=[mean(temp.FR_all(:,17:24,:),2),mean(temp.FR_all(:,29:36,:),2),mean(temp.FR_all(:,45:48,:),2)]; 
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.FR_Normailized=temp.FR;
    else
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.FR_Normailized=temp.FR(:,9:56,:);%[-1,11]
    end   
    
%     temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
    trials=temp.trials(temp.trials(:,9)~=-1&all(temp.trials(:,end-1:end),2),:);    
    out=cat(1,out,temp.FR_Normailized(:,:,temp.trials(:,9)~=-1&all(temp.trials(:,end-1:end),2)));
    
    
    idx=idx+ncell;
    clear temp
end

end

function [out,NeuronSet]=NeuronSetForPCA_DPA(homedir,NeuronSet,trial_type,single_neuron)
out=[];
idx=1;
trial_min=10;
while idx <= length(NeuronSet.id)
    if ~ispch
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    out0=[];
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));  
    
    temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
    temp.FR_Normailized=temp.FR;
%     temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
    
    
    [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','DPA','trials',trial_type);
    if single_neuron
        out=cat(1,out,[temp.FR_Normailized(:,:,sel_S1(randperm(numel(sel_S1),trial_min))),temp.FR_Normailized(:,:,sel_S2(randperm(numel(sel_S2),trial_min)))]);
    else
        out=cat(1,out,[mean(temp.FR_Normailized(:,:,sel_S1),3),mean(temp.FR_Normailized(:,:,sel_S2),3)]);
    end
    idx=idx+ncell;
    clear temp
end

end

function pev=cal_pev(x,group)
[p,tbl]=anova1(x,group,'off');
pev=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
end

function plotLine(Data,c,l)
    m=nanmean(Data,1);
    plot(1:size(m,2),m,c,'LineStyle',l);
%     s=nanstd(Data)/sqrt(size(Data,1));
%     fill([1:size(m,2),fliplr(1:size(m,2))],[m+s,fliplr(m-s)],c,'EdgeColor','none','FaceAlpha',0.2);
    ci=bootci(1000,{@(x) nanmean(x,1),Data},'type','normal');
    fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.2);
end