% calculated plane for AIopto

clc;
clear;
rpts=100;
%% PEV of dualtask distractorNo sample (naive-DPA)
homedir='F:\pixel-DPA8s\';
NeuronSet=util.ChooseNeuronSet_DPA('sel','sel');
FR=NeuronSetForPCA_DPA(homedir,NeuronSet,"correct",true); 
pev=[];
for rpt=1:10
    for i=1:size(FR,1)
        for bin=1:size(FR,2)/2
            [~,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),repmat([0,1],1,size(FR,3)),'off');
            pev(i,bin)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
        end
    end
end

%% PEV of dualtask distractorNo sample
homedir='F:\pixel-dualtask\';
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');
for rpt=1:rpts    
    FR=NeuronSetForPCA_dualtask(homedir,NeuronSet,"distractorNo-error",true,true);
    for bin=1:size(FR,2)/2
        if ndims(FR)==2
            [p,tbl]=anova1(reshape(FR(:,[bin,bin+size(FR,2)/2]),1,[]),kron([0,1],ones(1,size(FR,1))),'off');
            pev_distractorNo(1,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
        else
            for i=1:size(FR,1)
                [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),repmat([0,1],1,size(FR,3)),'off');
                pev_distractorNo(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
                
                [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),rem(randperm(size(FR,3)*2,size(FR,3)*2),2),'off');
                pev_shuffle1(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            end
        end
    end
end
% 
%% PEV of dualtask distractor sample
for rpt=1:rpts
    FR=NeuronSetForPCA_dualtask(homedir,NeuronSet,"distractor-error",true,true);
    for bin=1:size(FR,2)/2
        if ndims(FR)==2
            [p,tbl]=anova1(reshape(FR(:,[bin,bin+size(FR,2)/2]),1,[]),kron([0,1],ones(1,size(FR,1))),'off');
            pev_distractor(1,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
        else
            for i=1:size(FR,1)
                [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),repmat([0,1],1,size(FR,3)),'off');
                pev_distractor(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
                
                [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),rem(randperm(size(FR,3)*2,size(FR,3)*2),2),'off');
                pev_shuffle2(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            end
        end
    end
end


%% PEV of dualtask DR sample
for rpt=1:rpts
    FR=NeuronSetForPCA_dualtask(homedir,NeuronSet,"InnerTask-all-correct",true,true);
    for bin=1:size(FR,2)/2
        if ndims(FR)==2
            [p,tbl]=anova1(reshape(FR(:,[bin,bin+size(FR,2)/2]),1,[]),kron([0,1],ones(1,size(FR,1))),'off');
            pev_DR(1,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            [p,tbl]=anova1(reshape(FR(:,[bin,bin+size(FR,2)/2]),1,[]),kron([0,1],ones(1,size(FR,1))),'off');
            pev_shuffle(1,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
        else
            for i=1:size(FR,1)
                [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),repmat([0,1],1,size(FR,3)),'off');
                pev_DR(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
                
                [p,tbl]=anova1(reshape(FR(i,[bin,bin+size(FR,2)/2],:),1,[]),rem(randperm(size(FR,3)*2,size(FR,3)*2),2),'off');
                pev_shuffle(i,bin,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            end
        end
    end
end
% for bin=1:size(FR,2)/2
%     [p,tbl]=anova1(reshape(FR(:,[bin,bin+size(FR,2)/2]),1,[]),reshape([zeros(size(FR,1),1),ones(size(FR,1),1)],1,[]),'off');
%     pev(bin)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
% end
% plot(1:length(pev),pev)

save(fullfile(homedir,'pev_diff_delay.mat'),'pev_distractor','pev_distractorNo','pev_DR','pev_shuffle','pev_shuffle1','pev_shuffle2','NeuronSet')
save(fullfile(homedir,'pev.mat'),'pev_distractor','pev_distractorNo','pev_DR','pev_shuffle','pev_shuffle1','pev_shuffle2','NeuronSet')
load(fullfile(homedir,'pev.mat'),'pev_shuffle2')
%%

fh=figure('Color','w','Position',[50,50,450,215]);
hold on


plot(1:size(pev_distractorNo,2),mean(nanmean(pev_distractorNo,1),3),'k-')
plot(1:size(pev_distractor,2),mean(nanmean(pev_distractor,1),3),'r-')
plot(1:size(pev_shuffle,2),mean(nanmean(pev_shuffle,1),3),'k:')
plot(1:size(pev_DR,2),mean(nanmean(pev_DR,1),3),'b-')
arrayfun(@(x)fill([x,x+4,x+4,x],[0,0,0.1,0.1],'r','FaceAlpha',0.5,'EdgeColor','none'),[12.5,24.5,48.5])
%% Main Function
function [out,NeuronSet]=NeuronSetForPCA_dualtask(homedir,NeuronSet,trial_type,single_neuron,Isdiff_phase)
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

function [out,NeuronSet]=NeuronSetForPCA_DPA(homedir,NeuronSet,trial_type,single_neuron)
out=[];
idx=1;
trial_min=5;
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