clear
homedir= 'F:/pixel-dualtask';

%% selectivity index of outer task
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');

trial_type=["distractorNo-correct","distractor-correct"];

for i=1:length(trial_type)
    sel{i}=NeuronSetForPCA_dualtask(homedir,NeuronSet,trial_type{i},'Delay');
end

sel=cellfun(@(x)abs(x),sel,'UniformOutput',false);

fh=figure('Color','w','Position',[50,50,165,145]);
hold on
ci=bootci(1000,{@(x) mean(x),sel{1}},'type','normal');
errorbar(1:8,mean(sel{1},1),mean(sel{1},1)-ci(1,:),ci(2,:)-mean(sel{1},1),'k-')
ci=bootci(1000,{@(x) mean(x),sel{2}},'type','normal');
errorbar(1:8,mean(sel{2},1),mean(sel{2},1)-ci(1,:),ci(2,:)-mean(sel{2},1),'r-')

%% selectivity index of outer task vs inner task
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','d2-mixed-sel');
trial_type=["distractor-correct","InnerTask-all-correct"];
for i=1:length(trial_type)
    sel(:,i)=NeuronSetForPCA_dualtask(homedir,NeuronSet,trial_type{i},'MD');
end

sel(sel<0)=-sel(sel<0);
[r,p]=corr(sel(:,1),sel(:,2),'type','pearson');

fh=figure('Color','w','Position',[50,50,165,145]);
hold on
scatter(sel(:,1),sel(:,2),5,'filled','MarkerFaceColor','k')
plot([0,1],[0,1],'r--','LineWidth',1)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1,'FontSize',6)
xlabel('selectivity index','FontSize',6)
ylabel('selectivity index','FontSize',6)
exportgraphics(fh,fullfile(homedir,'selectivity_index_d2_mixed.pdf'),'ContentType','vector')

%%
function out=NeuronSetForPCA_dualtask(homedir,NeuronSet,trial_type,phase)
out=[];
idx=1;

while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    out0=[];
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));  
    if strcmp(phase','MD')
        temp.FR=mean(temp.FR_all(:,29:36,:),2);
    elseif strcmp(phase,'Delay')    
        temp.FR=temp.FR_all(:,5:12,:);
    end
    temp.FR=permute(temp.FR(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
    temp.FR_Normailized=temp.FR+0.5;
    
    
%     temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
    
    [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type);

    out=cat(1,out,(mean(temp.FR_Normailized(:,:,sel_S1),3)-mean(temp.FR_Normailized(:,:,sel_S2),3))./(mean(temp.FR_Normailized(:,:,sel_S1),3)+mean(temp.FR_Normailized(:,:,sel_S2),3)));

    idx=idx+ncell;
    clear temp
end

end
