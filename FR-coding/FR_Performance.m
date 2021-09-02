% Compare normalized FR and AUC in correct and error trials in different group (laser, selective etc)
close all
clear
homedir='D:\pixel-optogenetic';
typesel=["sustained";"transient_6"];
type_trial=["laseron","laseroff"];
f=figure('Color','w','Position',[100,100,750,352]);
sgtitle(sprintf('averaged cross-trial'));
for i=1:2 
    NeuronSet=util.ChooseNeuronSet('learning','L','sel','transient_6','opto','none');
%     NeuronSet=util.ChooseNeuronSet('learning','L','sel',typesel(i),'opto','none');
    [FR,~]=NeuronSetForPCA(homedir,NeuronSet,'trial',type_trial(i));
    
    %%
    s(i)=subplot(2,3,3*i-2);
    hold on
    histogram(FR(:,1),-2:0.1:2,'FaceColor','r');    
    histogram(FR(:,2),-2:0.1:2,'FaceColor','k');
    xlim([-1.5,1.5]);
    title('Correct trials')
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')
    xline(nanmean(FR(:,1)),'--r','LineWidth',1);
    xline(nanmean(FR(:,2)),'--k','LineWidth',1);
    box off
    
    subplot(2,3,3*i-1)
    hold on
    histogram(FR(:,3),-2:0.1:2,'FaceColor','r');    
    histogram(FR(:,4),-2:0.1:2,'FaceColor','k');
    xline(nanmean(FR(:,3)),'--r','LineWidth',1);
    xline(nanmean(FR(:,4)),'--k','LineWidth',1);
    xlim([-1.5,1.5]);
    title('Error trials')
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')
    box off
    
    %% auc
    onecolumn=size(FR,1);
    [xc,yc,~,aucc]=perfcurve((1:2*onecolumn)>onecolumn,[FR(:,1);FR(:,2)],0);
    [xe,ye,~,auce]=perfcurve((1:2*onecolumn)>onecolumn,[FR(:,3);FR(:,4)],0);
    
    subplot(2,3,3*i);
    hold on;
    hc=plot(xc,yc,'-r','LineWidth',1);
    he=plot(xe,ye,'-k','LineWidth',1);
    legend([hc,he],{sprintf('Correct trials AUC=%0.3f',aucc),...
        sprintf('Error trials AUC=%0.3f',auce)},'Location','southeast');
    xlabel('False positive rate (fpr)');
    ylabel('True positive rate (tpr)');    
end

annotation(f,'textarrow',[0.05 0.05],[0.1+s(1).Position(2)+s(1).Position(4)/2 0.1+s(1).Position(2)+s(1).Position(4)/2]...
    ,'String','Laser-on trial','LineStyle','none','Color','w','TextColor','k','FontSize',14,'TextRotation',90);
annotation(f,'textarrow',[0.05 0.05],[0.1+s(2).Position(2)+s(2).Position(4)/2 0.1+s(2).Position(2)+s(2).Position(4)/2]...
    ,'String','Laser-off trial','LineStyle','none','Color','w','TextColor','k','FontSize',14,'TextRotation',90);
exportgraphics(f,fullfile(homedir,'FR_Performance','FR_Performance_laser.pdf'))
exportgraphics(f,fullfile(homedir,'FR_Performance','FR_Performance_laser.png'))
%% function
function [out,NeuronSet]=NeuronSetForPCA(homedir,NeuronSet,opt)
arguments
    homedir (1,:) char 
    NeuronSet (1,1) struct
    opt.trial (1,:) char  = 'all'
   
end
temp1.FR1=[];temp1.FR2=[];temp1.FR3=[];temp1.FR4=[];
idx=1;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/Trials');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_1000ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    
    temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
%     temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:2,:),ncell,[]),0,2);
    if strcmp(opt.trial,'laseron')
        temp.FR1=temp.FR(:,10,temp.trials(:,5)==4&temp.trials(:,11)==1&temp.trials(:,9)==2); %S1 and laseroff
        temp.FR2=temp.FR(:,10,temp.trials(:,5)==8&temp.trials(:,11)==1&temp.trials(:,9)==2); %S2 and laseroff
        temp.FR3=temp.FR(:,10,temp.trials(:,5)==4&temp.trials(:,11)==0&temp.trials(:,9)==2);
        temp.FR4=temp.FR(:,10,temp.trials(:,5)==8&temp.trials(:,11)==0&temp.trials(:,9)==2);
    elseif strcmp(opt.trial,'laseroff')
        temp.FR1=temp.FR(:,10,temp.trials(:,5)==4&temp.trials(:,11)==1&temp.trials(:,9)==-1); %S1 and laseroff
        temp.FR2=temp.FR(:,10,temp.trials(:,5)==8&temp.trials(:,11)==1&temp.trials(:,9)==-1); %S2 and laseroff
        temp.FR3=temp.FR(:,10,temp.trials(:,5)==4&temp.trials(:,11)==0&temp.trials(:,9)==-1);
        temp.FR4=temp.FR(:,10,temp.trials(:,5)==8&temp.trials(:,11)==0&temp.trials(:,9)==-1);
    end
    delaymm=mean(cat(3,temp.FR1,temp.FR2),3);
    delaystd=std(cat(3,temp.FR1,temp.FR2),0,3);
    if delaystd==0, continue;  end 
    temp1.FR1=cat(1,temp1.FR1,mean((temp.FR1-delaymm)./delaystd,3));
    temp1.FR2=cat(1,temp1.FR2,mean((temp.FR2-delaymm)./delaystd,3));
    temp1.FR3=cat(1,temp1.FR3,mean((temp.FR3-delaymm)./delaystd,3));
    temp1.FR4=cat(1,temp1.FR4,mean((temp.FR4-delaymm)./delaystd,3));
    
    idx=idx+ncell;
    clear temp    
end
temp1.FR=[temp1.FR1,temp1.FR2,temp1.FR3,temp1.FR4];
sus_trans=h5read(fullfile(homedir,'Selectivity_AIopto_0419.hdf5'),'/sus_trans_noPermutaion')';
clusterid=h5read(fullfile(homedir,'Selectivity_AIopto_0419.hdf5'),'/cluster_id');
per_sample=any(sus_trans(ismember(clusterid,NeuronSet.id),7:12)==1,2)+ 2*any(sus_trans(ismember(clusterid,NeuronSet.id),7:12)==2,2);
% per_sample=-((temp1.FR1-temp1.FR2)>0)+2;
for i=1:size(per_sample,1)    
    out(i,1)=temp1.FR(i,per_sample(i));
    out(i,2)=temp1.FR(i,2/per_sample(i));
    out(i,3)=temp1.FR(i,2+per_sample(i));
    out(i,4)=temp1.FR(i,2+2/per_sample(i));
end
end
