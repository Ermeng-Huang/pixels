clc;
clear;
% homedir='D:\pixel-optogenetic\';
homedir='F:\pixel-dualtask\';
%% Dataset
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');
[FR_S1,FR_S2]=NeuronSetForPCA_dualtask(homedir,NeuronSet); 
%% CD
cdMat=FR_S1-FR_S2;
cdDelay=mean(cdMat,2);
cdDelay=cdDelay/norm(cdDelay);
[qq,rr]=qr(cdDelay); % non-orthonormal
cdfree=[FR_S1,FR_S2].'*qq(:,2:end); %exclude vaiernce explained by CD
[~,score,~]=pca(cdfree,'NumComponents',20);
score=mat2cell(score,ones(4,1)*size(score,1)/4,20);
% [CD_s1off,CD_s2off,CD_s1on,CD_s2on]=deal(score{1},score{2},score{3},score{4});

%% PCA
% [~,score,~]=pca(FR_All','NumComponents',20);
% score=mat2cell(score,ones(4,1)*size(score,1)/4,20);
% [PCA_s1off,PCA_s2off,PCA_s1on,PCA_s2on]=deal(score{1},score{2},score{3},score{4});

%%
fh=figure('Color','w','Position',[100,100,450,300]);
c={'k','r','b'};
% subplot(1,2,1)
view(160,-10)
hold on
for i=1:2
    plot3(score{i}(:,1),score{i}(:,2),FR_S1(:,32*(i-1)+(1:32)).'*cdDelay,'-','Color',c{i},'LineWidth',1.5);
    plot3(score{i+2}(:,1),score{i+2}(:,2),FR_S2(:,32*(i-1)+(1:32)).'*cdDelay,'--','Color',c{i},'LineWidth',1.5);
end

xlabel('PC1')
ylabel('PC2')
zlabel('CD')

% subplot(1,2,2)
% view(168,-3)
% hold on
% h1=plot3(PCA_s1off(:,1),PCA_s1off(:,2),PCA_s1off(:,3),'-c','LineWidth',1.5);
% h2=plot3(PCA_s2off(:,1),PCA_s2off(:,2),PCA_s2off(:,3),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
% h3=plot3(PCA_s1on(:,1),PCA_s1on(:,2),PCA_s1on(:,3),'-k','LineWidth',1.5);
% h4=plot3(PCA_s1on(:,1),PCA_s2on(:,2),PCA_s2on(:,3),'-r','LineWidth',1.5);
% 
% %laserOnset
% % plot3(PCA_s1off(17:18,1),PCA_s1off(17:18,2),PCA_s1off(17:18,3),'-c*','LineWidth',1.5);
% % plot3(PCA_s2off(17:18,1),PCA_s2off(17:18,2),PCA_s2off(17:18,3),'-*','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
% plot3(PCA_s1on(17:18,1),PCA_s1on(17:18,2),PCA_s1on(17:18,3),'-b*','LineWidth',1.5);
% plot3(PCA_s1on(17:18,1),PCA_s2on(17:18,2),PCA_s2on(17:18,3),'-b*','LineWidth',1.5);
% 
% legend([h1,h2,h3,h4],{'S1-LaserOff trial','S2-LaserOff trial','S1-LaserOn trial','S2-LaserOn trial'}...
%     ,'FontSize',16,'Location','eastoutside','Orientation','vertical','AutoUpdate','off')
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
exportgraphics(fh,fullfile(homedir,'CD_LaserNonMod_delay1.pdf'),'ContentType','vector')
%%
fh=figure('Color','w','Position',[100,100,450,300]);
c={'k','r','b'};
% subplot(1,2,1)
hold on
for i=1:3
    plot(1:32,FR_S1(:,32*(i-1)+(1:32)).'*cdDelay,'-','Color',c{i},'LineWidth',1.5);
    plot(1:32,FR_S2(:,32*(i-1)+(1:32)).'*cdDelay,'--','Color',c{i},'LineWidth',1.5);
end
%%
v=VideoWriter(fullfile(homedir,'cd_pc_devp_movie.mp4'),'MPEG-4');
open(v);
fh=figure('Color','w','Position',[100,100,900,600]);
hold on;
for tt=1:68
    cla   
    h1=plot3(CD_s1off(1:tt),CD_s1off(1:tt),FR_S1_laseroff(:,1:tt).'*cdDelay,'-c','LineWidth',1.5);
    h2=plot3(CD_s2off(1:tt),CD_s2off(1:tt),FR_S2_laseroff(:,1:tt).'*cdDelay,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
    h3=plot3(CD_s1on(1:tt),CD_s1on(1:tt),FR_S1_laseron(:,1:tt).'*cdDelay,'-k','LineWidth',1.5);
    h4=plot3(CD_s2on(1:tt),CD_s2on(1:tt),FR_S2_laseron(:,1:tt).'*cdDelay,'-r','LineWidth',1.5);
    legend([h1,h2,h3,h4],{'S1-LaserOff trial','S2-LaserOff trial','S1-LaserOn trial','S2-LaserOn trial'}...
        ,'FontSize',16,'Location','eastoutside','Orientation','vertical','AutoUpdate','off')    
    th1=text(CD_s1off(1:tt),CD_s1off(1:tt),FR_S1_laseroff(:,1:tt).'*cdDelay,'S1-LaserOff','Color','c','FontSize',16);
    th2=text(CD_s2off(1:tt),CD_s2off(1:tt),FR_S2_laseroff(:,1:tt).'*cdDelay,'S2-LaserOff','Color',[0.9290 0.6940 0.1250],'FontSize',16);
    th3=text(CD_s1on(1:tt),CD_s1on(1:tt),FR_S1_laseron(:,1:tt).'*cdDelay,'S1-LaserOn','Color','k','FontSize',16);
    th4=text(CD_s2on(1:tt),CD_s2on(1:tt),FR_S2_laseron(:,1:tt).'*cdDelay,'S2-LaserOn','Color','r','FontSize',16);
    pause(0.5)
    if tt<13
        tag='ITI';
    elseif tt<17
        tag='Sample';
    elseif tt<41
        tag='Delay';
    else
        tag='Test';
    end
    title(sprintf('%s (%d msec)',tag,tt*250-4000),'FontSize',20)
    xlabel('PC1')
    ylabel('PC2')
    zlabel('CD')
    for frame=1:15
        writeVideo(v,getframe(fh))
    end    
end

close(v);

%% function
function [out1,out2]=NeuronSetForPCA_dualtask(homedir,NeuronSet)
out1=[];out2=[];
idx=1;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    [out1_temp,out2_temp]=deal([]);
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));  
    
    temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
    temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
    trial_type=["distractorNo-correct","distractor-correct"];
    for i=1:size(trial_type,2)
        [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));        
        out1_temp=cat(2,out1_temp,mean(temp.FR_Normailized(:,17:48,sel_S1),3));
        out2_temp=cat(2,out2_temp,mean(temp.FR_Normailized(:,17:48,sel_S2),3));
    end
    out1=cat(1,out1,out1_temp);
    out2=cat(1,out2,out2_temp);
    idx=idx+ncell;
    clear temp
end

end

function [trial,cluster_ids,FR]=pre_process(HomePath,Sessionfolder,sessionIdx)
folder=dir(fullfile(HomePath,'DataSum',Sessionfolder,'*','FR_All.hdf5'));
trial=util.markLPerf(h5read(fullfile(folder(1,1).folder,folder(1).name),'/Trials'),0.75,0.5,60);
cluster_ids=[]; FR=[];
for f=1:size(folder,1)    
    cluster_ids=cat(1,cluster_ids,h5read(fullfile(folder(f,1).folder,folder(1).name),'/SU_id')+sessionIdx*100000 ...
        +10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))); 
    FR=cat(3,FR,h5read(fullfile(folder(f,1).folder,folder(1).name),'/FR_All'));    
end

    basemm=mean(FR(:,1:12,:),[1,2]);
    stdmm=std(FR(:,1:12,:),0,[1,2]);
    
    [sel_S1,sel_S2]=util.ExtractTrial(trial,'task','AIopto','trials','laseroff'); 
    basemm=mean(FR([sel_S1;sel_S2],:,:),1);   
    stdmm=1;
    
    FR=(FR-basemm)./stdmm;
end

function out=mem_type(id,cluster_id,sus_trans)
all_su=unique(id(:));
for i=1:numel(all_su)
    if ~any(sus_trans(cluster_id==all_su(i),1:2)) %non-mem
        mem(i)=-1;
    else
        if sus_trans(cluster_id==all_su(i),1)  %sustained
            prefer=unique(sus_trans(cluster_id==all_su(i),5:12));
            if prefer(prefer~=0)==1  %s1
                mem(i)=1;
            else
                mem(i)=3;
            end
        elseif sus_trans(cluster_id==all_su(i),2) % transient
            prefer=unique(sus_trans(cluster_id==all_su(i),5:12));
            if prefer(prefer~=0)==1 %s1
                mem(i)=2;
            else
                mem(i)=4;
            end     
        end
    end
end
out=arrayfun(@(x)mem(all_su==x),id);
end