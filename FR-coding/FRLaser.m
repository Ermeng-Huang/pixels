clc
clear
close all

if ispc
    CodePath='D:\code';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end
load(fullfile(HomePath,'session_list.mat'))
path=util.Path_default;
%% Calculated the effect of laser for all neurons
if ~isfile(fullfile(HomePath,'FR_AIopto_0524.mat'))
    FR=[];
    P=[];
    for i=1:length(session)
        [p_temp,FR_temp]=CalculateEffectOfLaser(session{i},HomePath,'rank');
        p_temp(:,2)=cellfun(@(x)(x+i*100000),p_temp(:,2),'UniformOutput',false);
        FR_temp(:,2)=cellfun(@(x)(x+i*100000),FR_temp(:,2),'UniformOutput',false);
        P=[P;p_temp];
        FR=[FR;FR_temp];  
        
        
    end
%     h5read(fullfile(HomePath,'FR_AIopto_0331.mat'))
    save('D:\pixel-optogenetic\FR_AIopto_0524.mat','P','FR','-v7.3')
else
    load('D:\pixel-optogenetic\FR_AIopto_0524.mat','P')
end
return
%% Transfer to hdf5
p=cell2mat(cellfun(@(x)(x(2,:)),P(:,1),'UniformOutput',false));
p(p<0.05/4)=1;
p(p~=1)=0;

p=cell2mat(P(:,1));
p(p<0.05/4)=1;
p(p~=1)=0;
nnz(any(p(:,17:20),2))
p_anova=h5read(fullfile(HomePath,'LaserModulation.hdf5'),'/LaserModulation');
h5create(fullfile(HomePath,'LaserModulation.hdf5'),'/rank',size(p),'Datatype','double')
h5write(fullfile(HomePath,'LaserModulation.hdf5'),'/rank',p)
h5create(fullfile(HomePath,'LaserModulation.hdf5'),'/anova',size(p_anova),'Datatype','double')
h5write(fullfile(HomePath,'LaserModulation.hdf5'),'/anova',p_anova)
%% Function
function [out1,out2]=CalculateEffectOfLaser(session,HomePath,type)
out1=[];
out2=[];
folder=dir(fullfile(HomePath,'DataSum',session,'*','FR_All_500ms.hdf5'));
trials=h5read(fullfile(folder(1,1).folder,'FR_All_500ms.hdf5'),'/Trials');
for f=1:length(folder)
    [p,FR]=PlotOneTrack(folder(f,1).folder,trials,type);
    p(:,2)=cellfun(@(x)(x+(f-1)*10000),p(:,2),'UniformOutput',false);
    FR(:,2)=cellfun(@(x)(x+(f-1)*10000),FR(:,2),'UniformOutput',false);
    out1=[out1;p];
    out2=[out2;FR];
end
end


function  [out1,out2]=PlotOneTrack(rootpath,trials,type)
FR_All=h5read(fullfile(rootpath,'FR_All_500ms.hdf5'),'/FR_All');
for cid=1:size(FR_All,3)   
    FR=FR_All(:,:,cid);
    if strcmp(type,'anova')
        out1{cid,1}=CalculatePvalue(FR,trials,type);
%         out1{cid,1}=CalculatePvalue(squeeze(FR_All(:,:,cid)),trials);
       
    elseif strcmp(type,'glm') 
        B_laser=FitGlm(squeeze(FR_All(:,:,cid)),trials);
        
    elseif strcmp(type,'rank') 
        out1{cid,1}=CalculatePvalue(FR,trials,type);
    end
    
%     out2{cid,1}=squeeze(FR_All(:,:,cid));
    out2{cid,1}=FR;
    out1{cid,2}=cid;
    out2{cid,2}=cid;
    out2{cid,3}=trials;
end
end

function out=CalculatePvalue(FR,trials,type)
if strcmp(type,'anova')
    FR_S1_ON=FR(trials(:,9)==2&trials(:,5)==4,:);
    FR_S1_OFF=FR(trials(:,9)==-1&trials(:,5)==4,:);
    FR_S2_ON=FR(trials(:,9)==2&trials(:,5)==8,:);
    FR_S2_OFF=FR(trials(:,9)==-1&trials(:,5)==8,:);
    for i=1:size(FR_S2_OFF,2)
        
        x=[FR_S1_ON(:,i);FR_S1_OFF(:,i);FR_S2_ON(:,i);FR_S2_OFF(:,i)];
        g1=[ones(size(FR_S1_ON,1)+size(FR_S1_OFF,1),1);2*ones(size(FR_S2_ON,1)+size(FR_S2_OFF,1),1)];
        g2=[ones(size(FR_S1_ON,1),1);2*ones(size(FR_S1_OFF,1),1);ones(size(FR_S2_ON,1),1);2*ones(size(FR_S2_OFF,1),1)]; %laser
        out(:,i)=anovan(x,{g1,g2},'display','off');
    end
elseif strcmp(type,'rank')
    FR_ON=FR(trials(:,9)==2,:);
    FR_OFF=FR(trials(:,9)==-1,:);
     for i=1:size(FR_OFF,2) 
        out(:,i)=ranksum(FR_ON(:,i),FR_OFF(:,i));
    end
    
end
end

function out=FitGlm(FR,trials)
FR_S1_ON=FR(trials(:,9)==2&trials(:,5)==4,:);
FR_S1_OFF=FR(trials(:,9)==-1&trials(:,5)==4,:);
FR_S2_ON=FR(trials(:,9)==2&trials(:,5)==8,:);
FR_S2_OFF=FR(trials(:,9)==-1&trials(:,5)==8,:);
X1=[mean(FR_S1_ON(:,1:5),2);mean(FR_S1_OFF(:,1:5),2);mean(FR_S2_ON(:,1:5),2);mean(FR_S2_OFF(:,1:5),2)]; % baseline
X2=[ones(size(FR_S1_ON,1)+size(FR_S1_OFF,1),1);zeros(size(FR_S2_ON,1)+size(FR_S2_OFF,1),1)];
X3=[ones(size(FR_S1_ON,1),1);zeros(size(FR_S1_OFF,1),1);ones(size(FR_S2_ON,1),1);zeros(size(FR_S2_OFF,1),1)]; %laser
y=[FR_S1_ON(:,17);FR_S1_OFF(:,17);FR_S2_ON(:,17);FR_S2_OFF(:,17)];
b=fitglm([X1,X2,X3],y);
if table2array(b.Coefficients(4,4))<0.05
    out=table2array(b.Coefficients(4,1));
else
    out=0;
end

end  
