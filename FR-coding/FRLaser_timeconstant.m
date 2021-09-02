clc
clearvars -except FR
close all

if ispc
    CodePath='D:\code';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end
load(fullfile(HomePath,'session_list.mat'))
load('D:\pixel-optogenetic\FR_AIopto_0330.mat','P')
% load('D:\pixel-optogenetic\FR_AIopto_0330.mat','FR')
%% Calculate the laser time constant
p=cell2mat(cellfun(@(x)(x(2,:)),P(:,1),'UniformOutput',false));


cid{1,1}=P(p(:,17)<0.05&p(:,18)>=0.05&p(:,19)>=0.05&p(:,20)>=0.05,2);
cid{1,2}=P(p(:,17)<0.05&p(:,18)<0.05&p(:,19)>=0.05&p(:,20)>=0.05,2);
cid{1,3}=P(p(:,17)<0.05&p(:,18)<0.05&p(:,19)<0.05&p(:,20)>=0.05,2);
cid{1,4}=P(p(:,17)<0.05&p(:,18)<0.05&p(:,19)<0.05&p(:,20)<0.05,2);

cid{2,1}=P(p(:,17)>=0.05&p(:,18)<0.05&p(:,19)>=0.05&p(:,20)>=0.05,2);
cid{2,2}=P(p(:,17)>=0.05&p(:,18)<0.05&p(:,19)<0.05&p(:,20)>=0.05,2);
cid{2,3}=P(p(:,17)>=0.05&p(:,18)<0.05&p(:,19)<0.05&p(:,20)<0.05,2);

cid{3,1}=P(p(:,17)>=0.05&p(:,18)>=0.05&p(:,19)<0.05&p(:,20)>=0.05,2);
cid{3,2}=P(p(:,17)>=0.05&p(:,18)>=0.05&p(:,19)<0.05&p(:,20)<0.05,2);

cid{4,1}=P(p(:,17)>=0.05&p(:,18)>=0.05&p(:,19)>=0.05&p(:,20)<0.05,2);

for a=1:4
    for j=1:size(cid(a,:),2)
        for i=1:length(cid{a,j})
            temp.trials=FR{cell2mat(FR(:,2))==cid{a,j}{i},3};
            temp.FR=FR{cell2mat(FR(:,2))==cid{a,j}{i},1};
            l=zeros(1,size(cid(a,:),2));
            l(a:j+a-1)=(mean(temp.FR(temp.trials(:,9)==2,16+a:15+a+j),1)-mean(temp.FR(temp.trials(:,9)==-1,16+a:15+a+j),1))./...,
                (mean(temp.FR(temp.trials(:,9)==2,16+a:15+a+j),1)+mean(temp.FR(temp.trials(:,9)==-1,16+a:15+a+j),1));
            FRchanged{a,j}(i,:)=l;
            %     FRchanged{j}(i,l>0)=1;
            %     FRchanged{j}(i,l<0)=-1;
            clear temp 
            
        end
    end
end
save('temp_FRchanged.mat','FRchanged','cid')
FRchanged=cellfun(@(x)(-1*x),FRchanged,'UniformOutput',false);
FR_increased=[];FR_decreased=[];
% 1 bin
for a=1:4    
    [m,n]=sort(mean(FRchanged{a,1}(:,a),2),'descend');
    FR_increased=[FRchanged{a,1}(n(m>0),:);FR_increased];
    FR_decreased=[FR_decreased;FRchanged{a,1}(n(m<0),:)];
    clear m n
end
%2 bin
for a=1:3    
    [m,n]=sort(mean(FRchanged{a,2}(:,a:a+1),2),'descend');
    FR_increased=[FRchanged{a,2}(n(m>0),:);FR_increased];
    FR_decreased=[FR_decreased;FRchanged{a,2}(n(m<0),:)];
    clear m n
end
%3 bin
for a=1:2    
    [m,n]=sort(mean(FRchanged{a,3}(:,a:a+2),2),'descend');
    FR_increased=[FRchanged{a,3}(n(m>0),:);FR_increased];
    FR_decreased=[FR_decreased;FRchanged{a,3}(n(m<0),:)];
    clear m n
end
%4 bin
for a=1:1    
    [m,n]=sort(mean(FRchanged{a,4}(:,a:a+3),2),'descend');
    FR_increased=[FRchanged{a,4}(n(m>0),:);FR_increased];
    FR_decreased=[FR_decreased;FRchanged{a,4}(n(m<0),:)];
    clear m n
end
fh=figure('Color','w','Position',[100,100,250,800]);
colormap([[[0:0.05:0.9;0:0.05:0.9]' ones(19,1)];[1,1,1];[ones(19,1) [0.9:-0.05:0;0.9:-0.05:0]']]);
imagesc([FR_increased;FR_decreased])
caxis([-1 1])

saveas(gcf,fullfile(HomePath,'Laser-modulated-heatmap'),'png')




%% Main function
function out=CalculatetimeConstant(FR,cid)
trials=FR{cid,3};
out=(mean(FR{cid,1}(trials(:,9)==-1,:),1)-mean(FR{cid,1}(trials(:,9)==2,:),1));
end

function  out=PlotOneTrack(rootpath,trials)
FR_All=h5read(fullfile(rootpath,'FR_All.hdf5'),'/FR_All');
for cid=1:size(FR_All,3)      
    out(cid,:)=(mean(FR_All(trials(:,9)==-1,29:38,cid),1)-mean(FR_All(trials(:,9)==2,29:38,cid),1));
        
   
%     out2{cid,2}=FR_laserOFF;
end
end

function out=CalculatePvalue(FR,trials)

end