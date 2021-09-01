clc;
clear;
homedir='F:\pixel-dualtask\';

%% cosine of the angle between the different planes (4a)
% NeuronSet=util.ChooseNeuronSet('learning','L','sel','transient','opto','laser-modulated');% AI-opto
NeuronSet=util.ChooseNeuronSet_dualtask(homedir);

% FR=NeuronSetForPCA(homedir,NeuronSet,false);
FR=NeuronSetForPCA_dualtask(homedir,NeuronSet,false);

% [coeff,score]=pca(cell2mat(struct2cell(FR)')','NumComponents',5); % whole trial
% for i=1:14
%     planeID=cell2mat(arrayfun(@(x)x+[i*4-4:i*4],[1:68:272],'UniformOutput',false));
%     nDim = 3;
%     eigvecsUp = pca( score( planeID(1:8), 1:nDim ) );
%     eigvecsDn = pca( score( planeID(9:16), 1:nDim ) );
%     
%     %compute the angle between the upper and lower plane of best fit
%     cosTheta(i) = planeAngle(eigvecsUp(:,1),eigvecsUp(:,2),eigvecsDn(:,1),eigvecsDn(:,2));
%     fprintf(1,'cosine of the angle between the upper and lower planes is %.2f\n',cosTheta(i))
% end


[coeff,score]=pca(FR(:,[17:24,73:80,29:36,85:92])','NumComponents',3);
    
% plotPCA(score)

nDim = 3;
eigvecsUp = pca( score( 1:16, 1:nDim ) );
eigvecsDn = pca( score( 17:32, 1:nDim ) );

%compute the angle between the upper and lower plane of best fit
f = figure('units','inches','position',[1 1 5 5]);
cosTheta = planeAngle(eigvecsUp(:,1),eigvecsUp(:,2),eigvecsDn(:,1),eigvecsDn(:,2));
% fprintf(1,'cosine of the angle between the upper and lower planes is %.2f\n',cosTheta(i))
plotPlane(score([1:8,flip(9:16),17:24,flip(25:32)],:),eigvecsUp,eigvecsDn)
title(sprintf('delay%d, cosine of angle is %1.2f',i-4,cosTheta))

% for i=4:11
%     
% [coeff,score]=pca(cell2mat(cellfun(@(x)x(:,i*4-3:i*4),struct2cell(FR)','UniformOutput',false))','NumComponents',5);
% % plotPCA(score)
% 
% nDim = 3;
% eigvecsUp = pca( score( 1:8, 1:nDim ) );
% eigvecsDn = pca( score( 9:16, 1:nDim ) );
% 
% %compute the angle between the upper and lower plane of best fit
% subplot(2,4,i-3)
% cosTheta(i-3) = planeAngle(eigvecsUp(:,1),eigvecsUp(:,2),eigvecsDn(:,1),eigvecsDn(:,2));
% % fprintf(1,'cosine of the angle between the upper and lower planes is %.2f\n',cosTheta(i))
% plotPlane(score([1:4,flip(5:8),9:12,flip(13:16)],:),eigvecsUp,eigvecsDn)
% title(sprintf('delay%d, cosine of angle is %1.2f',i-4,cosTheta(i-3)))
% end
exportgraphics(f,fullfile(homedir,'plane','plane-3D.png'))
exportgraphics(f,fullfile(homedir,'plane','plane-3D.pdf'))
%% Project to laser-off subspace (4d)
f = figure('units','inches','position',[1 1 15 7.5]);
for i=4:11
[coeff,score]=pca(cell2mat(cellfun(@(x)x(:,i*4-3:i*4),struct2cell(FR)','UniformOutput',false))','NumComponents',5);
% plotPCA(score)

nDim = 3;
eigvecsUp = pca( score( 1:8, 1:nDim ) );
eigvecsDn = pca( score( 9:16, 1:nDim ) );


%compute the angle between the upper and lower plane of best fit
subplot(2,4,i-3)
hold on
project1=[score(1:4,1:nDim)*eigvecsDn(:,1:2);flip(score(5:8,1:nDim)*eigvecsDn(:,1:2))];
project2=[score(9:12,1:nDim)*eigvecsDn(:,1:2);flip(score(13:16,1:nDim)*eigvecsDn(:,1:2))];
arrayfun(@(x,y)plot([project1(x,1),project1(y,1)],[project1(x,2),project1(y,2)],'-','color',[.5 .5 .5]),1:8,circshift(1:8,-1))
arrayfun(@(x,y)plot([project2(x,1),project2(y,1)],[project2(x,2),project2(y,2)],'-','color',[.5 .5 .5]),1:8,circshift(1:8,-1))

cols=[[linspace(0,0.8,4)',linspace(0,0.8,4)',ones(4,1)];flip([ones(4,1),linspace(0,0.8,4)',linspace(0,0.8,4)'])];
arrayfun(@(x)plot(project1(x,1),project1(x,2),'o','color',cols(x,:),'MarkerFaceColor',cols(x,:)),1:8)    
arrayfun(@(x)plot(project2(x,1),project2(x,2),'v','color',cols(x,:),'MarkerFaceColor',cols(x,:)),1:8)    

title(sprintf('delay%d',i-4))
xlabel('laser off subspaces PC1')
ylabel('laser off subspaces PC2')
end
exportgraphics(f,fullfile(homedir,'plane','plane-2D.png'))
exportgraphics(f,fullfile(homedir,'plane','plane-2D.pdf'))

%% Function 


function [out,NeuronSet]=NeuronSetForPCA(homedir,NeuronSet,IsRpt)
temp1.FR1=cell(1,100);temp1.FR2=cell(1,100);temp1.FR3=cell(1,100);temp1.FR4=cell(1,100);
idx=1;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All.hdf5'),'/Trials');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    if IsRpt
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.T{1}=find(temp.trials(:,5)==4&temp.trials(:,9)==2);
        temp.T{2}=find(temp.trials(:,5)==8&temp.trials(:,9)==2);
        temp.T{3}=find(temp.trials(:,5)==4&temp.trials(:,9)==-1);
        temp.T{4}=find(temp.trials(:,5)==8&temp.trials(:,9)==-1);
        temp.T0=mat2cell(cell2mat(cellfun(@(y)cell2mat(arrayfun(@(x)y(x),randi([1,length(y)],10,100),'UniformOutput',false))...
            ,temp.T,'UniformOutput',false)'),40,ones(1,100));
        temp.FR_Normailized=cellfun(@(y)(y-mean(reshape(y(:,1:10,:),ncell,[]),2))./std(reshape(y(:,1:10,:),ncell,[]),0,2),cellfun(@(x)temp.FR(:,:,x),temp.T0,'UniformOutput',false),'UniformOutput',false);     
        temp1.FR1=cellfun(@(x,y)cat(1,x,y),temp1.FR1,cellfun(@(x)mean(x(:,:,1:10),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);
        temp1.FR2=cellfun(@(x,y)cat(1,x,y),temp1.FR2,cellfun(@(x)mean(x(:,:,11:20),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);
        temp1.FR3=cellfun(@(x,y)cat(1,x,y),temp1.FR3,cellfun(@(x)mean(x(:,:,21:30),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);        
        temp1.FR4=cellfun(@(x,y)cat(1,x,y),temp1.FR4,cellfun(@(x)mean(x(:,:,31:40),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);

    else
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
        out.FR1=cat(1,out.FR1,mean(temp.FR_Normailized(:,:,temp.trials(:,5)==4&temp.trials(:,9)==2),3)); %S1 and laseron
        out.FR2=cat(1,out.FR2,mean(temp.FR_Normailized(:,:,temp.trials(:,5)==8&temp.trials(:,9)==2),3)); %S2 and laseron
        out.FR3=cat(1,out.FR3,mean(temp.FR_Normailized(:,:,temp.trials(:,5)==4&temp.trials(:,9)==-1),3)); %S1 and laseroff
        out.FR4=cat(1,out.FR4,mean(temp.FR_Normailized(:,:,temp.trials(:,5)==8&temp.trials(:,9)==-1),3)); %S2 and laseroff
    end
    idx=idx+ncell;
    clear temp
end
out=cellfun(@(a,b,c,d)[a,b,c,d],temp1.FR1,temp1.FR2,temp1.FR3,temp1.FR4,'UniformOutput',false);
end

function [temp1,NeuronSet]=NeuronSetForPCA_dualtask(homedir,NeuronSet,IsRpt)
temp1.FR1=[];temp1.FR2=[];temp1.FR3=[];temp1.FR4=[];
% temp1.FR1=cell(1,100);temp1.FR2=cell(1,100);temp1.FR3=cell(1,100);temp1.FR4=cell(1,100);
idx=1;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'));
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));
    if IsRpt
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.T{1}=find(temp.trials(:,5)==4&temp.trials(:,9)==2);
        temp.T{2}=find(temp.trials(:,5)==8&temp.trials(:,9)==2);
        temp.T{3}=find(temp.trials(:,5)==4&temp.trials(:,9)==-1);
        temp.T{4}=find(temp.trials(:,5)==8&temp.trials(:,9)==-1);
        temp.T0=mat2cell(cell2mat(cellfun(@(y)cell2mat(arrayfun(@(x)y(x),randi([1,length(y)],10,100),'UniformOutput',false))...
            ,temp.T,'UniformOutput',false)'),40,ones(1,100));
        temp.FR_Normailized=cellfun(@(y)(y-mean(reshape(y(:,1:10,:),ncell,[]),2))./std(reshape(y(:,1:10,:),ncell,[]),0,2),cellfun(@(x)temp.FR(:,:,x),temp.T0,'UniformOutput',false),'UniformOutput',false);     
        temp1.FR1=cellfun(@(x,y)cat(1,x,y),temp1.FR1,cellfun(@(x)mean(x(:,:,1:10),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);
        temp1.FR2=cellfun(@(x,y)cat(1,x,y),temp1.FR2,cellfun(@(x)mean(x(:,:,11:20),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);
        temp1.FR3=cellfun(@(x,y)cat(1,x,y),temp1.FR3,cellfun(@(x)mean(x(:,:,21:30),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);        
        temp1.FR4=cellfun(@(x,y)cat(1,x,y),temp1.FR4,cellfun(@(x)mean(x(:,:,31:40),3),temp.FR_Normailized,'UniformOutput',false),'UniformOutput',false);

    else
        temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
        temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
        temp1.FR1=cat(1,temp1.FR1,mean(temp.FR_Normailized(:,:,temp.trials(:,5)==4&temp.trials(:,9)~=-1),3)); %S1 and laseron
        temp1.FR2=cat(1,temp1.FR2,mean(temp.FR_Normailized(:,:,temp.trials(:,5)==8&temp.trials(:,9)~=-1),3)); %S2 and laseron
    end
    idx=idx+ncell;
    clear temp
end
temp1=[temp1.FR1,temp1.FR2];

% temp1=cellfun(@(a,b,c,d)[a,b,c,d],temp1.FR1,temp1.FR2,temp1.FR3,temp1.FR4,'UniformOutput',false);
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
out(out(:,end)==0,:)=[];
end

function plotPCA(score)
fh=figure('Color','w','Position',[50,50,600,600]);
view([45,45])
i=reshape(1:size(score,1),[],4)';
hold on
plot3(score(i(1,:),1),score(i(1,:),2),score(i(1,:),3),'Color','r','LineStyle','-','Marker','.','LineWidth',2,'MarkerSize',5)
plot3(score(i(2,:),1),score(i(2,:),2),score(i(2,:),3),'Color','b','LineStyle','-','Marker','.','LineWidth',2,'MarkerSize',5')
plot3(score(i(3,:),1),score(i(3,:),2),score(i(3,:),3),'Color','r','LineStyle','--','Marker','.','LineWidth',2,'MarkerSize',5)
plot3(score(i(4,:),1),score(i(4,:),2),score(i(4,:),3),'Color','b','LineStyle','--','Marker','.','LineWidth',2,'MarkerSize',5')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3') 
end

function cosTheta = planeAngle(v1,v2,v3,v4)
%compute the angle between the planes spanned by v1&v2 and v3&v4
%
%MP 2019

cosTheta = blade(v1,v2,v3,v4) ./ sqrt( blade(v1,v2,v1,v2) .* blade(v3,v4,v3,v4) );

end

function plotPlane(Y,eigvecsUp,eigvecsDn)
%calculate colors for plotting data
% xe = linspace(0, 2*pi, 9);
% xc = xe(1:8) + diff(xe(1:2))/2;
% [unitX, unitY] = pol2cart(xc,30);
% cols = lab2rgb([70.*ones(numel(unitX),1),unitX',unitY']);
% cols = [cols; cols];
cols=[[linspace(0,0.8,size(Y,1)/4)',linspace(0,0.8,size(Y,1)/4)',ones(size(Y,1)/4,1)] ...
    ;[ones(size(Y,1)/4,1),linspace(0,0.8,size(Y,1)/4)',linspace(0,0.8,size(Y,1)/4)']];
cols = [cols; cols];
%camera position
cp = [
    306.7153 -362.0354   70.4390
    594.3579   33.1912  114.6168
    ];

% f = figure('units','inches','position',[1 1 10 5]);

%plot the marker for each condition 
for i = 1:size(Y,1)
    if i<=size(Y,1)/2 %upper
        ls = 'o';
        mfc = cols(i,:);
    else %lower
        ls = 'v';
        mfc = cols(i,:);
    end
    hold on
    plot3(Y(i,1),Y(i,2),Y(i,3),ls,'MarkerSize',10,'color',cols(i,:),'MarkerFaceColor',mfc);
end

%connect the markers at each location

arrayfun(@(x,y)plot3([Y(x,1),Y(y,1)],[Y(x,2) Y(y,2)],[Y(x,3) Y(y,3)],'-','color',[.5 .5 .5]),1:size(Y,1)/2,circshift(1:size(Y,1)/2,-1))
arrayfun(@(x,y)plot3([Y(x,1),Y(y,1)],[Y(x,2) Y(y,2)],[Y(x,3) Y(y,3)],'-','color',[.5 .5 .5]),size(Y,1)/2+1:size(Y,1),circshift(size(Y,1)/2+1:size(Y,1),-1))

%plot plane of best fit for each location
scaleFactor = 25;

mn = mean(Y(1:size(Y,1)/2,1:3),1);
x = eigvecsUp(1,1:2).*scaleFactor;
y = eigvecsUp(2,1:2).*scaleFactor;
z = eigvecsUp(3,1:2).*scaleFactor;
patch([x -x]+mn(1),[y -y]+mn(2),[z -z]+mn(3),.5.*[1 1 1],'EdgeColor','none','FaceAlpha',.5);

mn = mean(Y(size(Y,1)/2+1:size(Y,1),1:3),1);
x = eigvecsDn(1,1:2).*scaleFactor;
y = eigvecsDn(2,1:2).*scaleFactor;
z = eigvecsDn(3,1:2).*scaleFactor;
patch([x -x]+mn(1),[y -y]+mn(2),[z -z]+mn(3),.5.*[1 1 1],'EdgeColor','none','FaceAlpha',.5);

%etc
grid on
axis square
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
% set(gca,'xtick',[-20 0 20],'ytick',[-20 0 20],'ztick',[-20 0 20]);

campos(cp(1,:));
end
function b = blade(v1,v2,v3,v4)

mat = [dot(v1,v3) dot(v1,v4);
    dot(v2,v3) dot(v2,v4)];

b = det(mat);
end
function out=extractFR(homedir,filename)
onefile=dir(fullfile(homedir,'DataSum','*','*',filename));
out.FR1=[];out.FR2=[];out.FR3=[];out.FR4=[];
for f=onefile'
    FR_all=h5read(fullfile(f.folder,f.name),'/FR_All');
    trials=h5read(fullfile(f.folder,f.name),'/Trials');
    out.FR1=cat(1,out.FR1,permute(mean(FR_all(trials(:,5)==4&trials(:,9)==2,:,:),1),[3,2,1])); %S1 and laseron
    out.FR2=cat(1,out.FR2,permute(mean(FR_all(trials(:,5)==8&trials(:,9)==2,:,:),1),[3,2,1])); %S2 and laseron
    out.FR3=cat(1,out.FR3,permute(mean(FR_all(trials(:,5)==4&trials(:,9)==-1,:,:),1),[3,2,1])); %S1 and laseroff
    out.FR4=cat(1,out.FR4,permute(mean(FR_all(trials(:,5)==8&trials(:,9)==-1,:,:),1),[3,2,1])); %S2 and laseroff
end
end

