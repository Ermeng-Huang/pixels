clc;
clear;
% homedir='D:\pixel-optogenetic\';
homedir='F:\pixel-dualtask\';
%% Dataset
NeuronSet=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel','pref',true);
[FR,~]=NeuronSetForPCA_dualtask(homedir,NeuronSet);
NeuronSet_dist=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel','pref',true);

%% CD
cd_Sample=-1*CodingDimension(FR(:,17:48),FR(:,85:116));
% cd_distractor=CodingDimension(FR_norm(:,425:456),FR_norm(:,493:524));
cd_distractor1=CodingDimension(FR(:,425:432),FR(:,493:500)); %Early delay
cd_distractor2=CodingDimension(FR(:,433:444),FR(:,501:512)); % Middle delay
cd_distractor3=CodingDimension(FR(:,445:448),FR(:,513:516)); %late delay

% cd_distractor=CodingDimension(FR(:,433:456),FR(:,501:524));
cd_Sample=(1-2*(NeuronSet.pref-1)).*abs(cd_Sample);
% cd_distractor1=(1-2*(NeuronSet.prefer-1)).*abs(cd_distractor1);
% cd_distractor2=(1-2*(NeuronSet.prefer-1)).*abs(cd_distractor1);
% cd_distractor3=(1-2*(NeuronSet.prefer-1)).*abs(cd_distractor1);

%%
fh=figure('Color','w','Position',[100,100,900,300]);
c={'k','b','r'};
subplot(1,3,1)
view(340,40)
hold on
plot3(1:32,FR(:,17:48).'*cd_Sample,FR(:,17:48).'*cd_distractor1,'-','Color',c{1},'LineWidth',1.5);
plot3(1:32,FR(:,85:116).'*cd_Sample,FR(:,85:116).'*cd_distractor1,'--','Color',c{1},'LineWidth',1.5);
plot3(1:32,FR(:,153:184).'*cd_Sample,[FR(:,153:160).'*cd_distractor1;FR(:,161:172).'*cd_distractor2;FR(:,173:184).'*cd_distractor3],'-','Color',c{2},'LineWidth',1.5);
plot3(1:32,FR(:,221:252).'*cd_Sample,[FR(:,221:228).'*cd_distractor1;FR(:,229:240).'*cd_distractor2;FR(:,241:252).'*cd_distractor3],'--','Color',c{2},'LineWidth',1.5);
plot3(1:32,FR(:,289:320).'*cd_Sample,[FR(:,289:296).'*cd_distractor1;FR(:,297:308).'*cd_distractor2;FR(:,309:320).'*cd_distractor3],'-','Color',c{3},'LineWidth',1.5);
plot3(1:32,FR(:,357:388).'*cd_Sample,[FR(:,357:364).'*cd_distractor1;FR(:,365:376).'*cd_distractor2;FR(:,377:388).'*cd_distractor3],'--','Color',c{3},'LineWidth',1.5);

% plot3(1:32,FR(:,17:48).'*cd_Sample,[FR(:,17:48).'*cd_distractor],'-','Color',c{1},'LineWidth',1.5);
% plot3(1:32,FR(:,85:116).'*cd_Sample,FR(:,85:116).'*cd_distractor,'--','Color',c{1},'LineWidth',1.5);
% plot3(1:32,FR(:,153:184).'*cd_Sample,FR(:,153:184).'*cd_distractor,'-','Color',c{2},'LineWidth',1.5);
% plot3(1:32,FR(:,221:252).'*cd_Sample,FR(:,221:252).'*cd_distractor,'--','Color',c{2},'LineWidth',1.5);
% plot3(1:32,FR(:,289:320).'*cd_Sample,FR(:,289:320).'*cd_distractor,'-','Color',c{3},'LineWidth',1.5);
% plot3(1:32,FR(:,357:388).'*cd_Sample,FR(:,357:388).'*cd_distractor,'--','Color',c{3},'LineWidth',1.5);

% plot3(1:32,FR(:,17:48).'*cd_Sample,FR(:,17:48).'*zeros(size(cd_distractor)),'-','Color',c{1},'LineWidth',1.5);
% plot3(1:32,FR(:,85:116).'*cd_Sample,FR(:,85:116).'*zeros(size(cd_distractor)),'--','Color',c{1},'LineWidth',1.5);

xlabel('time')
ylabel('sample')
zlabel('distractor')
set(gca,'Xtick',0.5:4:33,'XtickLabel',0:1:8,'xlim',[0.5,32.5])

subplot(1,3,2)
hold on
plotLine(FR(:,17:48).*cd_Sample,'-',c{1})
plotLine(FR(:,85:116).*cd_Sample,':',c{1})
plotLine(FR(:,153:184).*cd_Sample,'-',c{2})
plotLine(FR(:,221:252).*cd_Sample,':',c{2})
plotLine(FR(:,289:320).*cd_Sample,'-',c{3})
plotLine(FR(:,357:388).*cd_Sample,':',c{3})
xlabel('time')
ylabel('sample (A.U.)')
set(gca,'Xtick',0.5:4:33,'XtickLabel',0:1:8,'xlim',[0.5,32.5])

subplot(1,3,3)
hold on
plotLine(FR(:,17:48).*cd_distractor1,'-',c{1})
plotLine(FR(:,85:116).*cd_distractor1,'--',c{1});
plotLine([FR(:,153:160).*cd_distractor1,FR(:,161:172).*cd_distractor2,FR(:,173:184).*cd_distractor3],'-',c{2});
plotLine([FR(:,221:228).*cd_distractor1,FR(:,229:240).*cd_distractor2,FR(:,241:252).*cd_distractor3],'--',c{2});
plotLine([FR(:,289:296).*cd_distractor1,FR(:,297:308).*cd_distractor2,FR(:,309:320).*cd_distractor3],'-',c{3});
plotLine([FR(:,357:364).*cd_distractor1,FR(:,365:376).*cd_distractor2,FR(:,377:388).*cd_distractor3],'--',c{3});
% plot(1:32,FR(:,17:48).'*cd_distractor,'-','Color',c{1},'LineWidth',1.5);
% plot(1:32,FR(:,85:116).'*cd_distractor,'--','Color',c{1},'LineWidth',1.5);
% plot(1:32,FR(:,153:184).'*cd_distractor,'-','Color',c{2},'LineWidth',1.5);
% plot(1:32,FR(:,221:252).'*cd_distractor,'--','Color',c{2},'LineWidth',1.5);
% plot(1:32,FR(:,289:320).'*cd_distractor,'-','Color',c{3},'LineWidth',1.5);
% plot(1:32,FR(:,357:388).'*cd_distractor,'--','Color',c{3},'LineWidth',1.5);
xlabel('time (s)')
ylabel('distractor (A.U.)')
set(gca,'Xtick',0.5:4:33,'XtickLabel',0:1:8,'xlim',[0.5,32.5])

exportgraphics(fh,'CD_trajectory.pdf','ContentType','vector')

%% distance
fh=figure('Color','w','Position',[100,100,300,300]);
hold on
plot(1:32,FR(:,17:48).'*cd_Sample-FR(:,85:116).'*cd_Sample,'-','Color',c{1},'LineWidth',1.5)
plot(1:32,FR(:,153:184).'*cd_Sample-FR(:,221:252).'*cd_Sample,'-','Color',c{2},'LineWidth',1.5)
plot(1:32,FR(:,289:320).'*cd_Sample-FR(:,357:388).'*cd_Sample,'-','Color',c{3},'LineWidth',1.5)
xlabel('time (s)')
ylabel('distance (A.U.)')
set(gca,'Xtick',0.5:4:33,'XtickLabel',0:1:8,'xlim',[0.5,32.5])
exportgraphics(fh,'dist_trajectory.pdf','ContentType','vector')

%% plane of late delay
fh=figure('Color','w','Position',[100,100,250,250]);
hold on
plane_no=Plane3D([(1:4)',FR(:,45:48).'*cd_Sample,FR(:,45:48).'*cd_distractor1;(1:4)',FR(:,113:116).'*cd_Sample,FR(:,113:116).'*cd_distractor1],1,'k');
plane_nogo=Plane3D([(1:4)',FR(:,181:184).'*cd_Sample,FR(:,181:184).'*cd_distractor3;(1:4)',FR(:,249:252).'*cd_Sample,FR(:,181:184).'*cd_distractor3],1,'b');
plane_go=Plane3D([(1:4)',FR(:,317:320).'*cd_Sample,FR(:,317:320).'*cd_distractor3;(1:4)',FR(:,385:388).'*cd_Sample,FR(:,385:388).'*cd_distractor3],1,'r');
angle1=rad2deg(acos(planeAngle(plane_no(:,1),plane_no(:,2),plane_nogo(:,1),plane_nogo(:,2))));
angle2=rad2deg(acos(planeAngle(plane_no(:,1),plane_no(:,2),plane_go(:,1),plane_go(:,2))));
angle3=rad2deg(acos(planeAngle(plane_nogo(:,1),plane_nogo(:,2),plane_go(:,1),plane_go(:,2))));
box off
grid off
exportgraphics(fh,'plane_trajectory.pdf','ContentType','vector')

%% function
function [out,NeuronSet]=NeuronSetForPCA_dualtask(homedir,NeuronSet)
out=[];
idx=1;
while idx <= length(NeuronSet.id)
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
    out0=[];
    temp.FR_all=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/FR_All');     %FR_All_100ms_norm.hdf5
    temp.trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    temp.SU_id=h5read(fullfile(homedir,'DataSum',NeuronSet.path{idx},'FR_All_250ms.hdf5'),'/SU_id');    
    ncell=nnz(ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000)));  
    
    temp.FR=permute(temp.FR_all(:,:,ismember(temp.SU_id,rem(NeuronSet.id(strcmp(NeuronSet.path,NeuronSet.path{idx})),10000))),[3,2,1]);
    temp.FR_Normailized=(temp.FR-mean(temp.FR,3))./std(temp.FR,0,3);
%     temp.FR_Normailized=(temp.FR-mean(reshape(temp.FR(:,1:10,:),ncell,[]),2))./std(reshape(temp.FR(:,1:10,:),ncell,[]),0,2);
    trial_type=["distractorNo-correct","distractorNoGo-correct","distractorGo-correct","InnerTask-all-correct"];
    for i=1:size(trial_type,2)
        [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',trial_type(i));        
        out0=cat(2,out0,cat(2,mean(temp.FR_Normailized(:,:,sel_S1),3),mean(temp.FR_Normailized(:,:,sel_S2),3)));
    end
    out=cat(1,out,out0);
    idx=idx+ncell;
    clear temp
end
out(isnan(out))=0;
end


function cdDelay=CodingDimension(in1,in2)
cdMat=in1-in2;
cdDelay=mean(cdMat,2);
cdDelay(isnan(cdDelay))=0;
cdDelay=cdDelay/norm(cdDelay);
% [qq,rr]=qr(cdDelay);
% cdfree=[in1,in2].'*qq(:,2:end);
end

function plotLine(in,l,c)
plot(1:size(in,2),nansum(in,1),'color',c,'LineStyle',l)
ci=bootci(1000,{@(x)nansum(x,1),in},'type','normal');
fill([1:size(in,2),fliplr(1:size(in,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.2);
end

function out=Plane3D(score,num_plane,color)
if num_plane==1
    eigvecsUp = pca(score);
    plotPlane(score,eigvecsUp,[],color)
    out=eigvecsUp;
else
    eigvecsUp = pca(score(1:size(score,1)/2,:));
    eigvecsDn = pca(score(size(score,1)/2+1:end,:));
    
    %compute the angle between the upper and lower plane of best fit
    out=planeAngle(eigvecsUp(:,1),eigvecsUp(:,2),eigvecsDn(:,1),eigvecsDn(:,2));
    
    % fprintf(1,'cosine of the angle between the upper and lower planes is %.2f\n',cosTheta(i))
    n=reshape(1:size(score,1),[],4);
    plotPlane(score([n(:,1),flip(n(:,2)),n(:,3),flip(n(:,4))],:),eigvecsUp,eigvecsDn)
end
end


function cosTheta = planeAngle(v1,v2,v3,v4)
%compute the angle between the planes spanned by v1&v2 and v3&v4
%MP 2019

cosTheta = blade(v1,v2,v3,v4) ./ sqrt( blade(v1,v2,v1,v2) .* blade(v3,v4,v3,v4) );
end

function plotPlane(Y,eigvecsUp,eigvecsDn,color)
%camera position
cp = [306.7153 -362.0354 70.4390;594.3579 33.1912 114.6168];

if isempty(eigvecsDn)
    if strcmp(color,'k')
        cols=repmat([linspace(0.1,0.8,size(Y,1)/2)',linspace(0.1,0.8,size(Y,1)/2)',linspace(0.1,0.8,size(Y,1)/2)'],2,1);
    elseif strcmp(color,'r')
        cols=repmat([ones(size(Y,1)/2,1),linspace(0.1,0.8,size(Y,1)/2)',linspace(0.1,0.8,size(Y,1)/2)'],2,1);
    elseif strcmp(color,'b')
        cols=repmat([linspace(0.1,0.8,size(Y,1)/2)',linspace(0.1,0.8,size(Y,1)/2)',ones(size(Y,1)/2,1)],2,1);
    end
    for i = 1:size(Y,1)
        if i<=size(Y,1)/2 %upper
            ls = 'o';
            mfc = cols(i,:);
        else %lower
            ls = 'v';
            mfc = cols(i,:);
        end
        hold on
        plot3(Y(i,1),Y(i,2),Y(i,3),ls,'MarkerSize',5,'color',cols(i,:),'MarkerFaceColor',mfc);
    end
    %connect the markers at each location
    
%     arrayfun(@(x,y)plot3([Y(x,1),Y(y,1)],[Y(x,2) Y(y,2)],[Y(x,3) Y(y,3)],'-','color',[.5 .5 .5]),1:size(Y,1)/2,circshift(1:size(Y,1)/2,-1))
%     arrayfun(@(x,y)plot3([Y(x,1),Y(y,1)],[Y(x,2) Y(y,2)],[Y(x,3) Y(y,3)],'-','color',[.5 .5 .5]),size(Y,1)/2+1:size(Y,1),circshift(size(Y,1)/2+1:size(Y,1),-1))
    
    %plot plane of best fit for each location
    scaleFactor = 25;
    
    mn = mean(Y,1);
    x = eigvecsUp(1,1:2).*scaleFactor;
    y = eigvecsUp(2,1:2).*scaleFactor;
    z = eigvecsUp(3,1:2).*scaleFactor;
    patch([x -x]+mn(1),[y -y]+mn(2),[z -z]+mn(3),.5.*[1 1 1],'EdgeColor','none','FaceAlpha',.5);
    
  
else
    cols=[[linspace(0,0.8,size(Y,1)/4)',linspace(0,0.8,size(Y,1)/4)',ones(size(Y,1)/4,1)] ...
        ;[ones(size(Y,1)/4,1),linspace(0,0.8,size(Y,1)/4)',linspace(0,0.8,size(Y,1)/4)']];
    
    cols = [cols; flip(cols)];
    
    
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
        plot3(Y(i,1),Y(i,2),Y(i,3),ls,'MarkerSize',5,'color',cols(i,:),'MarkerFaceColor',mfc);
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
end

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
mat = [dot(v1,v3) dot(v1,v4);dot(v2,v3) dot(v2,v4)];
b = det(mat);
end
