%%
clc
clear
close all
addpath('D:/code/npy-matlab/npy-matlab')
addpath('D:/code/fieldtrip-20200320')
ft_defaults
homedir= 'F:/pixel-dualtask';

reg=regexp(h5read(fullfile(homedir,'Selectivity_1108.hdf5'),'/reg'),'(\w|\\|-)*','match','once');
cluster_id=h5read(fullfile(homedir,'Selectivity_1108.hdf5'),'/cluster_id');
sus_trans=h5read(fullfile(homedir,'Selectivity_1108.hdf5'),'/sust_trans_noPermutaion');
path=regexp(h5read(fullfile(homedir,'Selectivity_0907.hdf5'),'/path'),'(\w|\\|-)*','match','once');

trial_type=["distractorNo","distractorGo","distractorNoGo"];

% for typeID=1:3
%     Pvalue=h5read(fullfile(homedir,'Selectivity_0925.hdf5'),sprintf('/Pvalue_%s',trial_type(typeID)));
%     d1(:,typeID)=any(abs(Pvalue(:,3:4))<0.05/2,2);
%     d2(:,typeID)=any(abs(Pvalue(:,6:7))<0.05/2,2);
%     d3(:,typeID)=any(abs(Pvalue(:,10))<0.05/2,2);
%     s1(:,typeID)=any(abs(Pvalue(:,2))<0.05/2,2);
%     s2(:,typeID)=any(abs(Pvalue(:,5))<0.05/2,2);
% end
OBM1_index=isOBM1Map(unique(reg(:,7)));

return
% id=cluster_id(all(d1,2)&all(d2,2)&ismember(reg(:,7),OBM1_index(:,1))&s2);
%% Quick and dirty
id=cluster_id(any(sus_trans(:,1:2)~=0,2)&ismember(reg(:,7),OBM1_index(:,1)));

set(0,'DefaultFigureVisible', 'on') 
% figure('visible','off')
c=["k","r","b"];
for i=1:length(id)
    [avail,FT_SPIKE,FR,trial]=pre_process(homedir,path{cluster_id==id(i)},id(i));
    mem=util.mem_type(homedir,id(i));
    if ~avail
        continue
    end
    FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id(i)},'FR_All_1000ms.hdf5'),'/FR_All');
    SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id(i)},'FR_All_1000ms.hdf5'),'/SU_id');
    for typeID=1:3
%         [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials','InnerTask-correct');
        [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
        p(typeID,:)=arrayfun(@(x)ranksum(FR_All(S1,x,SU_id==rem(id(i),10000)),FR_All(S2,x,SU_id==rem(id(i),10000)))...
            ,1:size(FR_All,2))<0.05/8;         
    end
    if nnz(p(:,5:10))==0
       continue 
    end
    for typeID=1:3
%         [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials','InnerTask-correct');
        [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
        hold on
        plotLine(FR(S1,:),c(typeID),'-')
        plotLine(FR(S2,:),c(typeID),'--')  
        plot(find(p(typeID,:))*4,1.2*max(mean((FR([S1;S2],:)),1))*ones(nnz(p(typeID,:)),1),c(typeID),'Marker','*','LineStyle','none')
        
    end
    title(sprintf('%d %s, mem type=%d',id(i),reg{cluster_id==id(i),7},mem))
    saveas(gcf,fullfile(homedir,'showcase','dirty',sprintf('%d_%s.png',id(i),reg{cluster_id==id(i),7})))
    close
end

%%
close
id=2100143;
c=["k","r","b"];

[avail,FT_SPIKE,FR,trial]=pre_process(homedir,path{cluster_id==id},rem(id,10000));
for i=1:size(FT_SPIKE.trialtime,1)
    spike{i,:}=FT_SPIKE.time{1}(1,FT_SPIKE.trial{1,1}==i);
end
mem=util.mem_type(homedir,id);

FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/FR_All');
SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/SU_id');
for typeID=1:3   
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
    p(typeID,:)=arrayfun(@(x)ranksum(FR_All(S1,x,SU_id==rem(id,10000)),FR_All(S2,x,SU_id==rem(id,10000)))...
        ,1:size(FR_All,2))<0.05;
end

fh=figure('Color','w','Position',[100,100,500,250]);
sgtitle(sprintf('%d %s, mem type=%d',id,reg{cluster_id==id,7},mem))

for typeID=1:3     
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
    Raster1=spike(S1,1);
    Raster2=spike(S2,1);
    
    subplot(2,3,typeID)
    hold on
    plotLine(FR(S1,9:48),c(typeID),'-')
    plotLine(FR(S2,9:48),c(typeID),'--')
    plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,3:12))*4,'UniformOutput',false))...
        ,35*ones(4*nnz(p(typeID,3:12)),1),c(typeID),'Marker','.','LineStyle','none')   
    cellfun(@(x)fill([x,flip(x)],[0,0,40,40],[0.9500 0.8500 0.0030],'EdgeColor','none'),{[4.5,8.5]})
    cellfun(@(x)fill([x,flip(x)],[0,0,40,40],[0.4660 0.6740 0.1880],'EdgeColor','none'),{[16.5,20.5],[28.5,32.5]})
    alpha(0.2)
    set(gca,'XTick',4.5:20:40,'XTickLabel',{'0','5'})
    xlabel('time(s)')
    ylabel('FR(Hz)')
    
    num=10; % trial plot num
    b=10;
    subplot(2,3,3+typeID)
    hold on
    cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'EdgeColor','none'),{[0,1]})
    cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'EdgeColor','none'),{[3,4],[6,7]})
    alpha(0.2)
    for itr =1:num % trial
        if isempty(Raster1{itr+b,:})
            Raster1{itr+b,:}=-3;
        end
        if isempty(Raster2{itr+b,:})
            Raster2{itr+b,:}=-3;
        end
        plot([Raster1{itr+b,:};Raster1{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'r');
        hold on
        plot([Raster2{itr+b,:};Raster2{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'b');
        hold on
    end
    
    set(gca,'Xlim',[-1 9],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    xlabel('time(s)')
    ylabel('Trial ID')
end
exportgraphics(fh,fullfile(homedir,'showcase',sprintf('%d_%s.pdf',id,reg{cluster_id==id,7})),'ContentType','vector')


%% function
function [avail,out,FR,trials]=pre_process(homedir,rootpath,id)
sps=30000;
trials=markLPerf(h5read(fullfile(homedir,'DataSum',rootpath,'FR_All_250ms.hdf5'),'/Trials'));
if isempty(trials)
    avail=false;
    out=[];
    return
end
FR_All=h5read(fullfile(homedir,'DataSum',rootpath,'FR_All_250ms.hdf5'),'/FR_All');
SU_id=h5read(fullfile(homedir,'DataSum',rootpath,'FR_All_250ms.hdf5'),'/SU_id');

FR=FR_All(:,:,SU_id==id);
spkTS=h5read(fullfile(homedir,'DataSum',rootpath,'spike_times.hdf5'),'/spkTS');
spkId=readNPY(fullfile(homedir,'DataSum',rootpath,'spike_clusters.npy'));

FT_SPIKE=struct();
FT_SPIKE.label=strtrim(cellstr(num2str(id)));
FT_SPIKE.timestamp=cell(1,numel(id));
for i=1:numel(id)
    FT_SPIKE.timestamp{i}=spkTS(spkId==id)';
end
%  continuous format F T struct file
cfg=struct();
cfg.trl=[trials(:,1)-1*sps,trials(:,1)+9*sps,zeros(size(trials,1),1)-1*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end

function OBM1_index=isOBM1Map(reg)
load('D:/code-hem/OBM1Map.mat','OBM1map')
OBM1_index=[];
for i=1:length(reg)
    try
        OBM1_index{end+1,2}=OBM1map(reg{i});
        OBM1_index{end,1}=reg{i};
    catch     
    end
end
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

end

function plotLine(Data,c,l)
    m=mean(Data,1);
    plot(1:size(m,2),m,c,'LineStyle',l);
    ci=bootci(1000,@(x) mean(x),Data);
    fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.1);
end