
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

%% state of all neurons
sus_trans=h5read(path.selectivity,'/sus_trans_noPermutaion')';
cluster_id=h5read(path.selectivity,'/cluster_id');
reg_temp=h5read(path.selectivity,'/reg')';
reg=regexp(reg_temp(:,7),'(\w|\\|-)*','match','once');

session_good=floor(cluster_id/100000);
sel_good = sus_trans(:,1)==1 | sus_trans(:,2)==1;
% sel_good = sus_trans(:,2)==1;

nonsel_good = sus_trans(:,1)==0 & sus_trans(:,2)==0 & sus_trans(:,3)==0 & sus_trans(:,4)==0;
reg_good=unique(reg);
reg_good(ismissing(regexp(reg_good,'[A-Z]','match','once')))=[];

load(path.performance)
learning=PerfList(:,2)==1;
WT=PerfList(:,3)==1;
phase=learning;

load(path.region)
reg_good_L=ismember(reg,reg_good) & learning;
reg_good_WT=ismember(reg,reg_W) & WT;

% load(fullfile(HomePath,'reg_ccfid_map'),'reg2tree')
% for i=1:length(reg_L)
%     reg_L(i,2:length(reg2tree(reg_L{i}))+1)=reg2tree(reg_L{i});   
% end
% p=h5read('D:\pixel-optogenetic\LaserModulation.hdf5','/rank');
p=h5read('D:\pixel-optogenetic\LaserModulation.hdf5','/anova');

% load('D:\pixel-optogenetic\FR_AIopto_0330.mat','P')
% p=cell2mat(cellfun(@(x)(x(2,:)),P(:,1),'UniformOutput',false));
laser_good=any(p(:,17:20),2); %5~5.5~6~6.5



%% find laser-modulated neurons in different sessions and regions
%%%session 
if true
for i=1:length(session)    
    suidx_s{i,1}=length(cluster_id((session_good==i) & phase & sel_good));
    suidx_s{i,2}=length(cluster_id((session_good==i) & phase & nonsel_good));
    for bin=1:4        
        suidx_s{i,3}(bin)=length(cluster_id((session_good==i) & p(:,bin+16) & phase & sel_good));
        suidx_s{i,4}(bin)=length(cluster_id((session_good==i) & p(:,bin+16) & phase & nonsel_good));
    end   
    for s=1:2
        r{s}(i,:)= suidx_s{i,s+2}./suidx_s{i,s};
    end
end
%plot
% fh=figure('Color','w','Position',[100,100,250,250]);
% m=mean(ratio_s(~isnan(ratio_s(:,1)),:),1);
% s=std(ratio_s(~isnan(ratio_s(:,1)),:),1)/sqrt(size(ratio_s(~isnan(ratio_s(:,1)),:),1));
% hold on
% bar(1:length(m),m*100,'w')
% errorbar(1:length(m),m*100,s*100,'Color','k','LineStyle','none')
% t={'All neuron','Selective neuron','Non-selective'};
% ylabel('% laser-modulated neurons')
% set(gca,'XTick',1:length(m),'XTicklabel',t(2:3))
% print laser-effect.eps -depsc2 -r300
% close

fh=figure('Color','w','Position',[100,100,250,250]);
color={'k-','k--'};
for i=1:2    
    m=mean(r{i}(~isnan(r{i}(:,1)),:),1);
    s=std(r{i}(~isnan(r{i}(:,1)),:),1)/sqrt(size(r{i}(~isnan(r{i}(:,1)),:),1));
    hold on
    errorbar(1:length(m),m*100,s*100,color{i})
end
xlim([0,5])
ylim([0 35])
xlabel('time from laser onset (s)')
ylabel('% laser-modulated neurons')
set(gca,'XTick',1:length(m),'XTicklabel',{'0','0.5','1','1.5','2'})
legend('memory neuron','non-memory neuron')
for i=1:4
ppp(i,:)=ranksum(r{1}(~isnan(r{1}(:,1)),i),r{2}(~isnan(r{1}(:,1)),i));
end

print(fullfile(HomePath,'anova-laser-effect-crosstime.png'),'-dpng','-r300')
close
end

%%% CTX
if false
Catelogy=unique(reg_L(:,5));
for i=1:length(Catelogy)
    suidx_s{i,1}=cluster_id(ismember(reg,reg_L(ismember(reg_L(:,5),Catelogy{i}),1)) & reg_good_L & learning & sel_good);
    suidx_s{i,2}=cluster_id(ismember(reg,reg_L(ismember(reg_L(:,5),Catelogy{i}),1)) & reg_good_L & learning & nonsel_good);
    suidx_s{i,3}=cluster_id(ismember(reg,reg_L(ismember(reg_L(:,5),Catelogy{i}),1)) & reg_good_L & laser_good & learning & sel_good);
    suidx_s{i,4}=cluster_id(ismember(reg,reg_L(ismember(reg_L(:,5),Catelogy{i}),1)) & reg_good_L & laser_good & learning & nonsel_good);
    for s=1:2
        ratio_r(i,s)= length(suidx_s{i,s+2})/length(suidx_s{i,s});
    end
end

[ratio_r_sorted,index]=sort(ratio_r,'descend');
reg_sorted=Catelogy(index);

%plot
fh=figure('Color','w','Position',[100,100,600,600]);
for i=1:2
% m=mean(ratio_r_sorted(~isnan(ratio_r_sorted(:,i)),i),1);
% s=std(ratio_s(~isnan(ratio_s(:,1)),:),1)/sqrt(size(ratio_s(~isnan(ratio_s(:,1)),:),1));
r=ratio_r_sorted(~isnan(ratio_r_sorted(:,i))& ratio_r_sorted(:,i)~=0 ,i);
subplot(2,1,i)
hold on
bar(1:size(r,1),r,'w')
ylabel('% laser-modulated neurons')
set(gca,'XTick',1:size(r,1),'XTicklabel',reg_sorted(:,i),'XTickLabelRotation',45)
t={'Selective neuron','Non-selective'};
title(t{i})
end
print laser-effect-per-region.eps -depsc2 -r300
close
end

if false
for i=1:length(reg_L)
    suidx_s{i,1}=cluster_id(strcmp(reg,reg_L{i}) & learning & phase);
    suidx_s{i,2}=cluster_id(strcmp(reg,reg_L{i}) & learning & sel_good & phase);
    suidx_s{i,3}=cluster_id(strcmp(reg,reg_L{i}) & learning & nonsel_good & phase);
    suidx_s{i,4}=cluster_id(strcmp(reg,reg_L{i}) & laser_good & learning & phase);
    suidx_s{i,5}=cluster_id(strcmp(reg,reg_L{i}) & laser_good & learning & sel_good & phase);
    suidx_s{i,6}=cluster_id(strcmp(reg,reg_L{i}) & laser_good & learning & nonsel_good & phase);
    for s=1:3
        ratio_r(i,s)= length(suidx_s{i,s+3})/length(suidx_s{i,s});
    end
end


[ratio_r_sorted,index]=sort(ratio_r,'descend');
reg_sorted=reg_L(index);

%plot
fh=figure('Color','w','Position',[100,100,1000,800]);
for i=1:3
% m=mean(ratio_r_sorted(~isnan(ratio_r_sorted(:,i)),i),1);
% s=std(ratio_s(~isnan(ratio_s(:,1)),:),1)/sqrt(size(ratio_s(~isnan(ratio_s(:,1)),:),1));
r=ratio_r_sorted(~isnan(ratio_r_sorted(:,i))& ratio_r_sorted(:,i)~=0 ,i);
subplot(3,1,i)
hold on
bar(1:size(r,1),r,'w')
ylabel('% laser-modulated neurons')
set(gca,'XTick',1:size(r,1),'XTicklabel',reg_sorted(:,i),'XTickLabelRotation',45)
t={'All neuron','Selective neuron','Non-selective'};
title(t{i})
end
% print laser-effect-per-region.eps -depsc2 -r300
% close


fh=figure('Color','w','Position',[100,100,1000,500]);

for i=1:2
    subplot(1,2,i)
projection_density=h5read(fullfile(HomePath,'AllenData_Projection_20210413.hdf5'),'/projection_density');
[r,p]=corrcoef(mean(projection_density(~isnan(ratio_r(:,1)),:),2),ratio_r(~isnan(ratio_r(:,1)),i+1))
hold on
plot(mean(projection_density,2),ratio_r(:,i+1),'k.','LineStyle','none','MarkerSize',10)
text(mean(projection_density,2),ratio_r(:,i+1),reg_L,'%')
% text(0.5,0.6,num2str(r))
xlabel('projection density (data from AllenInstitute)','FontSize',12);
ylabel('fraction of laser-modulated neruons','FontSize',12);
box off
end
exportgraphics(fh,fullfile(HomePath,'laser-projection-correlation.pdf'),'ContentType','vector')

print laser-projection-correlation.png -dpng -r300

end
%% plot singlecase of laser-modulated neurons 
clc
clear
close all

if ispc
    CodePath='D:\code-hem';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
path=util.Path_default;

p_rank=h5read('D:\pixel-optogenetic\LaserModulation.hdf5','/rank');
p_anova=h5read('D:\pixel-optogenetic\LaserModulation.hdf5','/anova');
cluster_id=h5read(path.selectivity,'/cluster_id');
path=regexp(h5read(path.selectivity,'/path'),'(\w|\\|-)*','match','once');
id_rank=cluster_id(any(p_rank(:,17:20),2));
id_anova=cluster_id(any(p_anova(:,17:20),2));
id_rank1=id_rank(~ismember(id_rank,intersect(id_rank,id_anova)));
id_anova1=id_rank(~ismember(id_anova,intersect(id_rank,id_anova)));
if true    
%learning+sel
id=id_rank1;
% id=cluster_id(reg_good_L & learning & p(:,17)<0.05);
for i=1:length(id)    
    plotSU(id(i),path{cluster_id==id(i)},HomePath)
    close all
end
end

%% Heatmap
id=cluster_id(reg_good_L & learning  & p(:,17)<0.05);

function plotSU(idx,path,HomePath)
FR_All=h5read(fullfile(HomePath,'DataSum',path,'FR_All.hdf5'),'/FR_All');
SU_id=h5read(fullfile(HomePath,'DataSum',path,'FR_All.hdf5'),'/SU_id');
trials=h5read(fullfile(HomePath,'DataSum',path,'FR_All.hdf5'),'/Trials');
cluster_id=rem(idx,10000);
FR{1}=FR_All(trials(:,9)==-1&trials(:,5)==4,:,SU_id==cluster_id);
FR{2}=FR_All(trials(:,9)==2&trials(:,5)==4,:,SU_id==cluster_id);
FR{3}=FR_All(trials(:,9)==2&trials(:,5)==8,:,SU_id==cluster_id);
FR{4}=FR_All(trials(:,9)==-1&trials(:,5)==8,:,SU_id==cluster_id);
Sel(1,:)=abs(mean(FR{1},1)-mean(FR{4},1))./(mean(FR{1},1)+mean(FR{4},1));
Sel(2,:)=abs(mean(FR{2},1)-mean(FR{3},1))./(mean(FR{2},1)+mean(FR{3},1));

%plot
c=[1,0,0;1,0.5,0.5;0.5,0.5,1;0,0,1];
fh=figure('visible','off','Color','w','Position',[100,100,1000,500]);
subplot(1,2,1)
for i=1:4
    hold on
    time=size(FR{i},2);
%     [CI, ~] = bootci(2000, @mean, FR{i});
%     value = [CI(2,:), CI(1,:)];   
    Time = [1:time, fliplr(1:time)];
    Highervalue = mean(FR{i},1) + std(FR{i},0,1)/sqrt(size(FR{i},1));
    Lowervalue = mean(FR{i},1) - std(FR{i},0,1)/sqrt(size(FR{i},1));
    value = [Highervalue, fliplr(Lowervalue)];
    a = fill(Time, value, c(i,:), 'edgecolor','none');
    alpha(a,0.2);    
    plot(1:time,mean(FR{i},1),'Color',c(i,:))
    Max(i)=max(mean(FR{i},1));
    Min(i)=min(mean(FR{i},1));
end

plot([12.5 12.5],[min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)],'k--')
plot([16.5 16.5],[min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)],'k--')
plot([32.5 32.5],[min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)],'k--')
plot([34.5 34.5],[min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)],'k--')
plot([40.5 40.5],[min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)],'k--')
plot([44.5 44.5],[min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)],'k--')
ylim([min(Min)-0.4*min(Min) max(Max)+0.2*max(Max)])
legend('S1-laser off','S1-laser on','S2-laser on','S2-laser off')

subplot(1,2,2)
hold on
plot(Sel(1,:),'k-')
plot(Sel(2,:),'b-')
plot([12.5 12.5],[min(Sel(:)) max(Sel(:))*1.2],'k--')
plot([16.5 16.5],[min(Sel(:)) max(Sel(:))*1.2],'k--')
plot([32.5 32.5],[min(Sel(:)) max(Sel(:))*1.2],'k--')
plot([34.5 34.5],[min(Sel(:)) max(Sel(:))*1.2],'k--')
plot([40.5 40.5],[min(Sel(:)) max(Sel(:))*1.2],'k--')
plot([44.5 44.5],[min(Sel(:)) max(Sel(:))*1.2],'k--')
ylim([min(Sel(:)) max(Sel(:))*1.2])
legend('laser off','laser on')

saveas(fh,fullfile(HomePath,'showcase','laser0526',num2str(idx)),'png')

end
