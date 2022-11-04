%%
clc
clear
close all
addpath('D:/code/npy-matlab/npy-matlab')
addpath('D:/code/fieldtrip-20200320')
ft_defaults
homedir= 'F:/pixel-dualtask';

cluster_id=h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/cluster_id');
path=regexp(h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/path'),'(\w|\\|-)*','match','once');

trial_type=["distractorNo","distractor"];
load(fullfile(homedir,'regs_w2g.mat'),'reg_new');
reg=reg_new;
return
%% Quick and dirty
% id=cluster_id(any(sus_trans(:,1:2)~=0,2)&ismember(reg(:,7),OBM1_index(:,1)));
% load(fullfile(homedir,'cells_task.mat'));
% sust_trans=h5read(fullfile(homedir,'Selectivity_0606.hdf5'),'/sust_trans_noPermutaion');
% dist_sel=h5read(fullfile(homedir,'distractor_0606.hdf5'),'/sust_trans_noPermutaion');
% sust_trans(~(cells_task),1:2)=0;
% dist_sel(~(cells_task),1:2)=0;
% id=cluster_id(any(sust_trans(:,1:2),2)&ismember(reg(:,4),'CTX'));
load(fullfile(homedir,'difftype_ByPEV.mat'))
id=lost_all.id;


set(0,'DefaultFigureVisible', 'off') 
% figure('visible','off')
c=["k","r","b"];
for i=1:length(id)
    mem=util.mem_type(homedir,id(i));

    FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id(i)},'FR_All_1000ms.hdf5'),'/FR_All');
    SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id(i)},'FR_All_1000ms.hdf5'),'/SU_id');
    FR=h5read(fullfile(homedir,'DataSum',path{cluster_id==id(i)},'FR_All_250ms.hdf5'),'/FR_All');
    trial=util.markLPerf(h5read(fullfile(homedir,'DataSum',path{cluster_id==id(i)},'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
    for typeID=1:length(trial_type) 
%         [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials','InnerTask-correct');
        [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
        p(typeID,:)=arrayfun(@(x)ranksum(FR_All(S1,x,SU_id==rem(id(i),10000)),FR_All(S2,x,SU_id==rem(id(i),10000)))...
            ,1:size(FR_All,2))<0.05;         
    end
    if nnz(p(:,5:10))==0
       continue 
    end
    for typeID=1:length(trial_type) 
%         [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials','InnerTask-correct');
        [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
        hold on
        plotLine(FR(S1,:,SU_id==rem(id(i),10000)),c(typeID),'-')
        plotLine(FR(S2,:,SU_id==rem(id(i),10000)),c(typeID),'--')  
        plot(find(p(typeID,:))*4,1.2*max(mean((FR([S1;S2],:,SU_id==rem(id(i),10000))),1))*ones(nnz(p(typeID,:)),1),c(typeID),'Marker','*','LineStyle','none')
        
    end
    title(sprintf('%d %s, mem type=%d',id(i),reg{cluster_id==id(i),7},mem))
    saveas(gcf,fullfile(homedir,'showcase','lost',sprintf('%d_%s.png',id(i),reg{cluster_id==id(i),7})))
    close
end

%% sel-neuron showcase
close
id=1320174;
c={[0,0,0]/255,[112,112,112]/255;[255,0,0]/255,[255,153,153]/255;[0,0,255]/255,[153,153,255]/255};
trial_type=["distractorNo","distractor"];

[avail,FT_SPIKE,FR,trial]=pre_process(homedir,path{cluster_id==id},rem(id,10000));
for i=1:size(FT_SPIKE.trialtime,1)
    spike{i,:}=FT_SPIKE.time{1}(1,FT_SPIKE.trial{1,1}==i);
end
mem=util.mem_type(homedir,id);

FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/FR_All');
SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/SU_id');
for typeID=1:length(trial_type)      
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
    p(typeID,:)=arrayfun(@(x)ranksum(FR_All(S1,x,SU_id==rem(id,10000)),FR_All(S2,x,SU_id==rem(id,10000)))...
        ,1:size(FR_All,2))<0.05;
end
set(0,'DefaultFigureVisible', 'on') 
time_plot=[25,48];  %%%%%%
y_plot=15;
p_posi=14;
fh=figure('Color','w','Position',[100,100,500,250]);
sgtitle(sprintf('%d %s, mem type=%d',id,reg{cluster_id==id,7},mem))

for typeID=1:length(trial_type)     
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
    Raster1=spike(S1,1);
    Raster2=spike(S2,1);
    
    subplot(2,3,typeID)
    hold on
    plotLine(FR(S1,time_plot(1):time_plot(2)),c{typeID,1},'-')
    plotLine(FR(S2,time_plot(1):time_plot(2)),c{typeID,2},'-')
   
    if all(ismember(time_plot,[9,48]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,3:12))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,3:12)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[4.5,8.5]})
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[16.5,20.5],[28.5,32.5]})
        set(gca,'XTick',4.5:20:40,'XTickLabel',{'0','5'},'ylim',[0,y_plot])
    elseif all(ismember(time_plot,[41,56]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,11:14))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,11:14)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),8.5)
        set(gca,'XTick',4.5:4:12.5,'XTickLabel',-1:1:1)
    elseif all(ismember(time_plot,[25,48]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,7:12))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,7:12)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[0.5,4.5],[12.5,16.5]})
        set(gca,'XTick',0.5:20:20.5,'XTickLabel',{'0','5'},'Ylim',[0,y_plot])
    end

    
    xlabel('time(s)')
    ylabel('FR(Hz)')
    
    num=10; % trial plot num
    b=1;
    subplot(2,3,3+typeID)
    hold on
    if all(ismember(time_plot,[9,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[0,1]})
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[-1 9],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[9,10]})
        set(gca,'Xlim',[7 11],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[25,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[3 9],'Ylim',[0 20],'XTick',3:5:9,'XTickLabel',{'0','5'})
    end    
    for itr =1:num % trial
        if isempty(Raster1{itr+b,:})
            Raster1{itr+b,:}=-3;
        end
        if isempty(Raster2{itr+b,:})
            Raster2{itr+b,:}=-3;
        end
        plot([Raster1{itr+b,:};Raster1{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'color',c{typeID,1});        
        plot([Raster2{itr+b,:};Raster2{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'color',c{typeID,2});
       
    end

    xlabel('time(s)')
    ylabel('Trial ID')
end
exportgraphics(fh,fullfile(homedir,'showcase',sprintf('%d_%s.pdf',id,reg{cluster_id==id,7})),'ContentType','vector')
%% dist-neuron showcase
close
id=2100143;
c={[0,0,0]/255,[112,112,112]/255;[255,0,0]/255,[255,153,153]/255;[0,0,255]/255,[153,153,255]/255};


trial_type=["distractorGo","distractorNoGo"];

[avail,FT_SPIKE,FR,trial]=pre_process(homedir,path{cluster_id==id},rem(id,10000));
for i=1:size(FT_SPIKE.trialtime,1)
    spike{i,:}=FT_SPIKE.time{1}(1,FT_SPIKE.trial{1,1}==i);
end
mem=util.mem_type(homedir,id);

FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/FR_All');
SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/SU_id');

[go{1},go{2}]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(1)));
[nogo{1},nogo{2}]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(2)));

set(0,'DefaultFigureVisible', 'on') 
time_plot=[21,44];  %%%%%%
y_plot=40;
p_posi=38;
fh=figure('Color','w','Position',[100,100,333,250]);
sgtitle(sprintf('%d %s, mem type=%d',id,reg{cluster_id==id,7},mem))

for typeID=1:2
    p=arrayfun(@(x)ranksum(FR_All(go{typeID},x,SU_id==rem(id,10000)),FR_All(nogo{typeID},x,SU_id==rem(id,10000)))...
        ,1:size(FR_All,2))<0.05;
    
    Raster1=spike(go{typeID},1);
    Raster2=spike(nogo{typeID},1);
    
    subplot(2,2,typeID)
    hold on
    plotLine(FR(go{typeID},time_plot(1):time_plot(2)),c{2,1},'-')
    plotLine(FR(nogo{typeID},time_plot(1):time_plot(2)),c{3,1},'-')
   
    if all(ismember(time_plot,[9,48]))
         plot(cell2mat(arrayfun(@(x)x-3:x,find(p(:,3:12))*4,'UniformOutput',false))...
        ,11*ones(4*nnz(p(:,3:12)),1),c(:),'Marker','.','LineStyle','none')   
        cellfun(@(x)fill([x,flip(x)],[0,0,12,12],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[4.5,8.5]})
        cellfun(@(x)fill([x,flip(x)],[0,0,12,12],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[16.5,20.5],[28.5,32.5]})
        set(gca,'XTick',4.5:20:40,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(:,11:14))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(:,11:14)),1),'color',c{1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),8.5)
        set(gca,'XTick',4.5:4:12.5,'XTickLabel',-1:1:1)
    elseif all(ismember(time_plot,[21,44]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(:,6:11))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(:,6:11)),1)','color',c{1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),[4.5,16.5])
        set(gca,'XTick',4.5:8:20.5,'XTickLabel',0:2:8)
        
    end

    
    xlabel('time from distractor onset(s)')
    ylabel('FR(Hz)')
    
    num=10; % trial plot num
    b=1;
    subplot(2,2,2+typeID)
    hold on
    if all(ismember(time_plot,[9,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[0,1]})
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[-1 9],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[9,10]})
        set(gca,'Xlim',[7 11],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[21,44]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[2 8],'Ylim',[0 20],'XTick',3:2:12,'XTickLabel',0:2:8)
    end    
    for itr =1:num % trial
        if isempty(Raster1{itr+b,:})
            Raster1{itr+b,:}=-3;
        end
        if isempty(Raster2{itr+b,:})
            Raster2{itr+b,:}=-3;
        end
        plot([Raster1{itr+b,:};Raster1{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'color',c{2,1});        
        plot([Raster2{itr+b,:};Raster2{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'color',c{3,1});
       
    end

    xlabel('time from distractor onset (s)')
    ylabel('Trial ID')
end
exportgraphics(fh,fullfile(homedir,'showcase',sprintf('dist-%d_%s.pdf',id,reg{cluster_id==id,7})),'ContentType','vector')

%% mixed-neuron showcase
close
id=2100143;
% c={[0,0,0]/255,[112,112,112]/255;[255,0,0]/255,[255,153,153]/255;[0,0,255]/255,[153,153,255]/255};
c={'r','b';'r','b'};

[avail,FT_SPIKE,FR,trial]=pre_process(homedir,path{cluster_id==id},rem(id,10000));
for i=1:size(FT_SPIKE.trialtime,1)
    spike{i,:}=FT_SPIKE.time{1}(1,FT_SPIKE.trial{1,1}==i);
end
mem=util.mem_type(homedir,id);
FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/FR_All');
SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/SU_id');
set(0,'DefaultFigureVisible', 'on') 

fh=figure('Color','w','Position',[100,100,666,250]);
sgtitle(sprintf('%d %s, mem type=%d',id,reg{cluster_id==id,7},mem))

trial_type=["distractorNo","distractor"];
for typeID=1:length(trial_type)      
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
    p(typeID,:)=arrayfun(@(x)ranksum(FR_All(S1,x,SU_id==rem(id,10000)),FR_All(S2,x,SU_id==rem(id,10000)))...
        ,1:size(FR_All,2))<0.05;
end
time_plot=[9,48];  %%%%%%
y_plot=30;
p_posi=29;
for typeID=1:length(trial_type)     
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(typeID)));
    Raster1=spike(S1,1);
    Raster2=spike(S2,1);
    
    num=10; % trial plot num
    b=1;
    subplot(2,4,typeID)
    hold on
    if all(ismember(time_plot,[9,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[0,1]})
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[-1 9],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[9,10]})
        set(gca,'Xlim',[7 11],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[25,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[3 9],'Ylim',[0 20],'XTick',3:5:9,'XTickLabel',{'0','5'})
    end    
    for itr =1:num % trial
        if isempty(Raster1{itr+b,:})
            Raster1{itr+b,:}=-3;
        end
        if isempty(Raster2{itr+b,:})
            Raster2{itr+b,:}=-3;
        end
        plot([Raster1{itr+b,:};Raster1{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'color',c{typeID,1});        
        plot([Raster2{itr+b,:};Raster2{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'color',c{typeID,2});
       
    end

    xlabel('time(s)')
    ylabel('Trial ID')
    
    subplot(2,4,4+typeID)
    hold on
    plotLine(FR(S1,time_plot(1):time_plot(2)),c{typeID,1},'-')
    plotLine(FR(S2,time_plot(1):time_plot(2)),c{typeID,2},'-')
   
    if all(ismember(time_plot,[9,48]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,3:12))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,3:12)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[4.5,8.5]})
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[16.5,20.5],[28.5,32.5]})
        set(gca,'XTick',4.5:20:40,'XTickLabel',{'0','5'},'ylim',[0,y_plot])
    elseif all(ismember(time_plot,[41,56]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,11:14))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,11:14)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),8.5)
        set(gca,'XTick',4.5:4:12.5,'XTickLabel',-1:1:1)
    elseif all(ismember(time_plot,[25,48]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,7:12))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,7:12)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[0.5,4.5],[12.5,16.5]})
        set(gca,'XTick',0.5:20:20.5,'XTickLabel',{'0','5'},'Ylim',[0,y_plot])
    end

    
    xlabel('time(s)')
    ylabel('FR(Hz)')
        
   
end
trial_type=["distractorGo","distractorNoGo"];


[go{1},go{2}]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(1)));
[nogo{1},nogo{2}]=util.ExtractTrial(trial,'task','dualtask','trials',sprintf('%s-correct',trial_type(2)));

time_plot=[21,44];  
y_plot=30;
p_posi=29;

for typeID=1:2
    p=arrayfun(@(x)ranksum(FR_All(go{typeID},x,SU_id==rem(id,10000)),FR_All(nogo{typeID},x,SU_id==rem(id,10000)))...
        ,1:size(FR_All,2))<0.05;
    
    Raster1=spike(go{typeID},1);
    Raster2=spike(nogo{typeID},1);
     num=10; % trial plot num
    b=1;
    subplot(2,4,2+typeID)
    hold on
    if all(ismember(time_plot,[9,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[0,1]})
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[-1 9],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[9,10]})
        set(gca,'Xlim',[7 11],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[21,44]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[2 8],'Ylim',[0 20],'XTick',3:2:12,'XTickLabel',0:2:8)
    end    
    for itr =1:num % trial
        if isempty(Raster1{itr+b,:})
            Raster1{itr+b,:}=-3;
        end
        if isempty(Raster2{itr+b,:})
            Raster2{itr+b,:}=-3;
        end
        plot([Raster1{itr+b,:};Raster1{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'color',c{typeID,1});        
        plot([Raster2{itr+b,:};Raster2{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'color',c{typeID,2});
       
    end

    xlabel('time from distractor onset (s)')
    ylabel('Trial ID')
    
    subplot(2,4,6+typeID)
    hold on
    plotLine(FR(go{typeID},time_plot(1):time_plot(2)),c{typeID,1},'-')
    plotLine(FR(nogo{typeID},time_plot(1):time_plot(2)),c{typeID,2},'-')
   
    if all(ismember(time_plot,[9,48]))
         plot(cell2mat(arrayfun(@(x)x-3:x,find(p(:,3:12))*4,'UniformOutput',false))...
        ,11*ones(4*nnz(p(:,3:12)),1),c(:),'Marker','.','LineStyle','none')   
        cellfun(@(x)fill([x,flip(x)],[0,0,12,12],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[4.5,8.5]})
        cellfun(@(x)fill([x,flip(x)],[0,0,12,12],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[16.5,20.5],[28.5,32.5]})
        set(gca,'XTick',4.5:20:40,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(:,11:14))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(:,11:14)),1),'color',c{1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),8.5)
        set(gca,'XTick',4.5:4:12.5,'XTickLabel',-1:1:1)
    elseif all(ismember(time_plot,[21,44]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(:,6:11))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(:,6:11)),1)','color',c{1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),[4.5,16.5])
        set(gca,'XTick',4.5:8:20.5,'XTickLabel',0:2:8)
        
    end

    
    xlabel('time from distractor onset(s)')
    ylabel('FR(Hz)')
    
   
end
exportgraphics(fh,fullfile(homedir,'showcase',sprintf('%d_%s.pdf',id,reg{cluster_id==id,7})),'ContentType','vector')
%% reactivated neuron showcase
close
id=2700481;
c={[0,0,0]/255,[112,112,112]/255;[255,0,0]/255,[255,153,153]/255;[0,0,255]/255,[153,153,255]/255};
trial_type=["distractorNo-correct","distractor-correct","distractor-error"];

[avail,FT_SPIKE,FR,trial]=pre_process(homedir,path{cluster_id==id},rem(id,10000));
for i=1:size(FT_SPIKE.trialtime,1)
    spike{i,:}=FT_SPIKE.time{1}(1,FT_SPIKE.trial{1,1}==i);
end
mem=util.mem_type(homedir,id);

FR_All=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/FR_All');
SU_id=h5read(fullfile(homedir,'DataSum',path{cluster_id==id},'FR_All_1000ms.hdf5'),'/SU_id');
for typeID=1:length(trial_type)      
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',trial_type(typeID));
    p(typeID,:)=arrayfun(@(x)ranksum(FR_All(S1,x,SU_id==rem(id,10000)),FR_All(S2,x,SU_id==rem(id,10000)))...
        ,1:size(FR_All,2))<0.05;
end
set(0,'DefaultFigureVisible', 'on') 
time_plot=[41,56];  %%%%%%
y_plot=12;
p_posi=11;
fh=figure('Color','w','Position',[100,100,500,250]);
sgtitle(sprintf('%d %s, mem type=%d',id,reg{cluster_id==id,7},mem))

for typeID=1:length(trial_type)     
    [S1,S2]=util.ExtractTrial(trial,'task','dualtask','trials',trial_type(typeID));
    Raster1=spike(S1,1);
    Raster2=spike(S2,1);
    
    subplot(2,3,typeID+3)
    hold on
    plotLine(FR(S1,time_plot(1):time_plot(2)),c{typeID,1},'-')
    plotLine(FR(S2,time_plot(1):time_plot(2)),c{typeID,2},'-')
   
    if all(ismember(time_plot,[9,48]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,3:12))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,3:12)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[4.5,8.5]})
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[16.5,20.5],[28.5,32.5]})
        set(gca,'XTick',4.5:20:40,'XTickLabel',{'0','5'},'ylim',[0,y_plot])
    elseif all(ismember(time_plot,[41,56]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,11:14))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,11:14)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        arrayfun(@(x)fill([x,x+4,x+4,x],[-0.5,-0.5,y_plot,y_plot],'w','FaceColor',[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),8.5)
        set(gca,'XTick',4.5:4:12.5,'XTickLabel',-1:1:1)
    elseif all(ismember(time_plot,[25,48]))
        plot(cell2mat(arrayfun(@(x)x-3:x,find(p(typeID,7:12))*4,'UniformOutput',false))...
            ,p_posi*ones(4*nnz(p(typeID,7:12)),1),'color',c{typeID,1},'Marker','.','LineStyle','none')
        cellfun(@(x)fill([x,flip(x)],[0,0,y_plot,y_plot],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[0.5,4.5],[12.5,16.5]})
        set(gca,'XTick',0.5:20:20.5,'XTickLabel',{'0','5'},'Ylim',[0,y_plot])
    end

    
    xlabel('time(s)')
    ylabel('FR(Hz)')
    
    num=10; % trial plot num
    b=1;
    subplot(2,3,typeID)
    hold on
    if all(ismember(time_plot,[9,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[0,1]})
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[-1 9],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[41,56]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[9,10]})
        set(gca,'Xlim',[7 11],'Ylim',[0 20],'XTick',0:5:9,'XTickLabel',{'0','5'})
    elseif all(ismember(time_plot,[25,48]))
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
        set(gca,'Xlim',[3 9],'Ylim',[0 20],'XTick',3:5:9,'XTickLabel',{'0','5'})
    end    
    for itr =1:num % trial
        if isempty(Raster1{itr+b,:})
            Raster1{itr+b,:}=-3;
        end
        if isempty(Raster2{itr+b,:})
            Raster2{itr+b,:}=-3;
        end
        plot([Raster1{itr+b,:};Raster1{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'color',c{typeID,1});        
        plot([Raster2{itr+b,:};Raster2{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'color',c{typeID,2});
       
    end

    xlabel('time(s)')
    ylabel('Trial ID')
end
exportgraphics(fh,fullfile(homedir,'showcase',sprintf('%d_%s.pdf',id,reg{cluster_id==id,7})),'ContentType','vector')

%% function
function [avail,out,FR,trials]=pre_process(homedir,rootpath,id)
sps=30000;
trials=util.markLPerf(h5read(fullfile(homedir,'DataSum',rootpath,'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
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
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+14*sps,zeros(size(trials,1),1)-3*sps,trials];
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


function plotLine(Data,c,l)
    m=mean(Data,1);
    plot(1:size(m,2),m,'color',c,'LineStyle',l);
    ci=bootci(1000,{@(x) mean(x,1),Data},'type','normal');
    fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.1);
end