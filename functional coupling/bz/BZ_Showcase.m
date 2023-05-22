%%%% 
clearvars -Except sum 
clc

if ispc
    CodePath='D:\code-hem';
    HomePath='F:\pixel-dualtask';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end
addpath(fullfile('D:\code','npy-matlab','npy-matlab'))
addpath(fullfile('D:\code','fieldtrip-20200320'))
addpath(CodePath)
ft_defaults
%% quick and dirty
%%% FC showcase
% id=unique(floor(h5read(fullfile(HomePath,'Selectivity_1108.hdf5'),'/cluster_id')/100000));
% id=10;
% bzdata='10ms';
% pair_type='ext';
% set(0,'DefaultFigureVisible', 'off') 
% for fidx=id' 
%     disp(fidx)
%     BZ_SC_session(HomePath,bzdata,pair_type,fidx)
% end

%%% FCSP showcase
if ~(exist('sum','var')||exist('conn','var'))
    f=dir(fullfile(HomePath,'xcorr','coding_10ms','ext','fc_decoding_f*.mat'));
    sum=[];
    for i=1:size(f,1)
        fstr=load(fullfile(f(i).folder,f(i).name));
        sum=cat(1,sum,fstr.sums0);
    end
    sum(:,4)=mat2cell(util.mem_type(HomePath,cell2mat(sum(:,2))),ones(1,size(sum,1)),2);
    sum(:,5)=mat2cell(util.IDtoReg(HomePath,cell2mat(sum(:,2))),ones(1,size(sum,1)),2);
end

load(fullfile(HomePath,'difftype_ByPEV_1019.mat')')
% % NeuronSet_nonsel=util.ChooseNeuronSet_dualtask('sel','nonsel');
% % NeuronSet_sel=util.ChooseNeuronSet_dualtask('sel','sel');
% % NeuronSet_all=util.ChooseNeuronSet_dualtask('sel','all');
% % NeuronSet_dist=util.ChooseNeuronSet_dualtask('sel','dist');
% % load(fullfile(homedir,'dist_neuron.mat'),'NeuronSet_dist_act','NeuronSet_dist_coding')
% 
% id1=lost_all.id;
% id2=lost_all.id;
% conn=sum(cellfun(@(x)(ismember(x(1),id1)&ismember(x(2),id2))|(ismember(x(1),id2)&ismember(x(2),id1)),sum(:,2)),:);  
% FCSP_dirty(conn)

%%% FCSP coding
% id1=gain_all.id;
% id2=gain_all.id;
% conn=sum(cellfun(@(x)(ismember(x(1),id1)&ismember(x(2),id2))|(ismember(x(1),id2)&ismember(x(2),id1)),sum(:,2)),:);  
% FCSP_dirty(conn,'true')
%% single pair
% bzdata='10ms';
% pair_type='ext';
% cell1=1010202;
% cell2=1010212;
% BZ_SC(HomePath,bzdata,pair_type,cell1,cell2,)



%% FCSP showcase
% id1=2800137;
% id2=2800217;
% fh=FCSP_showcase(HomePath,sum(cellfun(@(x)ismember(x(1),id1)&ismember(x(2),id2),sum(:,2)),:));
% exportgraphics(fh,fullfile('F:\pixel-dualtask\xcorr\showcase\FCSP'...
%     ,sprintf('%d-%d.pdf',id1,id2)),'ContentType','vector')

%% FCSP coding showcase
id1=2800398;
id2=2800830;
fh=FCSP_showcase(HomePath,sum(cellfun(@(x)ismember(x(1),id1)&ismember(x(2),id2),sum(:,2)),:));
exportgraphics(fh,fullfile('F:\pixel-dualtask\xcorr\showcase\FCSP'...
    ,sprintf('%d-%d.pdf',id1,id2)),'ContentType','vector')
%% function
function BZ_SC_session(homedir,bzdata,pair_type,fidx)
load(fullfile(homedir,'xcorr',sprintf('bzdata_%s',bzdata),'bzdata_ZX',sprintf('%s_conn_w_reg_%d.mat',pair_type,fidx)),'sig_meta','pc_stem')

conn=double(sig_meta.suid(all((sig_meta.reg(:,1,:)==567|sig_meta.reg(:,1,:)==343)&sig_meta.reg(:,5,1)~=sig_meta.reg(:,5,2),3),:));
conn_reg=double(sig_meta.reg(all((sig_meta.reg(:,1,:)==567|sig_meta.reg(:,1,:)==343)&sig_meta.reg(:,5,1)~=sig_meta.reg(:,5,2),3),:,:));
SU_id=unique(conn);
[avail,FT_SPIKE,trials]=pre_process(fullfile(homedir,'DataSum',pc_stem),fidx,SU_id);    
idmap=load('reg_ccfid_map.mat');
ts_id_tagged=cell(0);
if strcmp(pair_type,'ext')
    load(fullfile(homedir,'xcorr',sprintf('bzdata_%s',bzdata),sprintf('BZ_XCORR_duo_f%d.mat',fidx)),'mono')
else
    load(fullfile(homedir,'xcorr',sprintf('bzdata_%s',bzdata),sprintf('BZ_XCORR_inh_f%d.mat',fidx)),'mono')
end
for i=1:size(conn,1)    
    for trlIdx=1:size(trials,1)
        ts1=FT_SPIKE.time{SU_id==conn(i,1)}(FT_SPIKE.trial{SU_id==conn(i,1)}==trlIdx)';
        ts2=FT_SPIKE.time{SU_id==conn(i,2)}(FT_SPIKE.trial{SU_id==conn(i,2)}==trlIdx)';
        ts1(:,2)=1;
        ts2(:,2)=2;
        ts_id=[ts1;ts2];
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged{trlIdx}=fc_tag(ts_id); %%???find spike pairs (>2ms, <10ms)
        clear ts1 ts2
    end
    trial_good=find(cellfun(@(x)any(x(:,3)==1),ts_id_tagged));
    if nnz(trial_good)<10
        continue
    end    
    CCH=squeeze(mono.ccgR(:,mono.completeIndex(mono.completeIndex(:,2)==conn(i,1),3),mono.completeIndex(mono.completeIndex(:,2)==conn(i,2),3)));   
    [M,I]=findpeaks(CCH,'MinPeakHeight',median(CCH)*1.5,'SortStr','descend','NPeaks',2);
    if ~isempty(I)&&I(1)<(251+25)&&I(1)>(251-25)
        if length(M)==1||(length(M)>1&&M(1)>M(2)*1.5)
            hold on
            plot(1:501,CCH,'r-');
            arrayfun(@(x)plot(repmat(x,1,2),[min(CCH),M(1)],'k--'),[251-25,251,251+25])
            title(sprintf('%d-%d,%s-%s',conn(i,1),conn(i,2),cell2mat(idmap.ccfid2reg(conn_reg(i,5,1))),cell2mat(idmap.ccfid2reg(conn_reg(i,5,2)))))
            set(gca,'XLim',[251-50,251+75])
            saveas(gcf,fullfile(homedir,'xcorr','showcase',sprintf('s%d',fidx)...
                ,sprintf('%d-%d-%s-%s.png',conn(i,1),conn(i,2),cell2mat(idmap.ccfid2reg(conn_reg(i,5,1))),cell2mat(idmap.ccfid2reg(conn_reg(i,5,2))))))
            close
        end
    end    
end
end

function BZ_SC(homedir,bzdata,pair_type,cell1,cell2)
fidx=floor(cell1/100000);
load(fullfile(homedir,'xcorr',sprintf('bzdata_%s',bzdata),'bzdata_ZX',sprintf('%s_conn_w_reg_%d.mat',pair_type,fidx)),'sig_meta','pc_stem')


fh=figure('Color','w','Position',[50,50,500,400]);
% sgtitle(sprintf('%d-%d,%s-%s',cell1,cell2,reg{cluster_id==cell1,7},reg{cluster_id==cell2,7}))

%% spike event raster
[avail,FT_SPIKE,trials]=pre_process(fullfile(homedir,'DataSum',pc_stem),fidx,[cell1;cell2],[1,2]);    
ts_id_tagged=cell(0);
for trlIdx=1:size(trials,1)
    ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
    ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';
    ts1(:,2)=1;
    ts2(:,2)=2;
    ts_id=[ts1;ts2];    
    [~,s]=sort(ts_id(:,1));
    ts_id=ts_id(s,:);
    ts_id_tagged{trlIdx}=fc_tag(ts_id); %%???find spike pairs (>2ms, <10ms)    
    clear ts1 ts2
end
trial_good=find(cellfun(@(x)any(x(:,3)==1),ts_id_tagged));  
subplot(2,2,[1,2])
num=4; % trial plot num
b=65;
hold on
% cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'EdgeColor','none'),{[0,1]})
% cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'EdgeColor','none'),{[3,4],[6,7]})
alpha(0.2)
for itr =1:num % trial
    spike1=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==1,1)';
    spike2=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==2,1)';
    spike1_fc=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==1&ts_id_tagged{itr+b}(:,3)==1,1)';
    spike2_fc=ts_id_tagged{itr+b}(ts_id_tagged{itr+b}(:,2)==2&ts_id_tagged{itr+b}(:,3)==1,1)';
    if isempty(spike1)
        spike1=-3;
    end
    if isempty(spike2)
        spike2=-3;
    end
    if isempty(spike1_fc)
        spike1_fc=-3;
    end
    if isempty(spike2_fc)
        spike2_fc=-3;
    end
    plot([spike1;spike1], [itr*2-1-0.25 itr*2-1+0.25],'color',[1,0.6,0.6]);
    plot([spike2;spike2], [itr*2-0.25 itr*2+0.25],'color',[0.6,0.6,1]);
    alpha(0.2)
    plot([spike1_fc;spike1_fc], [itr*2-1-0.25 itr*2-1+0.25],'r');
    plot([spike2_fc;spike2_fc], [itr*2-0.25 itr*2+0.25],'b');
    set(gca,'XLim',[1,2])
end
%% CCH
if strcmp(pair_type,'ext')
    load(fullfile(homedir,'xcorr',sprintf('bzdata_%s',bzdata),sprintf('BZ_XCORR_duo_f%d.mat',fidx)),'mono')
else
    load(fullfile(homedir,'xcorr',sprintf('bzdata_%s',bzdata),sprintf('BZ_XCORR_inh_f%d.mat',fidx)),'mono')
end
CCH=squeeze(mono.ccgR(:,mono.completeIndex(mono.completeIndex(:,2)==cell1,3),mono.completeIndex(mono.completeIndex(:,2)==cell2,3)));
subplot(2,2,3)
hold on
plot(1:501,CCH,'r-');
arrayfun(@(x)plot(repmat(x,1,2),[0,1500],'k--'),[251-25,251,251+25])
set(gca,'XLim',[251-50,251+75],'XTick',251-50:25:251+75,'YLim',[0,600])
xlabel('Latency (ms)')
ylabel('Splike count')

%% information
subplot(2,2,4)
hold on
idmap=load('reg_ccfid_map.mat');
text(0.1,0.1,idmap.ccfid2reg(sig_meta.reg(sig_meta.suid(:,1)==cell1&sig_meta.suid(:,2)==cell2,5,1)))
text(0.3,0.1,idmap.ccfid2reg(sig_meta.reg(sig_meta.suid(:,1)==cell1&sig_meta.suid(:,2)==cell2,5,2)))
text(0.1,0.3,num2str(cell1))
text(0.5,0.3,num2str(cell2))

exportgraphics(fh,fullfile(homedir,'xcorr','showcase',sprintf('%d_%d_%s_%s.pdf',cell1,cell2,bzdata,pair_type)),'ContentType','vector')

end

function FCSP_dirty(conn,IsSample)
set(0,'DefaultFigureVisible', 'off')

if IsSample
    for t=[1,7]
        FCSP_S1{1+(t-1)/6,1}=cell2mat(cellfun(@(x)permute(mean(x{2*t-1}(1,:,:),2),[2,3,1]),conn(:,3),'UniformOutput',false));        
        FCSP_S2{1+(t-1)/6,1}=cell2mat(cellfun(@(x)permute(mean(x{2*t}(1,:,:),2),[2,3,1]),conn(:,3),'UniformOutput',false));        
    end
else
    for t=[1,7]
        FCSP_sel{1+(t-1)/6,1}=cell2mat(cellfun(@(x)permute(mean(cat(2,x{2*t-1}(1,:,:),x{2*t}(1,:,:)),2),[2,3,1]),conn(:,3),'UniformOutput',false));        
    end
end

for i=1:size(conn,1)
    if ~all(cellfun(@(x)strcmp(x(3),'CH'),conn{i,5}))
        continue
    end
    fh=figure('Color','w','Position',[100,100,180,150]);
    
    hold on
    if IsSample
        plot(5:12,FCSP_S1{1}(i,5:12),'k.-')
        plot(5:12,FCSP_S1{2}(i,5:12),'r.-')       
        plot(5:12,FCSP_S2{1}(i,5:12),'k.--')
        plot(5:12,FCSP_S2{2}(i,5:12),'r.--') 
    else
        plot(5:12,FCSP_sel{1}(i,5:12),'k.-')
        plot(5:12,FCSP_sel{2}(i,5:12),'r.-')
    end
    set(gca,'XTick',[6.5,7.5,9.5,10.5])
    title(sprintf('%d-%d(%s-%s)',conn{i,2}(1),conn{i,2}(2),conn{i,5}{1}{7},conn{i,5}{2}{7}))
    saveas(fh,fullfile('F:\pixel-dualtask\xcorr\showcase\FCSP\gain_gain'...
        ,sprintf('%d-%d(%s-%s)-sample.png',conn{i,2}(1),conn{i,2}(2),conn{i,5}{1}{7},conn{i,5}{2}{7})))
end

set(0,'DefaultFigureVisible', 'on') 

end

function fh=FCSP_showcase(homedir,conn)
fh=figure('Color','w','Position',[50,50,150,180]);
sgtitle(sprintf('%d-%d(%s-%s)',conn{1,2}(1),conn{1,2}(2),conn{1,5}{1}{7},conn{1,5}{2}{7}),'FontSize',6)

fidx=floor(conn{2}(1)/100000);
load(fullfile(homedir,'xcorr',sprintf('bzdata_%s','10ms'),'bzdata_ZX',sprintf('%s_conn_w_reg_%d.mat','ext',fidx)),'pc_stem')

%%% spike event raster
[avail,FT_SPIKE,trials]=pre_process(fullfile(homedir,'DataSum',pc_stem),fidx,[conn{2}(1);conn{2}(2)],[3,9]);    
ts_id_tagged=cell(0);
for trlIdx=1:size(trials,1)
    ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
    ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';
    ts1(:,2)=1;
    ts2(:,2)=2;
    ts_id=[ts1;ts2];    
    [~,s]=sort(ts_id(:,1));
    ts_id=ts_id(s,:);
    ts_id_tagged{trlIdx}=fc_tag(ts_id); %%???find spike pairs (>2ms, <10ms)    
    clear ts1 ts2
end
trial{1}=find(trials(:,9)==-1&all(trials(:,14:15),2)&cellfun(@(x)any(x(:,3)==1),ts_id_tagged)');
trial{2}=find(trials(:,9)~=-1&all(trials(:,14:15),2)&cellfun(@(x)any(x(:,3)==1),ts_id_tagged)');
subplot(2,1,1)
color=[0,0,0;153,153,153;255,0,0;255,153,153]/255;
hold on
cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[3,4],[6,7]})
b=6;
num=4; % trial plot num
alpha(0.2)
for t=1:2
    for itr =1:num % trial
        spike1_fc=ts_id_tagged{trial{t}(itr+b)}(ts_id_tagged{trial{t}(itr+b)}(:,2)==1&ts_id_tagged{trial{t}(itr+b)}(:,3)==1,1)';
        spike2_fc=ts_id_tagged{trial{t}(itr+b)}(ts_id_tagged{trial{t}(itr+b)}(:,2)==2&ts_id_tagged{trial{t}(itr+b)}(:,3)==1,1)';
        
        if isempty(spike1_fc)
            spike1_fc=-3;
        end
        if isempty(spike2_fc)
            spike2_fc=-3;
        end
        plot([spike1_fc;spike1_fc], [8*(t-1)+itr*2-0.75-0.25 8*(t-1)+itr*2-0.75+0.25],'color',color(2*t-1,:));
        plot([spike2_fc;spike2_fc], [8*(t-1)+itr*2-0.25 8*(t-1)+itr*2+0.25],'color',color(2*t,:));
        cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.9500 0.8500 0.0030],'FaceAlpha',0.2,'EdgeColor','none'),{[0,1]})
        
        set(gca,'XLim',[3,9],'YLim',[0.5,16.5],'FontSize',6)
        xlabel('time','FontSize',6)
        ylabel('tiral','FontSize',6)
    end
end

for t=[1,7]
    FCSP_pertrial{1+(t-1)/6,1}=permute([conn{1,3}{2*t-1}(1,:,:),conn{1,3}{2*t}(1,:,:)],[2,3,1]);
end
ci=cellfun(@(x)bootci(1000,{@(y)mean(y(:,7:12),1),x},'type','normal'),FCSP_pertrial,'UniformOutput',false);
FCSP_mean=cell2mat(cellfun(@(x)mean(x(:,7:12),1),FCSP_pertrial,'UniformOutput',false));
subplot(2,1,2)
hold on
cellfun(@(x)fill([x,flip(x)],[0,0,20,20],[0.4660 0.6740 0.1880],'FaceAlpha',0.2,'EdgeColor','none'),{[0.5,1.5],[3.5,4.5]})

plot(1:size(FCSP_mean,2),FCSP_mean(1,:),'k');
plot(1:size(FCSP_mean,2),FCSP_mean(2,:),'r');
errorbar(1:size(FCSP_mean,2),FCSP_mean(1,:),FCSP_mean(1,:)-ci{1}(1,:),ci{1}(2,:)-FCSP_mean(1,:),'k')
errorbar(1:size(FCSP_mean,2),FCSP_mean(2,:),FCSP_mean(2,:)-ci{2}(1,:),ci{2}(2,:)-FCSP_mean(2,:),'r')
set(gca,'XLim',[0.5,6.5],'XTick',0.5:1:6.5,'XTickLabel',0:1:6,'Ylim',[0,2],'FontSize',6)
xlabel('time','FontSize',6)
ylabel('FCSP rate','FontSize',6)
end

function [avail,out,trials]=pre_process(path,sessionIdx,cluster_ids,time)
sps=30000;


folderIdx=unique(floor((cluster_ids-100000*sessionIdx)/10000));
spkTS=[];
spkId=[];

for f=1:size(folderIdx,1)  
    folder=dir(fullfile(path,sprintf('*imec%d*',folderIdx(f)),'spike_clusters.npy'));
    spkId_temp=[];
    spkId_temp=readNPY(fullfile(folder.folder,'spike_clusters.npy'))+sessionIdx*100000;
    spkId=[spkId;spkId_temp+10000*str2num(regexp(folder.folder,'(?<=imec)(\d)','match','once'))];
    
    spkTS_temp=[];
    spkTS_temp=h5read(fullfile(folder.folder,'spike_times.hdf5'),'/spkTS');
    spkTS=[spkTS;spkTS_temp];
end


FT_SPIKE=struct();
FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
% for i=1:numel(cluster_ids)
%     FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
% end
FT_SPIKE.timestamp=arrayfun(@(x)spkTS(spkId==x)',cluster_ids,'UniformOutput',false);    
%  continuous format F T struct file
trials=util.markLPerf(h5read(fullfile(folder.folder,'FR_All_250ms.hdf5'),'/Trials'),0.7,0.8,120,'dualtask');
if isempty(trials)
    avail=false;
    out=[];
    return
end

cfg=struct();
cfg.trl=[trials(:,1)+time(1)*sps,trials(:,1)+time(2)*sps,zeros(size(trials,1),1)+time(1)*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end

function [out]=markLPerf(facSeq)
trial_min=120;
lick_criteria=0.7;
Perf_criteria=0.65;
i=trial_min;
facSeq_WT=zeros(length(facSeq),1);
while i<=length(facSeq)
    goodOff=nnz(xor(facSeq(i-trial_min+1:i,5)==facSeq(i-trial_min+1:i,6) , facSeq(i-trial_min+1:i,7)>0));
    lickOff=nnz(facSeq(i-trial_min+1:i,5)~=facSeq(i-trial_min+1:i,6)&facSeq(i-trial_min+1:i,7)>0);
    if goodOff>=Perf_criteria*trial_min && lickOff>=0.5*lick_criteria*trial_min %.60 correct rate
        facSeq_WT(i-trial_min+1:i,1)=1;
    end
    i=i+1;
end
out=[facSeq,facSeq_WT];

end


function out=isOBM1Map(reg)
load('D:/code-hem/OBM1Map.mat','OBM1map')

out=[];
for i=1:length(reg)
    try
        out{end+1,2}=OBM1map(reg{i});
        out{end,1}=reg{i};
    catch     
    end
end
end

function out=fc_tag(in)

out=in;
out(:,3)=0;

for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.01 && in(j,1)>in(i,1)+0.002
                out(i,3)=1;
                out(j,3)=1;
            end
            j=j+1;
        end
    end
end
end