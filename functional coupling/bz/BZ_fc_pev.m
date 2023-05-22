clearvars -Except sum
homedir='F:/pixel-dualtask';
reg=["all"];
% reg=["within-region","cross-region"];

% fc_reactivated=cell2mat(sum(cellfun(@(x)all(ismember(x,NeuronSet.id)),sum(:,3)),3));
%%
if ~exist('sum','var')
    [sig,~]=load_sig_pair(fullfile(homedir,'xcorr','bzdata_10ms','bzdata_zx'),'ext');
    sum=load_data(homedir,sig.suid);
    sum(:,5)=mat2cell(util.mem_type(homedir,cell2mat(sum(:,2))),ones(1,size(sum,1)),2);
    sum(:,6)=mat2cell(util.IDtoReg(homedir,cell2mat(sum(:,2))),ones(1,size(sum,1)),2);

end
% reactivated=load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet');

NeuronSet_nonsel=util.ChooseNeuronSet_dualtask('sel','nonsel');
NeuronSet_sel=util.ChooseNeuronSet_dualtask('sel','sel');
NeuronSet_all=util.ChooseNeuronSet_dualtask('sel','all');
load(fullfile(homedir,'difftype_ByPEV_1019.mat'))
load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet_gain');



id1{1}=NeuronSet_gain.id;
id2{1}=NeuronSet_gain.id;
id1{2}=setdiff(NeuronSet_all.id,NeuronSet_gain.id);
id2{2}=setdiff(NeuronSet_all.id,NeuronSet_gain.id);
id1{3}=NeuronSet_sel.id;
id2{3}=NeuronSet_gain.id;
id1{4}=NeuronSet_sel.id;
id2{4}=setdiff(NeuronSet_all.id,NeuronSet_gain.id);
id1{5}=NeuronSet_sel.id;
id2{5}=NeuronSet_sel.id;
id1{6}=NeuronSet_nonsel.id;
id2{6}=NeuronSet_nonsel.id;

conn=[];
for i=1:length(id1)
    sums=sum(cellfun(@(x)(ismember(x(1),id1{i})&ismember(x(2),id2{i}))|(ismember(x(1),id2{i})&ismember(x(2),id1{i})),sum(:,2)),:);    
    for regidx=reg
        switch(regidx)
            case 'within-region'
                conn_reg=cellfun(@(x)strcmp(x{1}{7},x{2}{7}),sums(:,6))&cellfun(@(x)any(strcmp(x{1}{3},{'CH','BS'}))&any(strcmp(x{2}{3},{'CH','BS'})),sums(:,6))...
                    &~cellfun(@(x)isempty(x{1}{7})|isempty(x{2}{7}),sums(:,6));
            case 'cross-region'
                conn_reg=~cellfun(@(x)strcmp(x{1}{7},x{2}{7}),sums(:,6))&cellfun(@(x)any(strcmp(x{1}{3},{'CH','BS'}))&any(strcmp(x{2}{3},{'CH','BS'})),sums(:,6))...
                    &~cellfun(@(x)isempty(x{1}{7})|isempty(x{2}{7}),sums(:,6));
            case 'all'
                conn_reg=true(size(sums,1),1);
        end
        conn{end+1}=sums(conn_reg,:);
    end
end
% sig1=structfun(@(x)x(ismember(sig.suid(:,1),id1)&ismember(sig.suid(:,2),id2)&graysel & sig.reg(:,5,1)==sig.reg(:,5,2),:,:),sig,'UniformOutput',false);
% sig2=structfun(@(x)x(ismember(sig.suid(:,1),id1)&ismember(sig.suid(:,2),id2)&graysel & sig.reg(:,5,1)~=sig.reg(:,5,2),:,:),sig,'UniformOutput',false);

%%
for i=1:length(conn)
%     pev{i}=fc_pev(conn{i},[8,9,12]);
    pev{i}=fc_pev(conn{i},[13]);
end

% readme{1}='conn 1-2: gain; 3-4: maintain; 5-6: lost; 7-8 nonmemory(with-region/cross-region)';
% readme{2}='conn{n}: row1-1: trials without/with inner task; col1-2: real/shuffle data';
readme{1}='conn 1-4(all region): reac, non-reac, mem<->reac, mem<->non-reac,mem, non-mem';
readme{2}='conn{n}: row1-1: trials without/with inner task; col1-2: real/shuffle data';
save(fullfile(homedir,'xcorr','pev','reac_all.mat'),'pev','reg','id1','id2','readme','-v7.3')
%

%%
load(fullfile(homedir,'xcorr','pev','reac_all.mat'),'pev')

pev_m=cell2mat(cellfun(@(x)nanmean(nanmean(x{2,1},2)-nanmean(x{1,1},2)),pev,'UniformOutput',false));
pev_s=cell2mat(cellfun(@(x)bootci(1000,{@(m)nanmean(m),nanmean(x{2,1},2)-nanmean(x{1,1},2)},'type','normal')...
    ,pev,'UniformOutput',false));




% pev_m=cell2mat(cellfun(@(y)cellfun(@(x)nanmean(nanmean(x,2)),y(:,1)),pev,'UniformOutput',false));
% pev_s=cell2mat(cellfun(@(y)cell2mat(cellfun(@(x)bootci(1000,{@(m)nanmean(m),nanmean(x,2)},'type','normal')...
%     ,y(:,1),'UniformOutput',false)),pev,'UniformOutput',false));
%%
color=[1,0,0;0.5,0.5,0.5;0,0,1;1,1,1];

f=figure('Color','w','Position',[100,100,205,125]);
subplot(1,2,1)
hold on
bh=bar(pev_m(:,1:2:8)','FaceColor','flat');
bh.CData=color;
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-pev_s(1,1:2:8),pev_s(2,1:2:8)-[bh.YEndPoints],'k.');
[p,~,st]=anovan(cell2mat(cellfun(@(x)nanmean(x{2,1},2)-nanmean(x{1,1},2),pev(1:2:8),'UniformOutput',false)')...
    ,cell2mat(arrayfun(@(z,y)z*ones(y,1),1:4,cellfun(@(x)size(x{1},1),pev(1:2:8)),'UniformOutput',false)'),'display','off');

disp(p)
c=multcompare(st,'Display','off')
set(gca,'XTick',2.5,'XTickLabel',{'within region'},'ylim',[-0.1,0.05],'FontSize',6)
ylabel('PEV(FCSP)')

% set(gca,'XTick',1:2,'XtickLabel',{'LD','Test'})
% set(gca,'XTick',1:2,'XtickLabel',reg,'YLim',[0,1.5])
% ylabel('pev for ODPA sample (FCSP)')
% title(sprintf('within region,%0.3f,%0.3f',p(1),p(2)))
subplot(1,2,2)
hold on
bh=bar(pev_m(:,2:2:8)','FaceColor','flat');
bh.CData=color;
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-pev_s(1,2:2:8),pev_s(2,2:2:8)-[bh.YEndPoints],'k.');
[p,~,st]=anovan(cell2mat(cellfun(@(x)nanmean(x{2,1},2)-nanmean(x{1,1},2),pev(2:2:8),'UniformOutput',false)')...
    ,cell2mat(arrayfun(@(z,y)z*ones(y,1),1:4,cellfun(@(x)size(x{1},1),pev(2:2:8)),'UniformOutput',false)'),'display','off');
disp(p)
c=multcompare(st,'Display','off')
set(gca,'XTick',2.5,'XTickLabel',{'cross region'},'ylim',[-0.1,0.05],'FontSize',6)
ylabel('PEV(FCSP)')
exportgraphics(f,fullfile(homedir,'xcorr','pev','diff_type_memory_within_group.pdf'));

% set(gca,'XTick',1:2,'XtickLabel',{'LD','Test'})
% set(gca,'XTick',1:2,'XtickLabel',reg,'YLim',[0,1.5])
% ylabel('pev for ODPA sample (FCSP)')
% title(sprintf('cross region,%0.3f,%0.3f',p(3),p(4)))

p=cellfun(@(y,z)ranksum(y,z),cellfun(@(x)nanmean(x{2,1},2)-nanmean(x{1,1},2),pev(1:2:8),'UniformOutput',false)...
    ,cellfun(@(x)nanmean(x{2,2},2)-nanmean(x{1,2},2),pev(1:2:8),'UniformOutput',false))
%%
color=[1,0,0;0.5,0.5,0.5;0,0,1;1,1,1];

f=figure('Color','w','Position',[100,100,307,125]);
subplot(1,3,1)
hold on
load(fullfile(homedir,'xcorr','pev','diff_type_memory_within_group_all.mat'),'pev')
pev_m=cell2mat(cellfun(@(x)nanmean(nanmean(x{2,1},2)-nanmean(x{1,1},2)),pev,'UniformOutput',false));
pev_s=cell2mat(cellfun(@(x)bootci(1000,{@(m)nanmean(m),nanmean(x{2,1},2)-nanmean(x{1,1},2)},'type','normal')...
    ,pev,'UniformOutput',false));
bh=bar(pev_m','FaceColor','flat');
bh.CData=color;
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-pev_s(1,:),pev_s(2,:)-[bh.YEndPoints],'k.');
[p,~,st]=anovan(cell2mat(cellfun(@(x)nanmean(x{2,1},2)-nanmean(x{1,1},2),pev,'UniformOutput',false)')...
    ,cell2mat(arrayfun(@(z,y)z*ones(y,1),1:4,cellfun(@(x)size(x{1},1),pev),'UniformOutput',false)'),'display','off');

disp(p)
c=multcompare(st,'Display','off')
set(gca,'XTick',2.5,'XTickLabel',{'within region'},'ylim',[-0.1,0.05],'FontSize',6)
ylabel('PEV(FCSP)')


load(fullfile(homedir,'xcorr','pev','diff_type_memory_cross_group_all.mat'),'pev')
pev_m=cell2mat(cellfun(@(x)nanmean(nanmean(x{2,1},2)-nanmean(x{1,1},2)),pev,'UniformOutput',false));
pev_s=cell2mat(cellfun(@(x)bootci(1000,{@(m)nanmean(m),nanmean(x{2,1},2)-nanmean(x{1,1},2)},'type','normal')...
    ,pev,'UniformOutput',false));

subplot(1,3,2)
hold on
bh=bar(pev_m(1:3)','FaceColor','flat');
bh.CData=color(1:3,:);
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-pev_s(1,1:3),pev_s(2,1:3)-[bh.YEndPoints],'k.');
[p,~,st]=anovan(cell2mat(cellfun(@(x)nanmean(x{2,1},2)-nanmean(x{1,1},2),pev(1:3),'UniformOutput',false)')...
    ,cell2mat(arrayfun(@(z,y)z*ones(y,1),1:3,cellfun(@(x)size(x{1},1),pev(1:3)),'UniformOutput',false)'),'display','off');

disp(p)
c=multcompare(st,'Display','off')
set(gca,'XTick',2.5,'XTickLabel',{'within region'},'ylim',[-0.02,0.01],'FontSize',6)
ylabel('PEV(FCSP)')

subplot(1,3,3)
hold on
bh=bar(pev_m([2,4,6])','FaceColor','flat');
bh.CData=color(1:3,:);
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-pev_s(1,[2,4,6]),pev_s(2,[2,4,6])-[bh.YEndPoints],'k.');
[p,~,st]=anovan(cell2mat(cellfun(@(x)nanmean(x{2,1},2)-nanmean(x{1,1},2),pev([2,4,6]),'UniformOutput',false)')...
    ,cell2mat(arrayfun(@(z,y)z*ones(y,1),1:3,cellfun(@(x)size(x{1},1),pev([2,4,6])),'UniformOutput',false)'),'display','off');

disp(p)
c=multcompare(st,'Display','off')
set(gca,'XTick',2.5,'XTickLabel',{'within region'},'ylim',[-0.04,0.02],'FontSize',6)
ylabel('PEV(FCSP)')
exportgraphics(f,fullfile(homedir,'xcorr','pev','diff_type_memory_all.pdf'));


%%
pev_m=cell2mat(cellfun(@(z)cellfun(@(x,y)nanmean(nanmean(x,2)-nanmean(y,2)),z(:,1),z(:,2)),pev,'UniformOutput',false));
pev_s=cell2mat(cellfun(@(z)cell2mat(cellfun(@(x,y)bootci(1000,{@(m)nanmean(m),nanmean(x,2)-nanmean(y,2)},'type','normal')...
    ,z(:,1),z(:,2),'UniformOutput',false)),pev,'UniformOutput',false));

% pev_s=cell2mat(cellfun(@(z)cellfun(@(x,y)bootci(1000,{@(m)nanmean(m),nanmean(x,2)-nanmean(y,2)},'type','normal'),z(:,1),z(:,2)),pev,'UniformOutput',false));

p=cellfun(@(a)ranksum(a{1},a{2}),cellfun(@(z)cellfun(@(x,y)nanmean(x,2)-nanmean(y,2),z(:,1),z(:,2),'UniformOutput',false),pev,'UniformOutput',false));
p=cellfun(@(a)anova1([a{1};a{2}],[zeros(length(a{1}),1);ones(length(a{2}),1)],'off'),cellfun(@(z)cellfun(@(x,y)nanmean(x,2)-nanmean(y,2),z(:,1),z(:,2),'UniformOutput',false),pev,'UniformOutput',false));

f=figure('Color','w','Position',[100,100,300,188]);
subplot(1,2,1)
hold on
bh=bar(pev_m');
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-reshape(pev_s([1,3],1:2)',1,[]),reshape(pev_s([2,4],1:2)',1,[])-[bh.YEndPoints],'k.');
set(gca,'XTick',1:2,'XtickLabel',{'LD','Test'})
% set(gca,'XTick',1:2,'XtickLabel',reg,'YLim',[0,1.5])
ylabel('pev for ODPA sample (FCSP)')
title(sprintf('within region,%0.3f,%0.3f',p(1),p(2)))
subplot(1,2,2)
hold on
bh=bar(pev_m(:,3:4)');
errorbar([bh.XEndPoints],[bh.YEndPoints],[bh.YEndPoints]-reshape(pev_s([1,3],3:4)',1,[]),reshape(pev_s([2,4],3:4)',1,[])-[bh.YEndPoints],'k.');

set(gca,'XTick',1:2,'XtickLabel',{'LD','Test'})
% set(gca,'XTick',1:2,'XtickLabel',reg,'YLim',[0,1.5])
ylabel('pev for ODPA sample (FCSP)')
title(sprintf('cross region,%0.3f,%0.3f',p(3),p(4)))
return

p=cellfun(@(a)ranksum(a{1},a{2}),cellfun(@(z)cellfun(@(x)nanmean(x,2),z(:,1),'UniformOutput',false),pev,'UniformOutput',false));

cell2mat(cellfun(@(z)cellfun(@(x,y)ranksum(nanmean(x,2),nanmean(y,2),'tail','right'),z(:,1),z(:,2)),pev,'UniformOutput',false))


cvcorr2=fc_dec(sum(cellfun(@(x)all(ismember(x(1),sig2.suid(:,1))&ismember(x(2),sig2.suid(:,2)),2),sum(:,2)),:));

%% function
function sum0=load_data(homedir,in)
% persistent sum0 
if ~exist('sum0','var')|| isempty(sum0)
%     f=unique(floor(in));
    f=dir(fullfile(homedir,'xcorr','coding_10ms','ext','fc_decoding_f*.mat'));

    sum0=[];
    for i=1:size(f,1)
        fstr=load(fullfile(homedir,'xcorr','coding_10ms','ext',f(i).name));
%         fstr=load(fullfile(homedir,'xcorr','coding_10ms','ext',sprintf('fc_decoding_f%d.mat',f(i))));
        sum0=cat(1,sum0,fstr.sums0);
        %     sum=cat(1,sum,fstr.sums0(cellfun(@(x)all(ismember(x(1),in(:,1))&ismember(x(2),in(:,2)),2),fstr.sums0(:,2)),:));
    end
    sum0(:,4)=mat2cell(util.mem_type(homedir,cell2mat(sum0(:,2))),ones(1,size(sum0,1)),2);
end
sum=sum0(cellfun(@(x)all(ismember(x(1),in(:,1))&ismember(x(2),in(:,2)),2),sum0(:,2)),:);
end

function out=fc_pev(in,bin)

NTRIAL=5;
NRPT=500;
pev=cell(0);
trial_type=["distractorNo-correct","distractorNoGo-correct","distractorGo-correct"...
    ,"distractorNo-error","distractorNoGo-error","distractorGo-error"...
    ,"distractor-correct","InnerTask-correct","InnerTask-correct-correct"...
    ,"distractor-error","InnerTask-error","InnerTask-error-correct"];
%%
for tridx=[1,7]
    clear centt stdd per_trial_norm
    for i=1:size(in,1)
%         centt(i,:)=mean(cell2mat(cellfun(@(x)permute((x(1,:,bin)),[3,2,1]),in{i,3}(2*tridx-1:2*tridx),'UniformOutput',false)),[1,2]);
%         stdd(i,:)=std(cell2mat(cellfun(@(x)permute(x(1,:,bin),[3,2,1]),in{i,3}(2*tridx-1:2*tridx),'UniformOutput',false)),0,[1,2]);
        per_trial_norm{i,1}=(mean(permute(in{i,3}{2*tridx-1}(1,:,bin),[3,2,1]),1))';
        per_trial_norm{i,2}=(mean(permute(in{i,3}{2*tridx}(1,:,bin),[3,2,1]),1))';
%         per_trial_norm{i,1}=(mean(permute(in{i,3}{2*tridx-1}(1,:,bin),[3,2,1]),1)-centt(i,:))'./stdd(i,:);
%         per_trial_norm{i,2}=(mean(permute(in{i,3}{2*tridx}(1,:,bin),[3,2,1]),1)-centt(i,:))'./stdd(i,:);
        
        %             % normal-baseline
        %             centt(i,:)=mean(cell2mat(cellfun(@(x)squeeze(x(1,:,1:2)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)'),'all');
        %             stdd(i,:)=std(cell2mat(cellfun(@(x)squeeze(x(1,:,1:2)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)'),0,'all');  %
        %             per_trial_norm{i,1}=(squeeze(conn{midx}{i,1}{2*tridx-1}(1,:,12))-centt(i,:))'./stdd(i,:);
        %             per_trial_norm{i,2}=(squeeze(conn{midx}{i,1}{2*tridx}(1,:,12))-centt(i,:))'./stdd(i,:);
        %
    end
%     su_sel=find(min(cellfun(@(x) size(x,1),per_trial_norm(:,1:2)),[],2)>NTRIAL & all(centt>0&stdd>0,2));
    su_sel=find(min(cellfun(@(x) size(x,1),per_trial_norm(:,1:2)),[],2)>NTRIAL);
    for NSU=1:length(su_sel)   
        for rpt=1:NRPT
            s1=per_trial_norm{su_sel(NSU),1}(randperm(size(per_trial_norm{su_sel(NSU),1},1),NTRIAL),:);
            s2=per_trial_norm{su_sel(NSU),2}(randperm(size(per_trial_norm{su_sel(NSU),2},1),NTRIAL),:);
            [~,tbl]=anova1([s1;s2],[zeros(NTRIAL,1);ones(NTRIAL,1)],'off');
            pev{tridx,1}(NSU,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
            
            [~,tbl]=anova1([s1;s2],randsample([zeros(NTRIAL,1);ones(NTRIAL,1)],2*NTRIAL),'off');
            pev{tridx,2}(NSU,rpt)=(tbl{2,2}-tbl{2,3}*tbl{3,4})/(tbl{4,2}+tbl{3,4});
        end 
    end
end
out=pev([1,7],:);
end

function plotLine(Data,c,l)
    m=mean(Data,1);
    plot(1:size(m,2),m,c,'LineStyle',l);
%     s=std(Data)/sqrt(size(Data,1));
%     fill([1:size(m,2),fliplr(1:size(m,2))],[m+s,fliplr(m-s)],c,'EdgeColor','none','FaceAlpha',0.2);
    ci=bootci(1000,{@(x) mean(x),Data},'type','normal');
    fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.2);
end
function [sig,pair]=load_sig_pair(Path,type)
if strcmp(type,'ext')
    fl=dir(fullfile(Path,'ext_conn_w_reg_*.mat'));
elseif strcmp(type,'inh')
    fl=dir(fullfile(Path,'inh_conn_w_reg_*.mat'));
end

sig=struct();
[sig.suid,sig.reg,sig.sess,sig.mem_type]=deal(cell(0));
pair=sig;
for fidx=1:length(fl)
    fstr=load(fullfile(fl(fidx,1).folder,fl(fidx,1).name)); 
    fields={'suid','reg','mem_type'};
    for fi=fields
        sig.(fi{1}){fidx}=fstr.sig_meta.(fi{1});
        pair.(fi{1}){fidx}=fstr.pair_meta.(fi{1});
    end
    sig.sess{fidx}=repmat(9+fidx,size(fstr.sig_meta.suid,1),1);
    pair.sess{fidx}=repmat(9+fidx,size(fstr.pair_meta.suid,1),1);
end        

for fi=[fields,{'sess'}]
    sig.(fi{1})=cell2mat(sig.(fi{1})');
    pair.(fi{1})=cell2mat(pair.(fi{1})'); 
end
end
