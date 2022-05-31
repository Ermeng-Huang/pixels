% input data from
homedir='F:\pixel-dualtask';
[sig,pair]=load_sig_pair(fullfile(homedir,'xcorr','bzdata','bzdata_zx'),'ext');
allneuron=true;

if allneuron
    [hier_stats,fh,bh]=conn_prob_bars_hier(sig,pair);
    exportgraphics(fh,fullfile(homedir,'xcorr','conn_prob_bars_inh.pdf'));
else
    % load neuron id
    load(fullfile(homedir,'auc','auc_type.mat'),'sel')
%     NeuronSet=util.ChooseNeuronSet_dualtask('region','all','sel','mixed-sel');
    [hier_stats,fh,bh]=conn_prob_bars_hier(structfun(@(x)x(all(ismember(sig.suid,sel.reverse(1).id),2),:,:),sig,'UniformOutput',false)...
        ,structfun(@(x)x(all(ismember(pair.suid,sel.reverse(1).id),2),:,:),pair,'UniformOutput',false));
    exportgraphics(fh,fullfile(homedir,'xcorr','conn_prob_bars_mixed-sel.pdf'));
end

%% function 
function [hier_stats,fh,bh]=conn_prob_bars_hier(sig,pair,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    opt.dist (1,1) double = 5
    opt.per_region_within (1,1) logical = false
    opt.hierarchy (1,:)char = 'none'
end

sess_cnt=length(unique(sig.sess));
same_stats=struct();
[same_stats.nm_nm,same_stats.congr,same_stats.incon,same_stats.mem_nm,same_stats.nm_mem]...
    =deal(nan(sess_cnt,1));
[sig_diff,sig_same,sig_direction1,sig_direction2]=diff_at_level(sig.reg,'hierarchy',opt.hierarchy);
[pair_diff,pair_same,pair_direction1,pair_direction2]=diff_at_level(pair.reg,'hierarchy',opt.hierarchy);

same_stats=get_ratio(sig.mem_type(sig_same(:,opt.dist),:),pair.mem_type(pair_same(:,opt.dist),:));
% %% %%%%%%%%%per region local FC%%%%%%%%%%
if opt.per_region_within
    same_reg_pair=squeeze(pair.reg(pair_same(:,opt.dist),opt.dist,:));
    ureg=unique(same_reg_pair(:,1));
    for ridx=1:numel(ureg)
        onereg=ureg(ridx);
        reg_sig_type=sig.mem_type(sig_same(:,opt.dist) & sig.reg(:,opt.dist,1)==onereg,:);
        reg_pair_type=pair.mem_type(pair_same(:,opt.dist) & pair.reg(:,opt.dist,1)==onereg,:);
        reg_stats(ridx)=bz.bars_util.get_ratio(reg_sig_type,reg_pair_type);
    end
    [~,~,ratiomap]=ref.get_pv_sst();
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    
    congru_corr=cell2mat(arrayfun(@(x) [ratiomap(char(idmap.ccfid2reg(ureg(x)))),reg_stats(x).congr(1)],(1:numel(ureg)).','UniformOutput',false));
    incong_corr=cell2mat(arrayfun(@(x) [ratiomap(char(idmap.ccfid2reg(ureg(x)))),reg_stats(x).incon(1)],(1:numel(ureg)).','UniformOutput',false));
    nmnm_corr=cell2mat(arrayfun(@(x) [ratiomap(char(idmap.ccfid2reg(ureg(x)))),reg_stats(x).nm_nm(1)],(1:numel(ureg)).','UniformOutput',false));
    
    figure('Color','w')
    hold on;
    ch=scatter(congru_corr(:,1),congru_corr(:,2),100,'.','r');
    ih=scatter(incong_corr(:,1),incong_corr(:,2),100,'.','b');
    nh=scatter(nmnm_corr(:,1),nmnm_corr(:,2),100,'.','k');
    xlabel('Hierarchy Index');
    ylabel('Coupling rate');
    legend([ch,ih,nh],{'Congruent','Incongruent','Nonmemory'});
end

if strcmp(opt.hierarchy,'none')
    diff_stats=get_ratio(sig.mem_type(sig_diff(:,opt.dist),:),pair.mem_type(pair_diff(:,opt.dist),:));
    hier_stats=struct('same_stats',same_stats,'diff_stats',diff_stats);
    assignin('base','hier_stats',hier_stats);
    
    fh=figure('Color','w','Position',[32,32,235,235]);
    hold on
    bh=bar([same_stats.congr(1),same_stats.incon(1),same_stats.nm_nm(1);...
        diff_stats.congr(1),diff_stats.incon(1),diff_stats.nm_nm(1)].*100);
    [ci1,ci2]=deal([]);
    for f=["congr","incon","nm_nm"]
        ci1=[ci1,cellfun(@(x) x.(f)(2),{same_stats,diff_stats})];
        ci2=[ci2,cellfun(@(x) x.(f)(3),{same_stats,diff_stats})];
    end
    errorbar([bh.XEndPoints],[bh.YEndPoints],ci1.*100-[bh.YEndPoints],ci2.*100-[bh.YEndPoints],'k.');    
    
    bh(1).FaceColor='w';
    bh(2).FaceColor=[0.5,0.5,0.5];
    bh(3).FaceColor='k';
    
    legend(bh,{'Same memory','Diff. memory','Non-memory'})
    set(gca(),'XTick',1:2,'XTickLabel',{'Within reg.','Cross region'},'XTickLabelRotation',30)
    ylabel('Func. coupling probability (%)');

else
    %h2l
    h2l_stats=get_ratio(sig.mem_type(sig_direction1(:,opt.dist),:),pair.mem_type(pair_direction1(:,opt.dist),:));
    %l2h
    l2h_stats=get_ratio(sig.mem_type(sig_direction2(:,opt.dist),:),pair.mem_type(pair_direction2(:,opt.dist),:));
    
    hier_stats=struct('same_stats',same_stats,'l2h_stats',l2h_stats,'h2l_stats',h2l_stats);
    assignin('base','hier_stats',hier_stats);
    
    fh=figure('Color','w','Position',[32,32,235,235]);
    hold on
    bh=bar([same_stats.congr(1),same_stats.incon(1),same_stats.nm_nm(1);...
        l2h_stats.congr(1),l2h_stats.incon(1),l2h_stats.nm_nm(1);...
        h2l_stats.congr(1),h2l_stats.incon(1),h2l_stats.nm_nm(1)].*100);
    [ci1,ci2]=deal([]);
    for f=["congr","incon","nm_nm"]
        ci1=[ci1,cellfun(@(x) x.(f)(2),{same_stats,l2h_stats,h2l_stats})];
        ci2=[ci2,cellfun(@(x) x.(f)(3),{same_stats,l2h_stats,h2l_stats})];
    end
    errorbar([bh.XEndPoints],[bh.YEndPoints],ci1.*100-[bh.YEndPoints],ci2.*100-[bh.YEndPoints],'k.');
    
    
    bh(1).FaceColor='w';
    bh(2).FaceColor=[0.5,0.5,0.5];
    bh(3).FaceColor='k';
    
    legend(bh,{'Same memory','Diff. memory','Non-memory'})
    set(gca(),'XTick',1:3,'XTickLabel',{'Within reg.','sensory->motor','motor->sensory'},'XTickLabelRotation',30)
    ylabel('Func. coupling probability (%)');   
end
%% chisq test
p=0;
% p=structfun(@(x)chisq_3(x.congr(4),x.congr(5),x.incon(4),x.incon(5),x.nm_nm(4),x.nm_nm(5)),hier_stats);

% chisq_3(hier_stats.same_stats.congr(4),hier_stats.same_stats.congr(5),...
%     hier_stats.l2h_stats.congr(4),hier_stats.l2h_stats.congr(5),...
% hier_stats.h2l_stats.congr(4),hier_stats.h2l_stats.congr(5))
% 
% 
% chisq_3(hier_stats.same_stats.incon(4),hier_stats.same_stats.incon(5),...
%     hier_stats.l2h_stats.incon(4),hier_stats.l2h_stats.incon(5),...
% hier_stats.h2l_stats.incon(4),hier_stats.h2l_stats.incon(5))
% 
% 
% chisq_3(hier_stats.same_stats.nm_nm(4),hier_stats.same_stats.nm_nm(5),...
%     hier_stats.l2h_stats.nm_nm(4),hier_stats.l2h_stats.nm_nm(5),...
% hier_stats.h2l_stats.nm_nm(4),hier_stats.h2l_stats.nm_nm(5))

end

function [sig,pair]=load_sig_pair(Path,type)
if strcmp(type,'ext')
    fl=dir(fullfile(Path,'0902_conn_w_reg_*.mat'));
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

function [is_diff,is_same,h2l,l2h]=diff_at_level(reg,opt)
arguments
    reg
    opt.hierarchy (1,:) char = 'none'
    
end

switch(opt.hierarchy)
    case 'pv_sst'
    [~,~,ratiomap]=ref.get_pv_sst();
    case 'OBM1'
    fstr=load('OBM1map.mat','OBM1map');
    OBM1map=fstr.OBM1map;
end
idmap=load('reg_ccfid_map.mat');

if strcmp(opt.hierarchy,'none')
    is_diff=false(size(reg,1),6);
    is_same=false(size(reg,1),6);
    h2l=[];
    l2h=[];
    graysel=all(reg(:,1,:)==567 | reg(:,1,:)==343,3);
    for dep=1:6
        is_diff(:,dep)=graysel & reg(:,dep,1)~=reg(:,dep,2);
        is_same(:,dep)=graysel & reg(:,dep,1)==reg(:,dep,2);
    end
else
    is_diff=[];
    is_same=false(size(reg,1),6);
    h2l=false(size(reg,1),6);
    l2h=false(size(reg,1),6);
    graysel=all(reg(:,2,:)==688,3);
    is_same(:,5)=graysel & reg(:,5,1)==reg(:,5,2) & reg(:,5,1)>0;
    for ri=1:size(reg,1)
        % CTX and labeled at ccf depth 7
        if any(reg(ri,5,:)==0,'all') || any(reg(ri,2,:)~=688,'all') || reg(ri,5,1)==reg(ri,5,2)
            continue
        end
        if  ~all(arrayfun(@(x) OBM1map.isKey(char(idmap.ccfid2reg(x))),squeeze(reg(ri,5,:))))
            disp(arrayfun(@(x) char(idmap.ccfid2reg(x)),squeeze(reg(ri,5,:)),'UniformOutput',false));
            continue
        end
        switch(opt.hierarchy)
            case 'pv_sst'
                dhier=diff(arrayfun(@(x) ratiomap(char(idmap.ccfid2reg(x))),squeeze(reg(ri,5,:))));
            case 'OBM1'
                dhier=diff(arrayfun(@(x) OBM1map(char(idmap.ccfid2reg(x))),squeeze(reg(ri,5,:))));
        end       
        if dhier>0
        	l2h(ri,5)=true;
        else
            h2l(ri,5)=true;
        end
    end
end
end

function out=get_ratio(sig_type,pair_type,opt)
% NM=-1,S1=1-3,S2=4-6
arguments
    sig_type (:,2) int32
    pair_type (:,2) int32
    opt.nm_mem (1,1) logical = false
end
out=struct();

nm_p=nnz(all(pair_type==-1,2));
if nm_p>0
    sig=nnz(all(sig_type==-1,2));
    [phat,pci]=binofit(sig,nm_p);
    out.nm_nm=[phat,pci,sig,nm_p];
else
    out.nm_nm=[0,0,0,0];
end

congr_p=nnz(all(ismember(pair_type,1:3),2) | all(ismember(pair_type,4:6),2));

if congr_p>0
    sig=nnz(all(ismember(sig_type,1:3),2) | all(ismember(sig_type,4:6),2));
    [phat,pci]=binofit(sig,congr_p);
    out.congr=[phat,pci,sig,congr_p];
else
    out.congr=[0,0,0,0];
end
% 
incon_p=nnz(any(ismember(pair_type(:,1),1:3),2) & any(ismember(pair_type(:,2),4:6),2));

if incon_p>0
    sig=nnz((any(ismember(sig_type(:,1),1),2) & any(ismember(sig_type(:,2),4),2))...
    |(any(ismember(sig_type(:,1),2),2) & any(ismember(sig_type(:,2),5),2)));
    sig=nnz(any(ismember(sig_type(:,1),1:3),2) & any(ismember(sig_type(:,2),4:6),2));

    [phat,pci]=binofit(sig,incon_p);
    out.incon=[phat,pci,sig,incon_p];
else
    out.incon=[0,0,0,0];
end

% if opt.nm_mem
%     mem_nm_p=nnz(pair_type(:,1)>0 & pair_type(:,2)==0);
%     if mem_nm_p>0
%         out.mem_nm=nnz(sig_type(:,1)>0 & sig_type(:,2)==0)/mem_nm_p;
%     end
%     
%     nm_mem_p=nnz(pair_type(:,1)==0 & pair_type(:,2)>0);
%     if nm_mem_p>0
%         out.nm_mem=nnz(sig_type(:,1)==0 & sig_type(:,2)>0)/nm_mem_p;
%     end
% end
end

function p=chisq_3(pos1,cnt1,pos2,cnt2,pos3,cnt3)
vec1=[ones(cnt1,1);...
    2*ones(cnt2,1);...
    3*ones(cnt3,1)];
vec2=[(1:cnt1).'>pos1;...
    (1:cnt2).'>pos2;...
    (1:cnt3).'>pos3];

[~,~,p]=crosstab(vec1,vec2);

end