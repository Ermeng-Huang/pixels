clc
clear
% input data from
opt.pair='ext';
opt.type='dualtask';
xcorr_folder='bzdata_10ms';
homedir='F:\pixel-dualtask';
%%
allneuron=false;
%%%%%%%%%%%%%%%% load FC data
[sig,pair]=load_sig_pair(fullfile(homedir,'xcorr','bzdata_10ms','bzdata_zx'),'ext');
load(fullfile(homedir,'difftype_ByPEV_1019.mat'))
NeuronSet_sel=util.ChooseNeuronSet_dualtask('sel','sel');

NeuronSet_all=util.ChooseNeuronSet_dualtask('sel','all');
NeuronSet_nonsel=util.ChooseNeuronSet_dualtask('sel','nonsel');
load(fullfile(homedir,'reactived-neuron.mat'));

id1=NeuronSet_gain.id;
id2=NeuronSet_gain.id;
% id2=NeuronSet_gain.id;
% id3=NeuronSet.id;
% NeuronSet_nonsel=util.ChooseNeuronSet_dualtask('sel','nonsel');
% id4=setdiff(NeuronSet_nonsel.id,id3);
% sig1=structfun(@(x)x(ismember(sig.suid(:,1),id1)&ismember(sig.suid(:,2),id3),:,:),sig,'UniformOutput',false);
% pair1=structfun(@(x)x(ismember(pair.suid(:,1),id1)&ismember(pair.suid(:,2),id3),:,:),pair,'UniformOutput',false);
sig0=structfun(@(x)x(ismember(sig.suid(:,1),id1)&ismember(sig.suid(:,2),id2)...
    |ismember(sig.suid(:,1),id2)&ismember(sig.suid(:,2),id1),:,:),sig,'UniformOutput',false);
pair0=structfun(@(x)x(ismember(pair.suid(:,1),id1)&ismember(pair.suid(:,2),id2)...
    |ismember(pair.suid(:,1),id2)&ismember(pair.suid(:,2),id1),:,:),pair,'UniformOutput',false);

%%%%%%%%%%%%%%%% region and region order
idmap=load('reg_ccfid_map.mat');
grey_regs=util.getGreyRegs('range','CTX','mincount',100);
regs=grey_regs;
regs={{'PIR','TT','AON','DP'};{'EPd','EPv'};{'ACA','AI','ILA','ORB','PL'};{'HIP','RHP'};{'MO'};{'SS'}};
reg_label={'Olfactory','Cortical subplate','Association cortex','Hippocampus','Motor cortex','Somatosensory cortex'};
% 
% arrayfun(@(x))
% f=load(fullfile(homedir,'allen_index_zxx.mat'));
% reg_index(:,1)=f.intersect_regs;
% reg_index(:,2)=num2cell(f.index(2,:));
% idmap=load('reg_ccfid_map.mat');
% region_level='CTX';
% switch(region_level)
%     case 'CH'
%         intersect_regs=f.intersect_regs(cellfun(@(y)strcmp(y(:,3),'CH'),cellfun(@(x)idmap.reg2tree(x),f.intersect_regs,'UniformOutput',false)),:);
%         reg_index=reg_index(cellfun(@(y)strcmp(y(:,3),'CH'),cellfun(@(x)idmap.reg2tree(x),f.intersect_regs,'UniformOutput',false)),:);
%     case 'CTX'
%         intersect_regs=f.intersect_regs(cellfun(@(y)strcmp(y(:,4),'CTX'),cellfun(@(x)idmap.reg2tree(x),f.intersect_regs,'UniformOutput',false)),:);
%         reg_index=reg_index(cellfun(@(y)strcmp(y(:,4),'CTX'),cellfun(@(x)idmap.reg2tree(x),f.intersect_regs,'UniformOutput',false)),:);
% end
% 
% %%%
% [~,n]=sort(cell2mat(reg_index(:,2)),'ascend');
% regs=reg_index(n,1);
% regs(ismember(regs,{'AUD','VIS','RSP'}))=[];
% regs=intersect(util.getGreyRegs('range','CH','mincount',100),unique(gain_all.reg));

%%%%%%%%%%%%%%%%

%%
for r1=1:length(regs)
    for r2=1:length(regs)
        %     fra(r1,r2)=nnz(ismember(sig1.reg(:,5,1),idmap.reg2ccfid(regs{r1}))&ismember(sig1.reg(:,5,2),idmap.reg2ccfid(regs{r2})))...
        %         /nnz(ismember(pair1.reg(:,5,1),idmap.reg2ccfid(regs{r1}))&ismember(pair1.reg(:,5,2),idmap.reg2ccfid(regs{r2})));
        if nnz(ismember(pair0.reg(:,5,1),cellfun(@(x)idmap.reg2ccfid(x),regs{r1}))...
                &ismember(pair0.reg(:,5,2),cellfun(@(x)idmap.reg2ccfid(x),regs{r2})))<10
            fra0(r1,r2)=nan;
        else
            fra0(r1,r2)=nnz(ismember(sig0.reg(:,5,1),cellfun(@(x)idmap.reg2ccfid(x),regs{r1}))&ismember(sig0.reg(:,5,2),cellfun(@(x)idmap.reg2ccfid(x),regs{r2})))...
                /nnz(ismember(pair0.reg(:,5,1),cellfun(@(x)idmap.reg2ccfid(x),regs{r1}))&ismember(pair0.reg(:,5,2),cellfun(@(x)idmap.reg2ccfid(x),regs{r2})));
        end
        
    end
end
%%
fh=figure('Color','w','Position',[100,100,256,256]);
% ridx=2; %%1-AON,2-PIR
[m,s]=sort(diag(fra0),'descend');
h=imagesc(fra0(s,s)*100);
set(h,'alphadata',~isnan(fra0));
set(gca,'XDir','normal','YDir','normal','XTick',1:length(regs),'XTickLabel',reg_label(s),'YTick',1:length(regs),'XTickLabelRotation',30,'YTickLabel',reg_label(s),'FontSize',6)

colorbar
% caxis([0,0.1])
set(gca,'Position',[0.3750 0.2465 0.48 0.48])
exportgraphics(fh,fullfile(homedir,'xcorr','conn_prob_reg','reac_ext_10ms.pdf'));

%%
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
