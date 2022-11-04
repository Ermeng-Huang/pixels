clear
fh=connectivity_proportion_GLM;

%% function
function fh=connectivity_proportion_GLM(opt)
arguments    
    opt.skip_efferent (1,1) logical = true;
    opt.corr1 (1,1) logical =true;
    opt.plot1 (1,1) logical =false;
    opt.corr2 (1,1) logical =false;
    opt.plot2 (1,1) logical =false;
    opt.corr3 (1,1) logical =false;
    opt.plot3 (1,1) logical =false;
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'CH'
    opt.homedir(1,:) char ='F:\pixel-dualtask';
    opt.plot_all(1,1) logical =true;
end

%map_cells from K:\code\jpsth\+ephys\Both_either_reg_bars.m
% For duration, from: K:\code\jpsth\+ephys\duration_reg_bars.m
corr_type='PearsonLinearLog';
idmap=load('reg_ccfid_map.mat');
grey_regs=getGreyRegs('range','CTX','mincount',100);

sink_ccfid=h5read(fullfile('D:\code-hem','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('D:\code-hem','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('D:\code-hem','proj_mat.hdf5'),'/src_target_matrix');
% grey_regs assume in work space
src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

allen_src_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_src_regs,grey_regs);

%% -> feature_region_map entry point
one_reg_corr_list=[];
two_reg_corr_list=[];
featii=1; % sample, distarctor
epochii=1;
%%%% load data
phase='gain';
% map_cells=map_sel(opt.homedir,'p_criteria',0.05);
map_cells=map_sel_mod(opt.homedir,phase);
% map_cells=map_sel_lost_gain(opt.homedir);
%%%%

feat_reg_map=map_cells{featii}; %
intersect_regs=intersect(intersect_regs,keys(feat_reg_map));

feat_prop_cell=feat_reg_map.values(intersect_regs);
feat_prop=cellfun(@(x) x(1),feat_prop_cell);

allen_density=sink_src_mat(cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)))...
    ,cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs))));
[coeff,score,~,~,explained,~]=pca(allen_density,'NumComponents',3);
%% one region corr, bar plot
fh=figure('Color','w','Position',[100,100,600,225]);

for i=1:3    
    subplot(1,3,i)
    hold on
    [r,p]=corr(coeff(:,i),feat_prop,'type','Pearson');
    for ll=1:numel(feat_prop)
        c=util.getRegColor(intersect_regs{ll},'large_area',true);
        scatter(coeff(ll,i),feat_prop(ll)*100,4,c,'filled','o')
        text(coeff(ll,i),feat_prop(ll)*100,intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
    end
    xlabel(sprintf('PC %d, explained var.=%.1f%',i,explained(i)));
    ylabel('Proportion (%)')
    %                 set(gca,'XScale','log')
    title(sprintf(' r = %.3f, p = %.3f',r,p));
end

    exportgraphics(fh,fullfile(opt.homedir,'connectivity_proportion','frac_allen_scatter_PCA_%s.pdf'))
end




%% function
function map_cells=map_sel(homedir,opt)
arguments
    homedir (1,:) char
    opt.p_criteria (1,1) double = 0.01
    opt.reg_correct (1,1) logical = true
end
switch(opt.p_criteria)
    case 0.01                
        samp_sel=h5read(fullfile(homedir,'Selectivity_0602.hdf5'),'/sust_trans_noPermutaion');
        dist_sel=h5read(fullfile(homedir,'distractor_0605.hdf5'),'/sust_trans_noPermutaion');
    case 0.05
        load(fullfile(homedir,'cells_task.mat'));
        samp_sel=h5read(fullfile(homedir,'Selectivity_0606.hdf5'),'/sust_trans_noPermutaion');
        dist_sel=h5read(fullfile(homedir,'distractor_0606.hdf5'),'/sust_trans_noPermutaion'); 
        samp_sel(~(cells_task),1:2)=0;
        dist_sel(~(cells_task),1:2)=0;
    case 0
        dist_sel1=h5read(fullfile(homedir,'Selectivity_1119.hdf5'),'/sust_trans_noPermutaion'); 
end
if opt.reg_correct
    load(fullfile(homedir,'regs_w2g.mat'),'reg_new');
    reg=reg_new;
else
    reg=h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/reg');
end

r=unique(reg(~strcmp(reg(:,7),''),7));
for i=1:numel(r)
    samp(i)=nnz(strcmp(reg,r(i))&any(samp_sel(:,1:2),2))/nnz(strcmp(reg,r(i)));
    dist(i)=nnz(strcmp(reg,r(i))&any(dist_sel(:,1:2),2))/nnz(strcmp(reg,r(i)));
end
map_cells{1,1}=containers.Map(r,samp); %{ccfid:proj_dense_mat_idx}
map_cells{1,2}=containers.Map(r,dist); %{ccfid:proj_dense_mat_idx}

end

function map_cells=map_sel_mod(homedir,type)
load(fullfile(homedir,'difftype_ByPEV.mat'))
NeuronSet_sel=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');
NeuronSet_d2=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','d2');
NeuronSet_d3=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','d3');
NeuronSet_d2d3=structfun(@(x)x(ismember(NeuronSet_sel.id,union(NeuronSet_d2.id,NeuronSet_d3.id))),NeuronSet_sel,'UniformOutput',false);

switch(type)
    case 'lost-md'
        reg=unique(lost.md.reg);
        reg(cellfun(@(x)nnz(ismember(NeuronSet_d2.reg,x))<5,reg))=[];
        for r=1:length(reg)
            fra(r,1)=nnz(ismember(lost.md.reg,reg{r,1}))/nnz(ismember(NeuronSet_d2.reg,reg{r,1}));
        end
    case 'gain-md'
        reg=unique(gain.md.reg);
        reg(cellfun(@(x)nnz(ismember(NeuronSet_d2.reg,x))<5,reg))=[];
        for r=1:length(reg)
            fra(r,1)=nnz(ismember(gain.md.reg,reg{r,1}))/nnz(ismember(NeuronSet_d2.reg,reg{r,1}));
        end
    case 'lost-ld'
        reg=unique(lost.ld.reg);
        reg(cellfun(@(x)nnz(ismember(NeuronSet_d2.reg,x))<5,reg))=[];
        for r=1:length(reg)
            fra(r,1)=nnz(ismember(lost.ld.reg,reg{r,1}))/nnz(ismember(NeuronSet_d3.reg,reg{r,1}));
        end
    case 'gain-ld'
        reg=unique(gain.ld.reg);
        reg(cellfun(@(x)nnz(ismember(NeuronSet_d3.reg,x))<5,reg))=[];
        for r=1:length(reg)
            fra(r,1)=nnz(ismember(gain.ld.reg,reg{r,1}))/nnz(ismember(NeuronSet_d3.reg,reg{r,1}));
        end
    case 'stable-md'
        reg=union(unique(lost.md.reg),unique(gain.md.reg));
        reg(cellfun(@(x)nnz(ismember(NeuronSet_d2.reg,x))<5,reg))=[];
        for r=1:length(reg)
            fra(r,1)=(nnz(ismember(gain.md.reg,reg{r,1}))+nnz(ismember(lost.md.reg,reg{r,1})))/nnz(ismember(NeuronSet_d2.reg,reg{r,1}));
        end
        fra=1-fra;
    case 'stable-ld'
        reg=union(unique(lost.ld.reg),unique(gain.ld.reg));
        reg(cellfun(@(x)nnz(ismember(NeuronSet_d3.reg,x))<5,reg))=[];
        for r=1:length(reg)
            fra(r,1)=(nnz(ismember(gain.ld.reg,reg{r,1}))+nnz(ismember(lost.ld.reg,reg{r,1})))/nnz(ismember(NeuronSet_d3.reg,reg{r,1}));
        end
        fra=1-fra;
    case 'gain'
        reg=unique(gain_all.reg);        
        for r=1:length(reg)
            fra(r,1)=nnz(ismember(gain_all.reg,reg{r,1}))/nnz(ismember(NeuronSet_d2d3.reg,reg{r,1}));
        end
        
end
map_cells{1}=containers.Map(reg,fra); %{ccfid:proj_dense_mat_idx}
end

function map_cells=map_sel_reac(homedir)
load('F:\pixel-dualtask\pev-reactive-correct-error.mat')
for i=1:2912
    c=pev{1}(i,~isnan(pev{1}(i,:)));
    e=pev{2}(i,~isnan(pev{2}(i,:)));
    if length(c)<250||length(e)<250
        p(i)=100;
    else
        %         p(i)=ranksum(c,e,'tail','right');
        %         p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
        if mean(c)>mean(e)
            p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
        else
            p(i)=100;
        end
    end    
end
NeuronSet1=structfun(@(x)x(p<0.05),NeuronSet,'UniformOutput',false);

reg=unique(NeuronSet1.reg);
reg(cellfun(@(x)nnz(ismember(NeuronSet.reg,x))<5,reg))=[];
for r=1:length(reg)
    fra(r,1)=nnz(ismember(NeuronSet1.reg,reg{r,1}))/nnz(ismember(NeuronSet.reg,reg{r,1}));
end

map_cells{1}=containers.Map(reg,fra); %{ccfid:proj_dense_mat_idx}
end

function map_cells=map_sel_lost_gain(homedir)

homedir='F:/pixel-dualtask';
d2=load(fullfile(homedir,'pev_perform_d2.mat'));
d3=load(fullfile(homedir,'pev_perform.mat'));

load(fullfile(homedir,'difftype_ByPEV.mat'),'lost','gain')
% load(fullfile(homedir,'pev_diff_delay.mat'))
NeuronSet_sel=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');

lost_all=structfun(@(x)x(ismember(NeuronSet_sel.id,setdiff(union(lost.md.id,lost.ld.id),union(gain.md.id,gain.ld.id)))),NeuronSet_sel,'UniformOutput',false);
gain_all=structfun(@(x)x(ismember(NeuronSet_sel.id,setdiff(union(gain.md.id,gain.ld.id),union(lost.md.id,lost.ld.id)))),NeuronSet_sel,'UniformOutput',false);

for i=1:length(d2.NeuronSet.id)
    c=squeeze(d2.pev{3}(i,1,~isnan(d2.pev{3}(i,1,:))));
    e=squeeze(d2.pev{4}(i,1,~isnan(d2.pev{4}(i,1,:))));
    if length(c)<250||length(e)<250
        p1(i)=100;
    else
%         p(i)=ranksum(c,e,'tail','right');
%         p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
        if mean(c)>mean(e)
            p1(i)=anova1([c;e],[zeros(length(c),1);ones(length(e),1)],'off');
        else
            p1(i)=100;
        end
    end    
end
for i=1:length(d3.NeuronSet.id)
    c=squeeze(d2.pev{3}(i,1,~isnan(d2.pev{3}(i,1,:))));
    e=squeeze(d2.pev{4}(i,1,~isnan(d2.pev{4}(i,1,:))));
    if length(c)<250||length(e)<250
        p2(i)=100;
    else
%         p(i)=ranksum(c,e,'tail','right');
%         p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
        if mean(c)>mean(e)
            p2(i)=anova1([c;e],[zeros(length(c),1);ones(length(e),1)],'off');
        else
            p2(i)=100;
        end
    end    
end
d2.NeuronSet1=structfun(@(x)x(p1<0.05),d2.NeuronSet,'UniformOutput',false);
d3.NeuronSet1=structfun(@(x)x(p2<0.05),d3.NeuronSet,'UniformOutput',false); 

NeuronSet_d2d3=structfun(@(x)x(ismember(NeuronSet_sel.id,union(d2.NeuronSet1.id,d3.NeuronSet1.id))),NeuronSet_sel,'UniformOutput',false);
lost_perf=structfun(@(x)x(ismember(NeuronSet_d2d3.id,setdiff(union(lost.md.id,lost.ld.id),union(gain.md.id,gain.ld.id)))),NeuronSet_d2d3,'UniformOutput',false);
gain_perf=structfun(@(x)x(ismember(NeuronSet_d2d3.id,setdiff(union(gain.md.id,gain.ld.id),union(lost.md.id,lost.ld.id)))),NeuronSet_d2d3,'UniformOutput',false);

reg=unique(lost_perf.reg);
reg(cellfun(@(x)nnz(ismember(lost_all.reg,x))<10,reg))=[];
for r=1:length(reg)
    fra1(r,1)=nnz(ismember(lost_perf.reg,reg{r,1}))/nnz(ismember(lost_all.reg,reg{r,1}));
end
map_cells{1}=containers.Map(reg,fra1); %{ccfid:proj_dense_mat_idx}
reg=unique(gain_perf.reg);
reg(cellfun(@(x)nnz(ismember(gain_all.reg,x))<5,reg))=[];
for r=1:length(reg)
    fra2(r,1)=nnz(ismember(gain_perf.reg,reg{r,1}))/nnz(ismember(gain_all.reg,reg{r,1}));
end
map_cells{2}=containers.Map(reg,fra2); %{ccfid:proj_dense_mat_idx}




end


function grey_regs_=getGreyRegs(opt)
arguments
    opt.mincount (1,1) double = 100
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.reg_correct (1,1) logical = true
    opt.homedir(1,:) char ='F:\pixel-dualtask';
end
persistent grey_regs opt_
if isempty(grey_regs) || ~isequaln(opt_,opt)    
    if opt.reg_correct
        load(fullfile(opt.homedir,'regs_w2g.mat'),'reg_new');
        reg=reg_new;
    else
        reg=h5read(fullfile(opt.homedir,'Selectivity_1129.hdf5'),'/reg');
    end
    meta.reg_tree=reg(:,3:7)';
    
    BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
    CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
    CTXsel=strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),'');
    switch opt.range
        case 'grey'
            grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));
        case 'CH'
            grey_regs=unique(meta.reg_tree(5,CHsel));
        case 'CTX'
            grey_regs=unique(meta.reg_tree(5,CTXsel));
    end
    cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
    grey_regs=grey_regs(cnt>opt.mincount);
end
grey_regs_=grey_regs;
opt_=opt;
end