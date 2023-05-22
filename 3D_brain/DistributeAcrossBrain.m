% brain
close all

homedir='D:\pixel-dualtask\';
load('D:\code\BrainMesh-master\BrainMesh-master\977.mat')
v = v/10; % change to pixel/slice value
fh=figure('position',[729,387,320,240]) ;
hold on
view(7,-50)
% view(355,-25)
% patch('Vertices',v','Faces',F','EdgeColor','none','FaceColor',[235,164,100]/256,'FaceAlpha',0.1);
patch('Vertices',v','Faces',F','EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);

%%
% load('D:\histology\sucoords.mat')
% load(fullfile(homedir,'difftype_ByPEV_1019.mat'))
% cluster_id=h5read(fullfile(homedir,'Selectivity_0606.hdf5'),'/cluster_id');
% coord1=coord(ismember(cluster_id,lost_all.id),:);
% plot3(coord1(:,1),coord1(:,2),coord1(:,3),'b.')
% coord1=coord(ismember(cluster_id,gain_all.id),:);
% plot3(coord1(:,1),coord1(:,2),coord1(:,3),'r.')
% set(gca, 'ZDir', 'reverse')
% 
% axis(gca, 'equal');
% axis(gca, 'vis3d');
% axis(gca, 'off');
get(gca, 'Parent');
%%
id_map=load('reg_ccfid_map.mat');
% map_cells=map_sel(homedir,'p_criteria',0.05);% featii1-5:sample, distarctor ,mixed/all,mixed/sel,mixed/dist
map_cells=map_sel_lost_gain(homedir);% featii1-3:lost, gain, maintain 
% map_cells=map_sel_reac(homedir);
grey_regs=util.getGreyRegs('range','CH','mincount',100);
intersect_regs=intersect(grey_regs,keys(map_cells{1}));
feat_prop=cell2mat(map_cells{1}.values(intersect_regs))*100;
% color bar的范围，这里可能需要手动调整
Min=0;
Max=65;
map=colormap('jet'); %颜色范围
%
%%
for i=1:length(intersect_regs)
    [v,F]=util.loadawobj(['D:/code/BrainMesh-master/BrainMesh-master/Data/Allen_obj_files/',num2str(id_map.reg2ccfid(intersect_regs{i})),'.obj']);   
    patch('Vertices',v'/10,'Faces',F','EdgeColor','none','FaceColor',map(round(256*(feat_prop(i)-Min)/(Max-Min)),:),'FaceAlpha',0.3);
    
end
% colorbar
caxis([Min,Max])

colormap('jet')
axis('equal'); axis('off'); axis('tight'); axis('vis3d');
exportgraphics(fh,fullfile(homedir,'Distribution3D','Distribution3D_lost.pdf'),'ContentType','image')% 保存成vector内存过大


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
    mixed1(i)=nnz(strcmp(reg,r(i))&any(samp_sel(:,1:2),2)&any(dist_sel(:,1:2),2))/nnz(strcmp(reg,r(i)));
    mixed2(i)=nnz(strcmp(reg,r(i))&any(samp_sel(:,1:2),2)&any(dist_sel(:,1:2),2))/nnz(strcmp(reg,r(i))&any(samp_sel(:,1:2),2));
    mixed3(i)=nnz(strcmp(reg,r(i))&any(samp_sel(:,1:2),2)&any(dist_sel(:,1:2),2))/nnz(strcmp(reg,r(i))&any(dist_sel(:,1:2),2));
end
map_cells{1,1}=containers.Map(r,samp); %{ccfid:proj_dense_mat_idx}
map_cells{1,2}=containers.Map(r,dist); %{ccfid:proj_dense_mat_idx}
map_cells{1,3}=containers.Map(r,mixed1); %{ccfid:proj_dense_mat_idx}
map_cells{1,4}=containers.Map(r,mixed2); %{ccfid:proj_dense_mat_idx}
map_cells{1,5}=containers.Map(r,mixed3); %{ccfid:proj_dense_mat_idx}

end

function map_cells=map_sel_lost_gain(homedir)

% d2=load(fullfile(homedir,'pev_perform_d2.mat'));
% d3=load(fullfile(homedir,'pev_perform.mat'));

load(fullfile(homedir,'difftype_ByPEV_1019.mat'))
% load(fullfile(homedir,'pev_diff_delay.mat'))
NeuronSet_d2=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','d2');
NeuronSet_d3=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','d3');
NeuronSet_sel=util.ChooseNeuronSet_dualtask('region','all','regionLevel',7,'sel','sel');
NeuronSet_d2d3=structfun(@(x)x(ismember(NeuronSet_sel.id,union(NeuronSet_d2.id,NeuronSet_d3.id))),NeuronSet_sel,'UniformOutput',false);
reg=unique(NeuronSet_d2d3.reg);
for r=1:length(reg)
%     fra(r,1)=nnz(ismember(lost.md.reg,reg_index{r,1}))/nnz(ismember(NeuronSet_all.reg,reg_index{r,1}));
%     fra1(r,1)=nnz(ismember(lost.md.reg,reg_index{r,1}))/nnz(ismember(NeuronSet_d2.reg,reg_index{r,1}));
    fra1(r,1)=nnz(ismember(lost_all.reg,reg{r,1}))/nnz(ismember(NeuronSet_d2d3.reg,reg{r,1}));
end
for r=1:length(reg)
%     fra2(r,1)=nnz(ismember(gain.md.reg,reg_index{r,1}))/nnz(ismember(NeuronSet_d2.reg,reg_index{r,1}));
    fra2(r,1)=nnz(ismember(gain_all.reg,reg{r,1}))/nnz(ismember(NeuronSet_d2d3.reg,reg{r,1}));
end
for r=1:length(reg)
%     fra2(r,1)=nnz(ismember(gain.md.reg,reg_index{r,1}))/nnz(ismember(NeuronSet_d2.reg,reg_index{r,1}));
    fra3(r,1)=nnz(ismember(invarient_all.reg,reg{r,1}))/nnz(ismember(NeuronSet_d2d3.reg,reg{r,1}));
end


map_cells{1}=containers.Map(reg,fra1); %{ccfid:proj_dense_mat_idx}
map_cells{2}=containers.Map(reg,fra2); %{ccfid:proj_dense_mat_idx}
map_cells{3}=containers.Map(reg,fra3); %{ccfid:proj_dense_mat_idx}



end

function map_cells=map_sel_reac(homedir)
NeuronSet_all=util.ChooseNeuronSet_dualtask('sel','all');
load(fullfile(homedir,'reactived-neuron.mat'));
reg=unique(NeuronSet_gain.reg);

for r=1:length(reg)
    frac(r)=nnz(ismember(NeuronSet_all.id,NeuronSet_gain.id)&ismember(NeuronSet_all.reg,reg{r}))...
        /nnz(ismember(NeuronSet_all.reg,reg{r}));
end
map_cells{1}=containers.Map(reg,frac); %{ccfid:proj_dense_mat_idx}

end