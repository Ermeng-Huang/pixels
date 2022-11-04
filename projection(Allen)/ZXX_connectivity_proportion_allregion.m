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
only=true;
corr_type='PearsonLinearLog';
idmap=load('reg_ccfid_map.mat');
grey_regs=getGreyRegs('range','CH','mincount',100);

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
featii=1; 
epochii=1;
%%%% load data
% phase='gain';
% map_cells=map_sel(opt.homedir,'p_criteria',0.05);% featii1-5:sample, distarctor ,mixed/all,mixed/sel,mixed/dist
% map_cells=map_sel_mod(opt.homedir,phase);
% map_cells=map_sel_lost_gain(opt.homedir); % featii1-2:lost, gain ,mixed/all,mixed/sel,mixed/dist
map_cells=map_sel_reac(opt.homedir);
%%%%

feat_reg_map=map_cells{featii}; %
intersect_regs=intersect(intersect_regs,keys(feat_reg_map));
feat_prop_cell=feat_reg_map.values(intersect_regs);
feat_prop=cellfun(@(x) x(1),feat_prop_cell);

%%% 提取allen里面各脑区的投射密度
% if only %true：只plot实验记录到脑区 false：allen里有的脑区
%     sink_ccfid=intersect(sink_ccfid,cellfun(@(x)idmap.reg2ccfid(x),intersect_regs));
% end
idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs))); %%% 脑区对应的在sink_src_mat的列号

[glmxmat,glmxmeta]=deal([]);
% efferent_proj_dense(soruce) ~ feature_proportion
for ii=1:numel(sink_ccfid)
    %     allen_density=log10(sink_src_mat(ii,idx4corr)+1e-12); % from idx4corr, to one alternating target
    allen_density=sink_src_mat(ii,idx4corr);
    if nnz(allen_density~=0)<10
        continue
    end
    glmxmat=[glmxmat;allen_density];
    glmxmeta=[glmxmeta;1,ii,sink_ccfid(ii)];
end

% afferent_proj_dense(sink) ~ feature_proportion

idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
for ii=1:numel(src_ccfid)
    allen_density=sink_src_mat(idx4corr,ii); % from one alternating target, to idx4corr
    if nnz(allen_density~=0)<10
        disp(src_ccfid)
        continue
    end
    glmxmat=[glmxmat;allen_density.'];
    glmxmeta=[glmxmeta;2,ii,src_ccfid(ii)];
end

%% one region corr, bar plot
if opt.corr1
%     glmxmat=log(glmxmat);
    for regii=1:size(glmxmat,1)
        switch(corr_type)
            case 'Spearman'
                [r,p]=corr(glmxmat(regii,:).',feat_prop,'type','Spearman');
            case 'PearsonLogLog'
                vsel=(glmxmat(regii,:).')>0 & feat_prop>0;
                [r,p]=corr(reshape(log10(glmxmat(regii,vsel)),[],1),...
                    reshape(log10(feat_prop(vsel)),[],1),'type','Pearson');
            case 'PearsonLinearLog'
                vsel=(glmxmat(regii,:).')>0;
                [r,p]=corr(reshape(log10(glmxmat(regii,vsel)),[],1),...
                    reshape(feat_prop(vsel),[],1),'type','Pearson');
        end
        one_reg_corr_list=[one_reg_corr_list;featii,epochii,regii,double(glmxmeta(regii,:)),r,p];
        %====================================^^^1^^^^^^2^^^^^^3^^^^^^^^^^^4,5,6^^^^^^^^^^^^^7^8^^
    end
    if only
        s_list=sortrows(one_reg_corr_list(one_reg_corr_list(:,1)==featii...
            &ismember(one_reg_corr_list(:,6),cell2mat(idmap.reg2ccfid.values(intersect_regs))),:),7,'descend');
    else
        s_list=sortrows(one_reg_corr_list(one_reg_corr_list(:,1)==featii,:),7,'descend');
    end

    if opt.plot_all
        %         %                 s_list=s_list(isfinite(s_list(:,7)),:);
        fh=figure('Color','w','Position',[32,32,135,225]);
        for i=2
%             subplot(1,2,i)            
            s_list1=s_list(s_list(:,4)==i,:);
            s_list1=s_list1(1:10,:);
            xlbl=idmap.ccfid2reg.values(num2cell(s_list1(:,6)));
            xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);
            
            hold on
            bh=bar(s_list1(:,7),0.8,'FaceColor','flat','Horizontal','on');
            bh.CData(s_list1(:,4)==1,:)=repmat([0.5,0.5,0.5],nnz(s_list1(:,4)==1),1);
            bh.CData(s_list1(:,4)==2,:)=repmat([0,0,0],nnz(s_list1(:,4)==2),1);
            xlabel('Connectivity-coding proportion correlation (Pearson''s r)')
            set(gca(),'YTick',1:size(s_list1,1),'YTickLabel',xlbl,'YDir','reverse');
            xlim([0,1]);
%             legend([bhto,bhfrom],{'Connectivity to this region','Connectivity from this region'})
            if i==1
                title('Connectivity to this region') 
            else
               title('Connectivity from this region') 
            end
            %             exportgraphics(fh.(sprintf('bar1_feat%d',featii)),sprintf('frac_allen_mdls_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
        end
        exportgraphics(fh,fullfile(opt.homedir,'connectivity_proportion',sprintf('frac_allen_mdls_%s.pdf','reac')))
        for i=2            
            s_list1=s_list(s_list(:,4)==i,:);
            s_list1=s_list1([1:5,end-4:end],:);
            xlbl=idmap.ccfid2reg.values(num2cell(s_list1(:,6)));
            xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);
            
            fh=figure('Color','w','Position',[32,32,1020,420]);            
            for ii=1:10
                subplot(2,5,ii)
                hold on
                glmidx=find(glmxmeta(:,1)==s_list1(ii,4) & glmxmeta(:,2)==s_list1(ii,5)); %%%%%%%%%%%%%%%%
                for ll=1:numel(feat_prop)
                    c=util.getRegColor(intersect_regs{ll},'large_area',true);
                    scatter(log10(glmxmat(glmidx,ll)),feat_prop(ll)*100,4,c,'filled','o')
                    text(log10(glmxmat(glmidx,ll)),feat_prop(ll)*100,intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
                end
                xlabel(sprintf('%s projection density',xlbl{ii}));
                ylabel('Proportion (%)')
%                 set(gca,'XScale','log')
                title(sprintf(' r = %.3f, p = %.3f',s_list1(ii,7),s_list1(ii,8)));
                % xlim([0.15,0.5])
                % ylim([0.15,0.5])
                %             exportgraphics(fh.(sprintf('corr1_feat%d',featii)),sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
                
                %             close(fh)
            end
            if i==1
                sgtitle('Connectivity to this region')
%                 exportgraphics(fh,fullfile(opt.homedir,'connectivity_proportion',sprintf('frac_allen_scatter_ToRegion_%s.png',phase)))
            else
                sgtitle('Connectivity from this region')
                exportgraphics(fh,fullfile(opt.homedir,'connectivity_proportion',sprintf('frac_allen_scatter_FromRegion_%s.pdf','reac')))
            end
         end
    end
end
if opt.corr2
    comb2=nchoosek(1:size(glmxmat,1),2);
    mdlid_rsq_AICC=[];
    for jj=1:size(comb2,1)
        if rem(jj,1000)==0
            disp(jj)
        end
        allen_A=glmxmat(comb2(jj,1),:).'; % from one alternating target, to idx4corr
        allen_B=glmxmat(comb2(jj,2),:).';
        mdl=fitglm([allen_A,allen_B],feat_prop,'interactions'); % (1|2)jj => glmxmat==glmxmeta => (sink_ccfid|src_ccfid)
        two_reg_corr_list=[two_reg_corr_list;featii,epochii,jj,double(glmxmeta(comb2(jj,1),:)),double(glmxmeta(comb2(jj,2),:)),mdl.Coefficients.Estimate.',mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc];
        %============================================^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        mdlid_rsq_AICC=[mdlid_rsq_AICC;jj,mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc];
        %     if mdl.Rsquared.Ordinary>maxrsq
        %         maxrsq=mdl.Rsquared.Ordinary;
        %         maxidx=jj;
        %     end
    end
    if opt.plot2
        mdlid_rsq_AICC=sortrows(mdlid_rsq_AICC,3);
        % keyboard();
        ytk=cell(0,0);
        for kk=1:10
            maxidx=mdlid_rsq_AICC(kk,1);
            if glmxmeta(comb2(maxidx,1),1)==1
                reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,1),2)));
                reg1=['To ',reg1{1}];
            else
                reg1=(idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,1),2))));
                reg1=['From ',reg1{1}];
            end
            if glmxmeta(comb2(maxidx,2),1)==1
                reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,2),2)));
                reg2=['To ',reg2{1}];
            else
                reg2=idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,2),2)));
                reg2=['From ',reg2{1}];
            end
            disp([reg1,'-',reg2,' ',num2str(mdlid_rsq_AICC(kk,2:3))])
            ytk{end+1}=[reg1,'-',reg2];
        end
        
        fh=figure('Color','w','Position',[32,32,400,800]);
        bar(sqrt(mdlid_rsq_AICC(1:10,2)),'Horizontal','on')
        set(gca(),'YDir','reverse','YTick',1:10,'YTickLabel',ytk(1:10))
        xlabel('Person''s r')
        xlim([0,1])
        title(sprintf('selectivity %d - epoch %d',featii,epochii));
%         exportgraphics(fh,sprintf('frac_allen_mdls_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
%         close(fh)
        %%
        maxidx=mdlid_rsq_AICC(1,1);
        allen_A=glmxmat(comb2(maxidx,1),:).'; % from one alternating target, to idx4corr
        allen_B=glmxmat(comb2(maxidx,2),:).';
        mdl=fitglm([allen_A,allen_B],feat_prop,'interactions');
        
        
        disp(sqrt(mdl.Rsquared.Ordinary));
        if glmxmeta(comb2(maxidx,1),1)==1
            reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,1),2)));
            reg1=['To_',reg1{1}];
        else
            reg1=(idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,1),2))));
            reg1=['From_',reg1{1}];
        end
        if glmxmeta(comb2(maxidx,2),1)==1
            reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,2),2)));
            reg2=['To_',reg2{1}];
        else
            reg2=idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,2),2)));
            reg2=['From_',reg2{1}];
        end
        
        fh=figure('Color','w','Position',[32,32,320,320]);
        hold on
        for ll=1:numel(feat_prop)
            c=util.getRegColor(intersect_regs{ll},'large_area',true);
            scatter(mdl.Fitted.Response(ll),feat_prop(ll),4,c,'filled','o')
            text(mdl.Fitted.Response(ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
        end
        title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s] +%.3f [%s]', ...
            mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2,mdl.Coefficients.Estimate(4),'interaction'), ...
            'Interpreter','none')
        xlabel(sprintf('Model %d prediction',maxidx))
        ylabel('Proportion of coding neuron')
        set(gca,'XScale','log','YScale','log')
        text(min(xlim()),max(ylim()),sprintf(' r = %.3f, AIC = %.1f',sqrt(mdl.Rsquared.Ordinary),mdl.ModelCriterion.AIC),'HorizontalAlignment','left','VerticalAlignment','top');
        % xlim([0.15,0.5])
        % ylim([0.15,0.5])
%         exportgraphics(fh,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
%         close(fh)
        %         keyboard()
        
    end
    save('two_reg_corr_list.mat','two_reg_corr_list')
end
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
load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet');

reg=unique(NeuronSet.reg);
reg(cellfun(@(x)nnz(ismember(NeuronSet.reg,x))<5,reg))=[];
NeuronSet1=util.ChooseNeuronSet_dualtask('sel','all');
for r=1:length(reg)
    fra(r,1)=nnz(ismember(NeuronSet.reg,reg{r,1}))/nnz(ismember(NeuronSet1.reg,reg{r,1}));
end
map_cells{1}=containers.Map(reg,fra); %{ccfid:proj_dense_mat_idx}

% for r=1:length(reg_index)
%     fra(r,1)=nnz(ismember(NeuronSet.reg(ismember(NeuronSet.id,NeuronSet_nonsel.id)),reg_index{r,1}))...
%         /nnz(ismember(NeuronSet_nonsel.reg,reg_index{r,1}));
% end


% load('F:\pixel-dualtask\pev-reactive-correct-error.mat')
% for i=1:2912
%     c=pev{1}(i,~isnan(pev{1}(i,:)));
%     e=pev{2}(i,~isnan(pev{2}(i,:)));
%     if length(c)<250||length(e)<250
%         p(i)=100;
%     else
%         %         p(i)=ranksum(c,e,'tail','right');
%         %         p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
%         if mean(c)>mean(e)
%             p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
%         else
%             p(i)=100;
%         end
%     end    
% end
% NeuronSet1=structfun(@(x)x(p<0.05),NeuronSet,'UniformOutput',false);
% 
% reg=unique(NeuronSet1.reg);
% reg(cellfun(@(x)nnz(ismember(NeuronSet.reg,x))<5,reg))=[];
% for r=1:length(reg)
%     fra(r,1)=nnz(ismember(NeuronSet1.reg,reg{r,1}))/nnz(ismember(NeuronSet.reg,reg{r,1}));
% end
% 
% map_cells{1}=containers.Map(reg,fra); %{ccfid:proj_dense_mat_idx}
end

function map_cells=map_sel_lost_gain(homedir)

% d2=load(fullfile(homedir,'pev_perform_d2.mat'));
% d3=load(fullfile(homedir,'pev_perform.mat'));

load(fullfile(homedir,'difftype_ByPEV.mat'))
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

% lost_all=structfun(@(x)x(ismember(NeuronSet_sel.id,setdiff(union(lost.md.id,lost.ld.id),union(gain.md.id,gain.ld.id)))),NeuronSet_sel,'UniformOutput',false);
% gain_all=structfun(@(x)x(ismember(NeuronSet_sel.id,setdiff(union(gain.md.id,gain.ld.id),union(lost.md.id,lost.ld.id)))),NeuronSet_sel,'UniformOutput',false);
% 
% for i=1:length(d2.NeuronSet.id)
%     c=squeeze(d2.pev{3}(i,1,~isnan(d2.pev{3}(i,1,:))));
%     e=squeeze(d2.pev{4}(i,1,~isnan(d2.pev{4}(i,1,:))));
%     if length(c)<250||length(e)<250
%         p1(i)=100;
%     else
% %         p(i)=ranksum(c,e,'tail','right');
% %         p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
%         if mean(c)>mean(e)
%             p1(i)=anova1([c;e],[zeros(length(c),1);ones(length(e),1)],'off');
%         else
%             p1(i)=100;
%         end
%     end    
% end
% for i=1:length(d3.NeuronSet.id)
%     c=squeeze(d2.pev{3}(i,1,~isnan(d2.pev{3}(i,1,:))));
%     e=squeeze(d2.pev{4}(i,1,~isnan(d2.pev{4}(i,1,:))));
%     if length(c)<250||length(e)<250
%         p2(i)=100;
%     else
% %         p(i)=ranksum(c,e,'tail','right');
% %         p(i)=anova1([c,e],[zeros(1,length(c)),ones(1,length(e))],'off');
%         if mean(c)>mean(e)
%             p2(i)=anova1([c;e],[zeros(length(c),1);ones(length(e),1)],'off');
%         else
%             p2(i)=100;
%         end
%     end    
% end
% d2.NeuronSet1=structfun(@(x)x(p1<0.05),d2.NeuronSet,'UniformOutput',false);
% d3.NeuronSet1=structfun(@(x)x(p2<0.05),d3.NeuronSet,'UniformOutput',false); 
% 
% NeuronSet_d2d3=structfun(@(x)x(ismember(NeuronSet_sel.id,union(d2.NeuronSet1.id,d3.NeuronSet1.id))),NeuronSet_sel,'UniformOutput',false);
% lost_perf=structfun(@(x)x(ismember(NeuronSet_d2d3.id,setdiff(union(lost.md.id,lost.ld.id),union(gain.md.id,gain.ld.id)))),NeuronSet_d2d3,'UniformOutput',false);
% gain_perf=structfun(@(x)x(ismember(NeuronSet_d2d3.id,setdiff(union(gain.md.id,gain.ld.id),union(lost.md.id,lost.ld.id)))),NeuronSet_d2d3,'UniformOutput',false);

% reg=unique(lost_perf.reg);
% reg(cellfun(@(x)nnz(ismember(lost_all.reg,x))<10,reg))=[];
% for r=1:length(reg)
%     fra1(r,1)=nnz(ismember(lost_perf.reg,reg{r,1}))/nnz(ismember(lost_all.reg,reg{r,1}));
% end
% reg=unique(gain_perf.reg);
% reg(cellfun(@(x)nnz(ismember(gain_all.reg,x))<5,reg))=[];
% for r=1:length(reg)
%     fra2(r,1)=nnz(ismember(gain_perf.reg,reg{r,1}))/nnz(ismember(gain_all.reg,reg{r,1}));
% end
map_cells{1}=containers.Map(reg,fra1); %{ccfid:proj_dense_mat_idx}
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