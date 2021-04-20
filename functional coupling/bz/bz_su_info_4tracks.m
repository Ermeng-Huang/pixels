clear
CodePath='D:\code';
HomePath='D:\pixel-optogenetic';

%% Summary of results in xcorr(bz)
if ~isfile(fullfile(HomePath,'conn_bz_0314.mat'))
    f=dir(fullfile(HomePath,'conn_bz_*_*.mat'));
    sig_con=cell(0);
    
    for i=1:length(f)
        temp=load(fullfile(HomePath,f(i,1).name));
        sig_con(temp.range(1):temp.range(2),:)=temp.sig_con(temp.range(1):temp.range(2),:);
    end
    
    save(fullfile(HomePath,'conn_bz_0314.mat'),'sig_con')
else
    load(fullfile(HomePath,'conn_bz_0314.mat'))
end

%%
%TODO merge with load_sig_pair script
load(fullfile(HomePath,'session_list.mat'))
Path=util.Path_default;
sus_trans=h5read(Path.selectivity,'/sus_trans_noPermutaion')';
cluster_id=h5read(Path.selectivity,'/cluster_id');
reg=regexp(h5read(Path.selectivity,'/reg'),'(\w|\\|-)*','match','once');
path=regexp(h5read(Path.selectivity,'/path'),'(\w|\\|-)*','match','once');

for fidx=1:length(sig_con)
    disp(fidx);    
    sig_con{fidx,2}.pc_stem=session{fidx};      
    sig_con{fidx,2}.all_su=cluster_id(contains(path,sig_con{fidx,2}.pc_stem));   
    sig_con{fidx,2}.pair_comb=nchoosek(sig_con{fidx,2}.all_su,2); % all pairs combination    
    sig_con{fidx,2}.all_reg=reg(:,contains(path,sig_con{fidx,2}.pc_stem))';
    sig_con{fidx,2}.all_sel=sus_trans(contains(path,sig_con{fidx,2}.pc_stem),:); 
end

save(fullfile(HomePath,'states_conn_bz_0314.mat'),'sig_con','-v7.3','-nocompression')



function [sig,pair]=get_meta(sig_id,pair_id_one_dir,fpath)

idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
meta_str=ephys.util.get_meta();
reg_map=containers.Map('KeyType','int32','ValueType','any'); %reg_map(su_id)=reg
wrsp_map=containers.Map('KeyType','int32','ValueType','any'); 
selec_map=containers.Map('KeyType','int32','ValueType','any');
mem_type_map=containers.Map('KeyType','int32','ValueType','int32');
% pc_stem=replace(regexp(fpath,'(?<=SPKINFO/).*$','match','once'),'/','\');
sess_idx=find(startsWith(meta_str.allpath,fpath));
for suidx=reshape(sess_idx,1,[])
    suid=meta_str.allcid(suidx);
    acrontree=meta_str.reg_tree(:,suidx);
    ccfid=nan(1,6);
    for i=1:numel(acrontree), if isempty(acrontree{i}), ccfid(i)=0;else, ccfid(i)=idmap.reg2ccfid(acrontree{i});end;end
    reg_map(suid)=int32(ccfid);
    wrsp_map(suid)=meta_str.wrs_p(:,suidx);
    selec_map(suid)=meta_str.selec(:,suidx);
    mem_type_map(suid)=meta_str.mem_type(suidx);
end
sig=struct();
sig.suid=sig_id;
sig.reg=reshape(cell2mat(arrayfun(@(x) reg_map(x),sig_id,'UniformOutput',false)),[],6,2);
sig.wrsp=reshape(cell2mat(arrayfun(@(x) wrsp_map(x),sig_id,'UniformOutput',false)),[],14,2);
sig.selec=reshape(cell2mat(arrayfun(@(x) selec_map(x),sig_id,'UniformOutput',false)),[],14,2);
sig.mem_type=reshape(cell2mat(arrayfun(@(x) mem_type_map(x),sig_id,'UniformOutput',false)),[],2);

pair=struct();
pair.suid=pair_id_one_dir;
pair.reg=reshape(cell2mat(arrayfun(@(x) reg_map(x),pair_id_one_dir,'UniformOutput',false)),[],6,2);
pair.wrsp=reshape(cell2mat(arrayfun(@(x) wrsp_map(x),pair_id_one_dir,'UniformOutput',false)),[],14,2);
pair.selec=reshape(cell2mat(arrayfun(@(x) selec_map(x),pair_id_one_dir,'UniformOutput',false)),[],14,2);
pair.mem_type=reshape(cell2mat(arrayfun(@(x) mem_type_map(x),pair_id_one_dir,'UniformOutput',false)),[],2);
end
