% change data format for calculated zx code
clc
clear
close all

if ispc
    CodePath='D:\code';
    HomePath='F:\pixel-dualtask';
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
addpath('D:\code\pixels-master\jpsth')

opt.type='dualtask';
opt.data=[];
opt.poolsize=1;

reg_conn_bz(opt,HomePath)
%% Function

function reg_conn_bz(opt,HomePath)

   id=unique(floor(h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/cluster_id')/100000));
   load(fullfile(HomePath,'session_list.mat'))
   
    for fidx=id'
        disp(fidx);
        
%         if isfile(fullfile(HomePath,'xcorr','bzdata','bzdata_ZX',sprintf('%s_conn_w_reg_%d.mat',opt.prefix,fidx)))
%             continue
%         end            
        fstr=load(fullfile(HomePath,'xcorr','bzdata',sprintf('BZ_XCORR_inh_f%d.mat',fidx)));
        pc_stem=session{fidx};
 
        sig_con=fstr.inh_conn; % significant functional coupling
        
        if strcmp(opt.type,'dualtask')
            [sig_meta,pair_meta]=get_meta(sig_con,pc_stem); % assign meta info
        else
            pair_comb_one_dir=nchoosek(all_su,2); % all pairs combination
            [sig_meta,pair_meta]=bz.util.get_meta(sig_con,pair_comb_one_dir,pc_stem,opt); % assign meta info
        end
        
        
        fields={'suid','reg','mem_type'};
        for fi=fields
            pair_meta.(fi{1})=cat(1,pair_meta.(fi{1}),flip(pair_meta.(fi{1}),ndims(pair_meta.(fi{1}))));%uni-dir to bi-dir
        end        
        save(fullfile(HomePath,'xcorr','bzdata','bzdata_ZX',sprintf('inh_conn_w_reg_%d.mat',fidx)),'sig_meta','pair_meta','pc_stem','-v7.3','-nocompression')
        
    end

end

function [sig_meta,pair_meta]=get_meta(sig_con,pc_stem)
Path=util.Path_default('task','dualtask');
sus_trans=h5read(Path.selectivity,'/sust_trans_noPermutaion');
cluster_id=h5read(Path.selectivity,'/cluster_id');
reg=regexp(h5read(Path.selectivity,'/reg'),'(\w|\\|-)*','match','once');
path=regexp(h5read(Path.selectivity,'/path'),'(\w|\\|-)*','match','once');
all_su=cluster_id(contains(path,pc_stem));

pair_comb_one_dir=nchoosek(all_su,2);

for i=1:numel(all_su)
    if ~any(sus_trans(cluster_id==all_su(i),1:2)) %non-mem
        mem(i)=-1;
    else
        if sus_trans(cluster_id==all_su(i),1)==1  %sustained
            prefer=unique(sus_trans(cluster_id==all_su(i),5:12));
            if prefer(prefer~=0)==1  %s1
                mem(i)=1;
            else
                mem(i)=4;
            end
        elseif sus_trans(cluster_id==all_su(i),2)==1 || sus_trans(cluster_id==all_su(i),2)==3
            prefer=unique(sus_trans(cluster_id==all_su(i),5:12));
            if prefer(prefer~=0)==1 %s1
                mem(i)=2;
            else
                mem(i)=5;
            end
        elseif sus_trans(cluster_id==all_su(i),2)==2
            prefer=unique(sus_trans(cluster_id==all_su(i),13:20));
            if prefer(prefer~=0)==1 %s1
                mem(i)=3;
            else
                mem(i)=6;
            end            
        end
    end
end
reg=cellfun(@(x)Reg2Ccfid(x),reg(ismember(cluster_id,all_su),3:8));

sig_meta.suid=int32(sig_con);
pair_meta.suid=int32(pair_comb_one_dir);
sig_meta.mem_type=int32(arrayfun(@(x)mem(all_su==x),sig_con));
pair_meta.mem_type=int32(arrayfun(@(x)mem(all_su==x),pair_comb_one_dir));
sig_meta.reg=int32(cell2mat(arrayfun(@(x)reg(all_su==x,:),sig_con,'UniformOutput',false)));
sig_meta.reg=cat(3,sig_meta.reg(:,1:end/2),sig_meta.reg(:,end/2+1:end));
pair_meta.reg=int32(cell2mat(arrayfun(@(x)reg(all_su==x,:),pair_comb_one_dir,'UniformOutput',false)));
pair_meta.reg=cat(3,pair_meta.reg(:,1:end/2),pair_meta.reg(:,end/2+1:end));
end

function out=Reg2Ccfid(in)
load('D:\code-hem\reg_ccfid_map','reg2ccfid')

if isempty(in)
    out=0;
else
    out=reg2ccfid(in);
end

end