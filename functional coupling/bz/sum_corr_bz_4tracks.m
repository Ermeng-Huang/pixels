% change data format for calculated zx code
clc
clear
close all

if ispc
    CodePath='D:\code';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
addpath('D:\code\pixels-master\jpsth')

opt.type='AIOPTO';
opt.data=[];
opt.prefix= '0504';
opt.poolsize=1;
opt.type='sums'; 
opt.dtype='AIOPTO';

if ~isfile(fullfile(HomePath,'bz0313_sums_conn.mat'))
    load(fullfile(HomePath,'conn_bz_0314.mat'))
    load(fullfile(HomePath,'session_list.mat'))
    for i=1:length(session)
        sums_conn_str(i).sig_con=sig_con{i};
        sums_conn_str(i).folder=session{i};
    end
    save(fullfile(HomePath,'bz0313_sums_conn.mat'),'sums_conn_str')

end

reg_conn_bz(opt,HomePath)
%% Function

function reg_conn_bz(opt,HomePath)
%TODO merge with load_sig_pair script
if isempty(opt.data)
    if strcmp(opt.type,'neupix')
        load('sums_conn.mat','sums_conn_str')
    else
        load(fullfile(HomePath,'bz0313_sums_conn.mat'),'sums_conn_str')
    end
else
    sums_conn_str=opt.data;
end

for fidx=1:length(sums_conn_str)
    disp(fidx);
    fpath=sums_conn_str(fidx).folder; %session data folder
    if strcmp(opt.type,'neupix')
        pc_stem=replace(regexp(fpath,'(?<=SPKINFO/).*$','match','once'),'/','\');
        inputf=fullfile('K:','neupix','SPKINFO',pc_stem,'FR_All_1000.hdf5');
        all_su=int32(h5read(inputf,'/SU_id'));
    else
        pc_stem=fpath;
        inputf=fullfile('D:','pixel-optogenetic','DataSum',fpath,'FT_SPIKE.mat');
        fstr=load(inputf);
        all_su=int32(cellfun(@(x) str2double(x),fstr.FT_SPIKE.label));
    end
    
    sig_con=int32(sums_conn_str(fidx).sig_con); % significant functional coupling
    pair_comb_one_dir=nchoosek(all_su,2); % all pairs combination
    [sig_meta,pair_meta]=bz.util.get_meta(sig_con,pair_comb_one_dir,pc_stem,opt); % assign meta info
    
%     fields={'suid','reg','wrsp','selec','mem_type'};
    fields={'suid','reg','mem_type'};
    for fi=fields
        %TODO online genenrate session tag
%         sig.(fi{1})=cat(1,sig.(fi{1}),sig_meta.(fi{1}));
        pair_meta.(fi{1})=cat(1,pair_meta.(fi{1}),flip(pair_meta.(fi{1}),ndims(pair_meta.(fi{1}))));%uni-dir to bi-dir
%         pair.(fi{1})=cat(1,pair.(fi{1}),pair_meta.(fi{1}));
    end
    tic
    save(fullfile('D:','pixel-optogenetic','xcorr','bzdata',sprintf('%s_conn_w_reg_%d.mat',opt.prefix,fidx)),'sig_meta','pair_meta','pc_stem','-v7.3','-nocompression')
    toc
end
end
