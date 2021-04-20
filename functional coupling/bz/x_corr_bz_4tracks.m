clear
clc

if ispc
    CodePath='D:\code';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end

addpath(fullfile(CodePath,'npy-matlab','npy-matlab'))
addpath(fullfile(CodePath,'pixels-bleeding','jpsth'))
addpath(genpath(fullfile(CodePath,'buzcode')))

prefix='0313';
load(fullfile(HomePath,'session_list.mat'))
range=[21,82]
for i=range(1):range(2)    
    tic
%     if isfile([savepath sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,i)])
%         continue
%     end
    folder=session{i};  
    disp(folder) 
    if ~isfile(fullfile(HomePath,'xcorr',sprintf('bz%s',prefix),sprintf('BZ_XCORR_duo_f%d.mat',i)))
        [spkTS,spkID]=pre_process(i,folder,HomePath);
        mono=bz.sortSpikeIDz(spkTS,spkID);
        save(fullfile(HomePath,'xcorr',sprintf('bz%s',prefix),sprintf('BZ_XCORR_duo_f%d.mat',i)),'mono','-v7.3','folder')
    else
        load(fullfile(HomePath,'xcorr',sprintf('bz%s',prefix),sprintf('BZ_XCORR_duo_f%d.mat',i)),'mono')
    end
    for n=1:size(mono.sig_con,1)
        sig_con{i,1}(n,:)=[mono.completeIndex(mono.sig_con(n,1),2),mono.completeIndex(mono.sig_con(n,2),2)];
    end
    toc
    clear mono
end

save(fullfile(HomePath,sprintf('conn_bz_%d_%d.mat',range(1),range(2))),'sig_con','range')

return

%% Function
function [spkTS,spkId]=pre_process(sessionIdx,Sessionfolder,HomePath)
folder=dir(fullfile(HomePath,'DataSum',Sessionfolder,'*','FR_All.hdf5'));
trial=h5read(fullfile(folder(1,1).folder,'FR_All.hdf5'),'/Trials');
cluster_ids=[];
for f=1:size(folder,1)
    cluster_ids_temp=[];
    cluster_ids_temp=h5read(fullfile(folder(f,1).folder,'FR_All.hdf5'),'/SU_id')+sessionIdx*100000;
    cluster_ids=[cluster_ids;cluster_ids_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];   
end

spkTS=[];
spkId=[];
for f=1:size(folder,1)
    spkId_temp=[];
    spkId_temp=double(readNPY(fullfile(folder(f,1).folder,'spike_clusters.npy'))+sessionIdx*100000);
    spkId=[spkId;spkId_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
    
    spkTS_temp=[];
    spkTS_temp=double(h5read(fullfile(folder(f,1).folder,'spike_times.hdf5'),'/spkTS'));
    n=folder(f,1).folder;
    spkTS=[spkTS;spkTS_temp];
end
spk_bad=~ismember(spkId,cluster_ids);
spkTS(spk_bad)=[];
spkId(spk_bad)=[];

end


