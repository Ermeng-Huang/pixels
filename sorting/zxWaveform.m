tic
addpath('D:\code\neuropixel-utils')
addpath('D:\code\npy-matlab\npy-matlab\')
addpath(genpath('D:\code\Kilosort2'))
channelMapFile='D:\code\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat';
% lwd=pwd();
% lwd='D:\Data\20200617_g0\20200617_g0_imec0';
% lwd=cd;
cd(lwd);
%% waveform
tic
try
    metaf=ls(fullfile(lwd,'*.ap.meta'));
    s1s=30000;
    FR_Th=1;
    fh=fopen(fullfile(lwd,metaf));
    ts=textscan(fh,'%s','Delimiter',{'\n'});
    nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
    spkNThresh=nSample/385/s1s/2*FR_Th;
    clusterInfo = readtable(fullfile(lwd,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
    waveformGood=strcmp(clusterInfo{:,4},'good');
    freqGood=clusterInfo{:,10}>spkNThresh;
    cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
    
    if numel(cluster_ids)>0
        waveform=cell(0,4);
        ks=Neuropixel.KiloSortDataset(lwd,'channelMap',channelMapFile);
        ks.load();
        
        for cidx=1:numel(cluster_ids)
            try
                snippetSetTet = ks.getWaveformsFromRawData('cluster_ids',cluster_ids(cidx),'num_waveforms', 100, 'best_n_channels', 4, 'car', true, ...
                    'subtractOtherClusters', false,'window', [-30 60]);
                
                snippetSetBest = ks.getWaveformsFromRawData('cluster_ids',cluster_ids(cidx),'num_waveforms', 500, 'best_n_channels', 1, 'car', true, ...
                    'subtractOtherClusters', false,'window', [-30 60]);
                
                waveform{end+1,1}=lwd;
                waveform{end,2}=cluster_ids(cidx);
                waveform{end,3}=snippetSetTet.data;
                waveform{end,4}=mean(snippetSetBest.data,3);
            catch ME
                
            end
            %                 if to_plot
            %                     fh=figure();
            %                     plot(cell2mat(arrayfun(@(x) mean(squeeze(snippetSet.data(x,:,:))'), 1:4,'UniformOutput',false)));
            %                     pause
            %                     close(fh);
            %                 end
        end
    end
    save([lwd '\waveform.mat'],'waveform');
    
catch ME
    
end

toc

