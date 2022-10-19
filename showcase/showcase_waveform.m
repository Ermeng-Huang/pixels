% extract waveform of certain neuron   
clc
clear
close all
addpath('D:\code\neuropixel-utils')
addpath(genpath('D:\code\Kilosort2'))
channelMapFile='D:\code\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat';
data_type='hem';
cidx=2700481;
if strcmp(data_type,'zx')
    homedir='F:\neupix\raw';
    filepath='F:\neupix\raw';
    folder='191204_49_learning8_g0_imec1_cleaned';
   
else
    homedir='F:\pixel-dualtask\DataSum';
    path=regexp(h5read('F:\pixel-dualtask\Selectivity_1129.hdf5','/path'),'(\w|\\|-)*','match','once');
    cluster_id=h5read('F:\pixel-dualtask\Selectivity_1129.hdf5','/cluster_id');
    folder=path{cluster_id==cidx};
    filepath='H:\Neuropixels';
end

fl=dir(fullfile(filepath,folder,'cluster_info.tsv'));
% savepath=fullfile('F:\pixel-dualtask\DataSum',folder);
savepath=fullfile(homedir,folder);
if isfile(fullfile(savepath,sprintf('waveform_%d.mat',cidx)))
    load(fullfile(savepath,sprintf('waveform_%d.mat',cidx)),'waveform');
else    
    waveform=cell(0,5);
    ks=Neuropixel.KiloSortDataset(fl.folder,'channelMap',channelMapFile);
    ks.load();    

    snippetSetTet = ks.getWaveformsFromRawData('cluster_ids',rem(cidx,10000),'num_waveforms', 100, 'best_n_channels', 50, 'car', true, ...
        'subtractOtherClusters', false,'window', [-30 60]);
    
    snippetSetBest = ks.getWaveformsFromRawData('cluster_ids',rem(cidx,10000),'num_waveforms', 500, 'best_n_channels', 1, 'car', true, ...
        'subtractOtherClusters', false,'window', [-30 60]);
    
    waveform{end+1,1}=folder;
    waveform{end,2}=cidx;
    waveform{end,3}=snippetSetTet.data;
    waveform{end,4}=mean(snippetSetBest.data,3);
    waveform{end,5}=snippetSetTet.channel_ids_by_cluster;

    save(fullfile(savepath,sprintf('waveform_%d.mat',cidx)),'waveform');
end
%% plot waveform
channel=waveform{1,5};
waveform=waveform{1,3};
f=figure('Color','w','Position',[100,100,80,800]);
channel_best=double(channel(1));
channel_plot=channel_best+(-6:1:6);
channel_plot=channel_plot(ismember(channel_plot,channel));
for i=1:9
    subplot(9,1,i)
    plot(mean(waveform(channel==channel_plot(find(channel_plot==channel_best)-5+i),:,:),3),'k-','LineWidth',2)
    set(gca,'XTick',1:30:91,'XTickLabel',{'-1','0','1','2'},'YLim',[-30 15])
    xlabel('time (ms)')
    ylabel('MicroVolt')
    title(sprintf('channel %d',channel(1)+i-5))
    box off
    if i==9
        hold on
        plot([1,31],[-20,-20],'k-','LineWidth',2)
        plot([1,1],[1,-20],'k-','LineWidth',2)
    end
end

exportgraphics(f,fullfile('F:\pixel-dualtask\showcase\waveform',sprintf('waveform-%d.pdf',cidx)))
% exportgraphics(f,fullfile(savepath,sprintf('waveform-%d.pdf',cidx)));