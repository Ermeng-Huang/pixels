tic
%% 
% lwd=cd;
cd(lwd)
addpath('D:\code\npy-matlab\npy-matlab\')
addpath(genpath('D:\code\Kilosort3')) % path to kilosort folder
addpath('D:\code\npy-matlab') % for converting to Phy
rootZ = lwd; % the raw data binary file is in this folder
rootH = 'D:\temp'; % path to temporary binary file (same size as data, should be on fast SSD)
main_kilosort3(rootZ,rootH)
toc

%%
tic
addpath('D:\code\neuropixel-utils')
channelMapFile='D:\code\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat';
fn=ls('*.ap.bin');
imec=Neuropixel.ImecDataset(fullfile(lwd,fn),'ChannelMap',channelMapFile);
sync6=imec.readSync();
% save('sync.mat','sync6');
if cleaned
    syncH5=fullfile(lwd,'sync.hdf5');
else
    syncH5=fullfile([lwd,'_cleaned'],'sync.hdf5');
end
h5create(syncH5,'/sync',size(sync6),'Datatype','int8')
h5write(syncH5,'/sync',int8(sync6))

toc

%%
