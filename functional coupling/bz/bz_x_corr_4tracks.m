% certain time
clear
clc

if ispc
    CodePath='D:\code-hem';
    HomePath='F:\pixel-dualtask';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end
all=true; % choose all recording time or no
tasktype='dualtask';
%%

sess=unique(floor(h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/cluster_id')/100000));
load(fullfile(HomePath,'session_list.mat'))
% f=dir(fullfile(HomePath,'DataSum','M*'));
% session=struct2cell(f)';
% session=session(:,1);
% sess=(1:length(session))';

% range=[10,54];
for i=sess'    
    tic
%     if isfile([savepath sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,i)])
%         continue
%     end
    folder=session{i};  
    disp(folder) 
   
    if all
        if isfile(fullfile(HomePath,'xcorr','bzdata_20ms',sprintf('BZ_XCORR_duo_f%d.mat',i)))
            continue
        end
        [avail,spkTS,spkID]=pre_process(i,folder,HomePath,tasktype,'true');

         mono=bz.sortSpikeIDz(spkTS,spkID,'20ms');
         
         save(fullfile(HomePath,'xcorr','bzdata_20ms',sprintf('BZ_XCORR_duo_f%d.mat',i)),'mono','folder','-v7.3')
         clear mono spkTS spkID
    else
        [avail,FT_SPIKE,trials]=pre_process(i,folder,HomePath,tasktype,'false');
        trial_type=["distractorNo-correct","distractorNoGo-correct","distractorGo-correct"];
        phase=["delay"];
        %     phase=["d3","delay"];
        for t=1:numel(trial_type)
            [sel_S1,sel_S2]=util.ExtractTrial(trials,'task','dualtask','trials',trial_type(t));
            for p=1:numel(phase)
                if exist(fullfile(HomePath,'xcorr','bz0212',trial_type(t),phase(p),sprintf('BZ_XCORR_duo_f%d.mat',i)),'file')
                    delete(fullfile(HomePath,'xcorr','bz0212',trial_type(t),phase(p),sprintf('BZ_XCORR_duo_f%d.mat',i)))
                end
                spkTS=cell2mat(cellfun(@(x,y)y(ismember(x,sel_S1)),FT_SPIKE.(phase(p)).trial,FT_SPIKE.(phase(p)).timestamp','UniformOutput',false))';
                spkID=cell2mat(cellfun(@(x,y)str2double(y)*ones(nnz(ismember(x,sel_S1)),1),FT_SPIKE.(phase(p)).trial,FT_SPIKE.(phase(p)).label'...
                    ,'UniformOutput',false)');
                mono.s2=bz.sortSpikeIDz(spkTS,spkID);
                
                spkTS=cell2mat(cellfun(@(x,y)y(ismember(x,sel_S2)),FT_SPIKE.(phase(p)).trial,FT_SPIKE.(phase(p)).timestamp','UniformOutput',false))';
                spkID=cell2mat(cellfun(@(x,y)str2double(y)*ones(nnz(ismember(x,sel_S2)),1),FT_SPIKE.(phase(p)).trial,FT_SPIKE.(phase(p)).label'...
                    ,'UniformOutput',false)');
                mono.s2=bz.sortSpikeIDz(spkTS,spkID);
                save(fullfile(HomePath,'xcorr','bzdata_20ms',trial_type(t),phase(p),sprintf('BZ_XCORR_duo_f%d.mat',i)),'mono','folder','-v7.3')
                clear mono spkTS spkID
            end            
        end
    end
end
return

%% Function
function [avail,out1,out2]=pre_process(sessionIdx,Sessionfolder,HomePath,task,allrecording)
addpath(fullfile('D:\code','fieldtrip-20200320'))
addpath(fullfile('D:\code','npy-matlab','npy-matlab'))
ft_defaults
folder=dir(fullfile(HomePath,'DataSum',Sessionfolder,'*','FR_All_250ms.hdf5'));
[cluster_ids,spkTS,spkId]=deal([]);
for f=1:size(folder,1)    
    cluster_ids_temp=h5read(fullfile(folder(f,1).folder,folder(1).name),'/SU_id')+sessionIdx*100000;
    cluster_ids=[cluster_ids;cluster_ids_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
    spkId_temp=double(readNPY(fullfile(folder(f,1).folder,'spike_clusters.npy'))+sessionIdx*100000);
    spkId=[spkId;spkId_temp+10000*str2num(regexp(folder(f).folder,'(?<=imec)(\d)','match','once'))];
    spkTS_temp=double(h5read(fullfile(folder(f,1).folder,'spike_times.hdf5'),'/spkTS'));    
    spkTS=[spkTS;spkTS_temp];
end

if allrecording
    out1=spkTS(ismember(spkId,cluster_ids));
    out2=spkId(ismember(spkId,cluster_ids));  
else
    
    FT_SPIKE=struct();
    FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
    FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
    
    FT_SPIKE.timestamp=arrayfun(@(x)spkTS(spkId==x)',cluster_ids,'UniformOutput',false);
    %  continuous format F T struct file
    trials=util.markLPerf(h5read(fullfile(folder(1,1).folder,folder(1,1).name),'/Trials'),0.7,0.8,120,task);
    if isempty(trials)
        avail=false;
        out1=[];
        out2=[];
        return
    end
    
    sps=30000;
    out1=trials;
    phase={'d1','d2','d3','delay'};
    t=[1,3;4,6;8,9;1,9];
    for i=1:numel(phase)
        cfg=struct();
        cfg.trl=[trials(:,1)+t(i,1)*sps,trials(:,1)+t(i,2)*sps,zeros(size(trials,1),1)+t(i,1)*sps,trials];
        cfg.trlunit='timestamps';
        cfg.timestampspersecond=sps;
        
        FT_SPIKE.(phase{i})=ft_spike_maketrials(cfg,FT_SPIKE);
    end
    out1=FT_SPIKE;
end
avail=true;
end


