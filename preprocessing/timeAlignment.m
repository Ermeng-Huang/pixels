clear
clc
cd('D:\pixel-optogenetic')
% cd('D:\neupix\hem')
addpath('D:\code\npy-matlab\npy-matlab')
load('D:\pixel-optogenetic\session_list.mat')

for i=1:length(session)    
    file=dir(fullfile('D:\pixel-optogenetic\DataSum\',session{i,1},'*\events.hdf5'));
    for j=1:size(file,1)
        if exist(fullfile(file(j).folder,'spike_times.hdf5'),'file')
            delete(fullfile(file(j).folder,'spike_times.hdf5'))
        end
        disp(session{i,1})
        if j==1
            trials_reference=h5read(fullfile(file(j).folder,'events.hdf5'),'/trials');            
            spkTS_reference=readNPY(fullfile(file(j).folder,'spike_times.npy'));            
            h5create(fullfile(file(j).folder,'spike_times.hdf5'),'/spkTS',size(spkTS_reference),'Datatype','uint64')
            h5write(fullfile(file(j).folder,'spike_times.hdf5'),'/spkTS',spkTS_reference)            
        else
            trials=h5read(fullfile(file(j).folder,'events.hdf5'),'/trials');
            spkTS=readNPY(fullfile(file(j).folder,'spike_times.npy'));
            spkTS_new=[];
            for t=1:size(trials,2)
                if t == 1
                  spkTS_new=spkTS(spkTS<trials(1,t+1))-double((trials(1,t)-trials_reference(1,t))); 
                elseif t < size(trials,2)
                    spkTS_new=[spkTS_new;spkTS(spkTS<trials(1,t+1)&spkTS>=trials(1,t))-double((trials(1,t)-trials_reference(1,t)))]; 
                else
                    spkTS_new=[spkTS_new;spkTS(spkTS>=trials(1,t))-double((trials(1,t)-trials_reference(1,t)))]; 
                end
            end
            h5create(fullfile(file(j).folder,'spike_times.hdf5'),'/spkTS',size(spkTS_new),'Datatype','uint64')
            h5write(fullfile(file(j).folder,'spike_times.hdf5'),'/spkTS',spkTS_new)
            clear trials spkTS
        end
    end
    clearvars -Except session
end


