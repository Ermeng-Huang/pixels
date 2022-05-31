function path=Path_default(opt)
arguments   
    opt.task (1,:) char        
end
if ispc
    switch (getenv('USERNAME'))
        case 'hem'
            if strcmp(task,'AIopto')
            CodePath='D:\code\AIopto';
            HomePath='I:\pixel-optogenetic';
            else
                
            end
        case 'zx'
            CodePath='D:\code-hem';
            if strcmp(opt.task,'AIopto')
                HomePath='D:\pixel-optogenetic';
            elseif strcmp(opt.task,'dualtask')
                HomePath='F:\pixel-dualtask';
            end
    end
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
if strcmp(opt.task,'AIopto')
    path.selectivity=fullfile(HomePath,'Selectivity_AIopto_0419.hdf5');
    path.performance=fullfile(HomePath,'Performance_allneuron_0308.mat');
    path.region=fullfile(HomePath,'reg_keep_40_0412.mat'); %% good region (cells > 40)
    path.opto=fullfile(HomePath,'LaserModulation.hdf5');
elseif strcmp(opt.task,'dualtask')
   path.selectivity=fullfile(HomePath,'Selectivity_1129.hdf5');  
end
end