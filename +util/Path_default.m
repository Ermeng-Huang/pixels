function path=Path_default

[~,path.home]=util.ComputerPath();

path.selectivity=fullfile(path.home,'Selectivity_AIopto_0419.hdf5');
path.performance=fullfile(path.home,'Performance_allneuron_0308.mat');
path.region=fullfile(path.home,'reg_keep_40_0412.mat'); %% good region (cells > 40)
path.opto=fullfile(path.home,'LaserModulation.hdf5');
end