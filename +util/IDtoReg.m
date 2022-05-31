function out = IDtoReg(homedir,in)
% extract region information from neuron ID
cluster_id=h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/cluster_id');
reg=regexp(h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/reg'),'(\w|\\|-)*','match','once');
if iscell(in)
    all_su=unique(cell2mat(in(:)));
else
    all_su=unique(in(:));
end
for i=1:numel(all_su)
    r(i,:)=reg(cluster_id==all_su(i),:);
end
out=arrayfun(@(x)r(all_su==x,:),in,'UniformOutput',false);
end
