function out=ChooseNeuronSet(opt)
arguments    
    opt.learning (1,1) char {mustBeMember(opt.learning,{'L','W','all'})};
    opt.region (1,:) char
    opt.regionLevel (1,1) double
    opt.sel (1,:) char
    opt.opto (1,:) char
end

Path=util.Path_default;
sus_trans=h5read(Path.selectivity,'/sus_trans_noPermutaion')';
cluster_id=h5read(Path.selectivity,'/cluster_id');
reg=regexp(h5read(Path.selectivity,'/reg'),'(\w|\\|-)*','match','once')';
path=regexp(h5read(Path.selectivity,'/path'),'(\w|\\|-)*','match','once');
% phase
if strcmp(opt.learning,'L')    
    load(Path.performance)
    phase=PerfList(:,2)==1;
elseif strcmp(opt.learning,'W')
    load(Path.performance)
    phase=PerfList(:,3)==1;
else
    phase=ones(size(cluster_id));
end 
%region
if strcmp(opt.region,'all')    
    reg_good=ismember(reg(:,3),{'CH','BS'}) & ~cell2mat(cellfun(@(x)isempty(x),reg(:,opt.regionLevel),'UniformOutput',false));    
%     load(Path.region)
%     if strcmp(opt.learning,'L') 
%         reg_good=ismember(reg,reg_L);
%     else strcmp(opt.learning,'W')
%         reg_good=ismember(reg,reg_W);
%     end
else
    
end
if strcmp(opt.sel,'sustained')
    sel=sus_trans(:,1)==1;
elseif strcmp(opt.sel,'transient')
    sel=sus_trans(:,2)==1;
elseif strcmp(opt.sel,'sel') %not include switch neuron
    sel=sus_trans(:,1)==1|sus_trans(:,2)==1;
elseif strcmp(opt.sel,'nonsel')
    sel=sus_trans(:,1)==0& sus_trans(:,2)==0 & sus_trans(:,3)==0;
else
    sel=ones(size(cluster_id));
end

if strcmp(opt.opto,'none')
    laser_good=ones(size(cluster_id));
elseif strcmp(opt.opto,'laser-nonmodulated')
    laser=h5read(Path.opto,'/rank');
    laser_good=any(laser(:,17:20),2);
elseif contains(opt.opto,'laser')    %laser-0/1/2/3/4
%     load(Path.opto) %%%%%%%%%
%     p=cell2mat(cellfun(@(x)(x(2,:)),P(:,1),'UniformOutput',false));
%     n=str2double(regexp(opt.opto,'(?<=laser-)(\d)','match','once'));    
%     if opt.opto0
%         laser_good=p(:,17+n)<0.05;
%     else
%         laser_good=p(:,17+n)<0.05 & p(:,17)>0.05;
%     end
end

out.id=cluster_id(phase & sel & laser_good & reg_good);
out.path=path(phase & sel & laser_good & reg_good);
out.reg=reg(phase & sel & laser_good & reg_good);
end