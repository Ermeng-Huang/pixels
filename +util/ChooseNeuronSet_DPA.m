function out=ChooseNeuronSet_DPA(homedir,opt)
arguments    
    homedir (1,:) char = 'F:/pixel-DPA8s';
    opt.region (1,:) char = 'all'
    opt.regionLevel (1,1) double = 7
    opt.sel (1,:) char = 'all'
    opt.opto (1,:) char ='none'
    opt.hemisphere (1,:) char ='all'
    opt.learning (1,:) char='none'
    opt.selfile (1,:) char='0327'
    opt.trial (1,:) char = 'distractorNo'
    opt.pref (1,:) logical = false
end
if strcmp(opt.selfile,'0327')
%     reg=regexp(h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/reg'),'(\w|\\|-)*','match','once');
    cluster_id=h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/cluster_id');
    path=regexp(h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/path'),'(\w|\\|-)*','match','once');
    sust_trans=h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/sust_trans_noPermutaion');
    if strcmp(opt.sel,'sel')
        sel=any(sust_trans(:,1:2),2);
    elseif strcmp(opt.sel,'d1')
        sel=any(sust_trans(:,1:2),2)&any(sust_trans(:,5:6),2);
    elseif strcmp(opt.sel,'d2')        
        sel=any(sust_trans(:,1:2),2)&any(sust_trans(:,8:9),2);
    elseif strcmp(opt.sel,'d3')
        sel=any(sust_trans(:,1:2),2)&any(sust_trans(:,12),2);
    elseif strcmp(opt.sel,'mixed-sel')
        distractor=h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5','1119')),'/sust_trans_noPermutaion');
        sel=any(sust_trans(:,1:2),2)&any(distractor(:,1:2),2); 
    elseif strcmp(opt.sel,'nonmixed-sel')
        distractor=h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5','1119')),'/sust_trans_noPermutaion');
        sel=any(sust_trans(:,1:2),2)&~any(distractor(:,1:2),2); 
   
    else
        sel=ones(size(cluster_id));
    end
   
%     if strcmp(opt.sel,'d1d2')
%         sel=any(abs(Pvalue(:,3:4))<0.05/2&abs(Pvalue(:,6:7))<0.05/2,2);
%     elseif strcmp(opt.sel,'d1')
%         sel=any(abs(Pvalue(:,3:4))<0.05/2&~(abs(Pvalue(:,6:7))<0.05/2),2);
%     elseif strcmp(opt.sel,'d2')
%         sel=any(~(abs(Pvalue(:,3:4))<0.05/2)&abs(Pvalue(:,6:7))<0.05/2,2);
%     elseif strcmp(opt.sel,'transinet_end')
%         sel=sust_trans(:,12)~=0;
%     elseif strcmp(opt.sel,'d1andd2')
%         sel=any(abs(Pvalue(:,3:4))<0.05/2|abs(Pvalue(:,6:7))<0.05/2,2);
%     end
elseif strcmp(opt.selfile,'1108')   
    reg=regexp(h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/reg'),'(\w|\\|-)*','match','once');
    cluster_id=h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/cluster_id');
    path=regexp(h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/path'),'(\w|\\|-)*','match','once');
    sust_trans=h5read(fullfile(homedir,sprintf('Selectivity_%s.hdf5',opt.selfile)),'/sust_trans_noPermutaion');
    if strcmp(opt.sel,'d1d2')
        
    elseif strcmp(opt.sel,'s2')       
        
    elseif strcmp(opt.sel,'d')
        sel=any(sust_trans(:,1:2),2);
    elseif strcmp(opt.sel,'d1')
        sel=any(sust_trans(:,1:2),2)&(any(sust_trans(:,5:6),2));
    elseif strcmp(opt.sel,'d2')        
        sel=any(sust_trans(:,1:2),2)&(any(sust_trans(:,8:9),2));
    elseif strcmp(opt.sel,'d3')
        sel=any(sust_trans(:,1:2),2)&(any(sust_trans(:,12),2));
    end

elseif strcmp(opt.selfile,'0925')
    Pvalue=h5read(fullfile(homedir,'Selectivity_0925.hdf5'),sprintf('/Pvalue_%s',opt.trial));
    
    reg=regexp(h5read(fullfile(homedir,'Selectivity_0925.hdf5'),'/reg'),'(\w|\\|-)*','match','once');
    cluster_id=h5read(fullfile(homedir,'Selectivity_0925.hdf5'),'/cluster_id');
    path=regexp(h5read(fullfile(homedir,'Selectivity_0925.hdf5'),'/path'),'(\w|\\|-)*','match','once');
    Pvalue_distractor=h5read(fullfile(homedir,'Selectivity_0925.hdf5'),'/Pvalue_distractor');
    ds_s2=abs(Pvalue_distractor(:,5))<0.05;
    ds_d2=any(abs(Pvalue_distractor(:,6:7))<0.05/2,2);
    ds_d3=abs(Pvalue_distractor(:,10))<0.05; 
    
    d1=any(abs(Pvalue(:,3:4))<0.05/2,2);
    d2=any(abs(Pvalue(:,6:7))<0.05/2,2);
    d3=any(abs(Pvalue(:,10))<0.05/2,2);
    if strcmp(opt.sel,'d1')
        sel=d1;
%         sel=d1&~d2&~d3;
    elseif strcmp(opt.sel,'d2')
        sel=d2;
%         sel=d2&~d1&~d3;
    elseif strcmp(opt.sel,'d3')
        sel=d3&ds_d3;
%         sel=d3&~d1&~d2;
    elseif strcmp(opt.sel,'sel')
        sel=d1|d2|d3;
    else
        sel=ones(size(cluster_id));
    end
    
end

%%% region
% if strcmp(opt.region,'all')    
%     reg_good=ismember(reg(:,3),{'CH','BS'}) & ~cell2mat(cellfun(@(x)isempty(x),reg(:,opt.regionLevel),'UniformOutput',false));    
% else
%    reg_good=ismember(reg(:,7),{opt.region});    
% end
% if strcmp(opt.hemisphere,'all') 
% else
%    reg_good=reg_good&ismember(reg(:,11),opt.hemisphere);  
% end
if ~exist('reg_good','var')
    reg_good=ones(size(cluster_id));
end
%%% phase
if strcmp(opt.learning,'L')    
    load(Path.performance)
    phase=PerfList(:,2)==1;
elseif strcmp(opt.learning,'W')
    load(Path.performance)
    phase=PerfList(:,3)==1;
else
    phase=ones(size(cluster_id));
end

%%% laser
if strcmp(opt.opto,'none')
    laser_good=ones(size(cluster_id));
elseif strcmp(opt.opto,'laser-modulated')
    laser=h5read(Path.opto,'/rank');
    laser_good=any(laser(:,17),2);
elseif strcmp(opt.opto,'laser-nonmodulated')
    laser=h5read(Path.opto,'/rank');
    laser_good=~any(laser(:,17:20),2);
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
% out.reg=reg(phase & sel & laser_good & reg_good,opt.regionLevel);

if opt.pref
   out.pref(any(sust_trans(ismember(cluster_id,out.id),5:end)==1,2),1)=1;
   out.pref(any(sust_trans(ismember(cluster_id,out.id),5:end)==2,2),1)=2;
end
end