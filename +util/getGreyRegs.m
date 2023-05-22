function grey_regs_=getGreyRegs(opt)
arguments
    opt.mincount (1,1) double = 100
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.reg_correct (1,1) logical = true
    opt.homedir(1,:) char ='F:\pixel-dualtask';
end
persistent grey_regs opt_
if isempty(grey_regs) || ~isequaln(opt_,opt)    
    if opt.reg_correct
        load(fullfile(opt.homedir,'regs_w2g.mat'),'reg_new');
        reg=reg_new;
    else
        reg=h5read(fullfile(opt.homedir,'Selectivity_1129.hdf5'),'/reg');
    end
    meta.reg_tree=reg(:,3:7)';
    
    BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
    CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
    CTXsel=strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),'');
    switch opt.range
        case 'grey'
            grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));
        case 'CH'
            grey_regs=unique(meta.reg_tree(5,CHsel));
        case 'CTX'            
            grey_regs=unique(meta.reg_tree(5,CTXsel));
    end
    cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
    grey_regs=grey_regs(cnt>opt.mincount);
end
grey_regs_=grey_regs;
opt_=opt;
end