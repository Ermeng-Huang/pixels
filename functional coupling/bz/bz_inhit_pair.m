clear
homedir='F:\pixel-dualtask';
file=dir(fullfile(homedir,'xcorr','bz0902','BZ_XCORR_duo_f*.mat'));
for i=1:size(file,1)
    load(fullfile(homedir,'xcorr','bz0902',file(i,1).name));
    binSize=mono.binSize;
    alpha=0.05;
    sig_con=[];
    
    for refcellID=1:numel(mono.n)
        
        for cell2ID=1:numel(mono.n)
            if refcellID==cell2ID
                continue
            end
            cch=mono.ccgR(:,refcellID,cell2ID);
            hiBound=mono.Bounds(:,refcellID,cell2ID,1);
            loBound=mono.Bounds(:,refcellID,cell2ID,2);
            sig = cch < loBound;
            
            prebins = round(length(cch)/2 - .0092/binSize):round(length(cch)/2);
            postbins = round(length(cch)/2 + .0008/binSize):round(length(cch)/2 + .01/binSize);
            cchud  = flipud(cch);
            sigud  = flipud(sig);
            sigpost=max(cch(postbins))>poissinv(1-alpha,max(cch(prebins)));
            sigpre=max(cchud(postbins))>poissinv(1-alpha,max(cchud(prebins)));
            
            if (any(sigud(prebins)) && sigpre)
                sig_con = [sig_con;cell2ID refcellID];
                
            end
            if any(sig(postbins)) && sigpost
                sig_con = [sig_con;refcellID cell2ID];
            end
        end
    end
    sig_con=unique(sig_con,'row');
    ccgR=cell2mat(cellfun(@(x)mono.ccgR(:,x(1),x(2)),mat2cell(sig_con,ones(size(sig_con,1),1),2),'UniformOutput',false)')';
    inh_conn=arrayfun(@(x)mono.completeIndex(mono.completeIndex(:,3)==x,2),sig_con);
    save(fullfile(homedir,'xcorr','bz0902',replace(file(i,1).name,'duo','inh')),'inh_conn','sig_con','ccgR')
end