%% SVM
clearvars -Except sums conn
clc
homedir='F:/pixel-dualtask';
if ~(exist('sums','var')||exist('conn','var'))
    f=dir(fullfile(homedir,'xcorr','coding','fc_decoding_f*.mat'));
    sums=[];
    for i=1:size(f,1)
        fstr=load(fullfile(f(i).folder,f(i).name));
        sums=cat(1,sums,fstr.sums0);
    end
    sums(:,4)=mat2cell(util.mem_type(homedir,cell2mat(sums(:,2))),ones(1,size(sums,1)),2);
end


%%
mem=["congru","incongru","non_mem","mem_nonmem","nonmem_mem"];
trial_type=["distractorNo-correct","distractorNoGo-correct","distractorGo-correct"...
    ,"distractorNo-error","distractorGo-error","distractorNoGo-error"];
NTRIAL=10;
NERR=5;
NRPT=100;
if ~exist('conn','var')
    conn=cell(0);
    sums_mem=cell2mat(sums(:,4));
    for midx=mem
        switch(midx)
            case 'congru'
                conn{end+1}=sums(all(ismember(sums_mem,1:3),2)|all(ismember(sums_mem,4:6),2),[3,2]);
            case 'incongru'
                conn{end+1}=sums(any(ismember(sums_mem,1:3),2)&any(ismember(sums_mem,4:6),2),[3,2]);
            case 'non_mem'
                conn{end+1}=sums(all(ismember(sums_mem,-1),2),[3,2]);
            case 'mem_nonmem'                
                conn{end+1}=sums(any(ismember(sums_mem(:,1),1:6),2)&any(ismember(sums_mem(:,2),-1),2),[3,2]);
            case 'nonmem_mem'
                conn{end+1}=sums(any(ismember(sums_mem(:,1),-1),2)&any(ismember(sums_mem(:,2),1:6),2),[3,2]);
        end
    end
% clear sums
end
%%    

for tridx=1:3
    for midx=1:5       
        clear centt stdd per_trial_norm
        for i=1:size(conn{midx},1)
            % normal-average trial
            centt(i,:)=mean(cell2mat(cellfun(@(x)squeeze(x(1,:,12)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)),2);
            stdd(i,:)=std(cell2mat(cellfun(@(x)squeeze(x(1,:,12)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)),0,2);
            per_trial_norm{i,1}=(squeeze(conn{midx}{i,1}{2*tridx-1}(1,:,12))-centt(i,:))'./stdd(i,:);
            per_trial_norm{i,2}=(squeeze(conn{midx}{i,1}{2*tridx}(1,:,12))-centt(i,:))'./stdd(i,:); 
            
            centt1(i,:)=mean(cell2mat(cellfun(@(x)squeeze(x(3,:,12)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)),2);       
            stdd1(i,:)=std(cell2mat(cellfun(@(x)squeeze(x(3,:,12)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)),0,2);          
            per_trial_norm1{i,1}=(squeeze(conn{midx}{i,1}{2*tridx-1}(3,:,12))-centt1(i,:))'./stdd1(i,:);
            per_trial_norm1{i,2}=(squeeze(conn{midx}{i,1}{2*tridx}(3,:,12))-centt1(i,:))'./stdd1(i,:); 
            
%             % normal-baseline
%             centt(i,:)=mean(cell2mat(cellfun(@(x)squeeze(x(1,:,1:2)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)'),'all');       
%             stdd(i,:)=std(cell2mat(cellfun(@(x)squeeze(x(1,:,1:2)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)'),0,'all');  %             
%             per_trial_norm{i,1}=(squeeze(conn{midx}{i,1}{2*tridx-1}(1,:,12))-centt(i,:))'./stdd(i,:);
%             per_trial_norm{i,2}=(squeeze(conn{midx}{i,1}{2*tridx}(1,:,12))-centt(i,:))'./stdd(i,:); 
%             
%             centt1(i,:)=mean(cell2mat(cellfun(@(x)squeeze(x(3,:,1:2)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)'),'all');       
%             stdd1(i,:)=std(cell2mat(cellfun(@(x)squeeze(x(3,:,1:2)),conn{midx}{i,1}(2*tridx-1:2*tridx),'UniformOutput',false)'),0,'all');              
%             per_trial_norm1{i,1}=(squeeze(conn{midx}{i,1}{2*tridx-1}(3,:,12))-centt1(i,:))'./stdd1(i,:);
%             per_trial_norm1{i,2}=(squeeze(conn{midx}{i,1}{2*tridx}(3,:,12))-centt1(i,:))'./stdd1(i,:); 
            
        end
        for NSU=1000
            [decdata.s1,decdata.s2]=deal(nan(NSU,NTRIAL,1,NRPT));
            su_sel=find(min(cellfun(@(x) size(x,1),per_trial_norm(:,1:2)),[],2)>NTRIAL & all(centt>0.5&stdd>0,2));
            conn_id=cell2mat(conn{midx}(:,2));
            for rpt=1:NRPT
                su_perm=randsample(su_sel,NSU);
                temp=arrayfun(@(x) per_trial_norm{x,1}(randperm(size(per_trial_norm{x,1},1),NTRIAL),:),su_perm,'UniformOutput',false);
                decdata.s1(:,:,:,rpt)=permute(cat(3,temp{:}),[3,1,2]);
                temp=arrayfun(@(x) per_trial_norm{x,2}(randperm(size(per_trial_norm{x,2},1),NTRIAL),:),su_perm,'UniformOutput',false);
                decdata.s2(:,:,:,rpt)=permute(cat(3,temp{:}),[3,1,2]);   
            end
            
            [decdata.train,decdata.test]=deal(decdata);
            decdata=rmfield(decdata,{'s1','s2'});
            cvcorr.(mem(midx)).(sprintf('%s_SU%d',replace(trial_type(tridx),'-','_'),NSU))=Decoding(decdata,'decoder','SVM');
            
            if true
                % FR
                id=cell2mat(conn{midx}(su_sel,2));
                if numel(unique(id(:,1)))>500
                    su_perm=randsample(unique(id(:,1)),500);
                else
                    su_perm=unique(id(:,1));
                end
                for rpt=1:NRPT
%                     su_perm=randsample(su_sel,NSU);
%                     temp=arrayfun(@(x) per_trial_norm{x,1}(randi(size(per_trial_norm{x,1},1),NTRIAL,1),:),su_perm,'UniformOutput',false);
%                     decdata1.s1(:,:,:,rpt)=permute(cat(3,temp{:}),[3,1,2]);
%                     temp=arrayfun(@(x) per_trial_norm{x,2}(randi(size(per_trial_norm{x,2},1),NTRIAL,1),:),su_perm,'UniformOutput',false);
%                     decdata1.s2(:,:,:,rpt)=permute(cat(3,temp{:}),[3,1,2]);
                    temp=arrayfun(@(x) per_trial_norm1{x,1}(randperm(size(per_trial_norm1{x,1},1),NTRIAL),:)...
                        ,arrayfun(@(y)find(conn_id(:,1)==y,1),su_perm),'UniformOutput',false);
                    decdata1.s1(:,:,:,rpt)=permute(cat(3,temp{:}),[3,1,2]);
                    temp=arrayfun(@(x) per_trial_norm1{x,2}(randperm(size(per_trial_norm1{x,2},1),NTRIAL),:)...
                        ,arrayfun(@(y)find(conn_id(:,1)==y,1),su_perm),'UniformOutput',false);
                    decdata1.s2(:,:,:,rpt)=permute(cat(3,temp{:}),[3,1,2]);
                end
                %             decdata=NeuronSetForDecoding(homedir,su_perm,'trial',trial_type(tridx),'rps',NRPT,'trialnum',NTRIAL);
                [decdata1.train,decdata1.test]=deal(decdata1);
                decdata1=rmfield(decdata1,{'s1','s2'});
                ctrlcorr.(mem(midx)).(sprintf('%s_SU%d',replace(trial_type(tridx),'-','_'),NSU))=Decoding(decdata1,'decoder','SVM');
           end
        end
    end
end
%%
% close
fh=figure('Color','w','Position',[100,100,500,250]);
c={[162,73,0.01]/255,[97,139,0.01]/255,[0.01,0.01,0.01]/255,[0.5,0.5,0.5],[0.5,0.75,0.75]};

for memidx=1:5
    for tridx=1:3
        for i=1:length(NSU)
            dec(memidx,tridx,i,:)=cell2mat(cvcorr.(sprintf('%s',mem(memidx)))...
                .(sprintf('%s_SU%d',replace(trial_type(tridx),'-','_'),NSU(i))).cvcorr);
            dec_shuffle(memidx,tridx,i,:)=cell2mat(cvcorr.(sprintf('%s',mem(memidx)))...
                .(sprintf('%s_SU%d',replace(trial_type(tridx),'-','_'),NSU(i))).shufcorr);
            
        end
    end
end
subplot(1,2,1)
hold on
bh=bar(mean(dec(:,1:3,:,:),4)',0.8);
for memidx=1:5
    bh(memidx).EdgeColor=c{memidx};
    bh(memidx).FaceColor='w';
    errorbar(bh(memidx).XEndPoints,bh(memidx).YEndPoints,std(dec(memidx,1:3,:,:),0,4)'/sqrt(NRPT),'LineStyle','none','Color',c{memidx},'LineWidth',1)
end
ChanceLevel=mean(dec_shuffle,'all');
plot([0,4],[ChanceLevel,ChanceLevel],'r--')
% legend(bh,{'Same memory','Diff. memory','Non-memory'})
set(gca(),'XTick',1:3,'XTickLabel',{'Distractor No','Distractor No-go','Distractor Go'},'XTickLabelRotation',30,'YLim',[0,1])
ylabel('Decoding accuracy');
title('FCSP')


for memidx=1:5
    for tridx=1:3
        for i=1:length(NSU)
            dec1(memidx,tridx,i,:)=cell2mat(ctrlcorr.(sprintf('%s',mem(memidx)))...
                .(sprintf('%s_SU%d',replace(trial_type(tridx),'-','_'),NSU(i))).cvcorr);
            dec1_shuffle(memidx,tridx,i,:)=cell2mat(ctrlcorr.(sprintf('%s',mem(memidx)))...
                .(sprintf('%s_SU%d',replace(trial_type(tridx),'-','_'),NSU(i))).shufcorr);
        end
    end
end
subplot(1,2,2)
hold on
bh=bar(mean(dec1(:,1:3,:,:),4)',0.8);
for memidx=1:5
    bh(memidx).EdgeColor=c{memidx};
    bh(memidx).FaceColor='w';
    errorbar(bh(memidx).XEndPoints,bh(memidx).YEndPoints,std(dec1(memidx,1:3,:,:),0,4)'/sqrt(NRPT),'LineStyle','none','Color',c{memidx},'LineWidth',1)
end
ChanceLevel=mean(dec1_shuffle(memidx,1:3,:,:),'all');
plot([0,4],[ChanceLevel,ChanceLevel],'r--')
legend(bh,{'Same memory','Diff. memory','Non-memory'})
set(gca(),'XTick',1:3,'XTickLabel',{'Distractor No','Distractor No-go','Distractor Go'},'XTickLabelRotation',30,'YLim',[0,1])
ylabel('Decoding accuracy');
title('FR')
exportgraphics(fh,fullfile(homedir,'xcorr','FC_decoding_inh.pdf'))
return
%%
homedir='F:/pixel-dualtask';
load(fullfile(homedir,'xcorr','coding','fc_SVM.mat'),'cvcorr')

NSU=500;
trial_type=["distractorNo_correct","distractorGo_correct","distractorNoGo_correct"];
mem=["congru","incongru","non_mem"];

for memidx=1:3
for tridx=1:3
    for i=1:length(NSU)
        dec_d1(memidx,tridx,i,:)=cell2mat(cvcorr.(sprintf('%s',mem(memidx))).(sprintf('%s_SU%d',trial_type(tridx),NSU(i))).cvcorr(5:6));
        dec_d2(memidx,tridx,i,:)=cell2mat(cvcorr.(sprintf('%s',mem(memidx))).(sprintf('%s_SU%d',trial_type(tridx),NSU(i))).cvcorr(8:9));
        dec_d3(memidx,tridx,i,:)=cell2mat(cvcorr.(sprintf('%s',mem(memidx))).(sprintf('%s_SU%d',trial_type(tridx),NSU(i))).cvcorr(12));
    end
end
end
%%
close
fh=figure('Color','w','Position',[100,100,500,250]);
c={'r','b','k'};
for tridx=1:3
    subplot(2,3,tridx)
    hold on
    X=[];g1=[];g2=[];
    for memidx=1:3        
        plotLine(squeeze(dec_d2(memidx,tridx,:,:)),c{memidx},'-')
        title(sprintf('%s',trial_type(tridx)))
%         set(gca,'XTick',1:6,'XTickLabel',50:50:300)
        ylabel('WM classification accuracy');
        xlabel('Number of FC pairs');
%         xlim([1,6])
        ylim([0.6,1])
        X=[X,squeeze(dec_d2(memidx,tridx,:,:))];
        g1=[g1,memidx*ones(6,100)];
        g2=[g2,kron([1:6]',ones(1,100))];
    end
    %anova
%     p{tridx,1}=anovan(reshape(X,1,[]),{reshape(g1,1,[]),reshape(g2,1,[])});
    
    subplot(2,3,tridx+3)
    hold on
    X=[];g1=[];g2=[];
    for memidx=1:3        
        plotLine(squeeze(dec_d3(memidx,tridx,:,:)),c{memidx},'-')
        set(gca,'XTick',1:6,'XTickLabel',50:50:300)
        ylabel('WM classification accuracy');
        xlabel('Number of FC pairs');
        xlim([1,6])
        ylim([0.6,1])
        X=[X,squeeze(dec_d3(memidx,tridx,:,:))];
        g1=[g1,memidx*ones(6,50)];
        g2=[g2,kron([1:6]',ones(1,50))];
    end
    %anova
%     p{tridx,2}=anovan(reshape(X,1,[]),{reshape(g1,1,[]),reshape(g2,1,[])});
end
exportgraphics(fh,fullfile(homedir,'xcorr','FC_decoding.pdf'))



%%

if size(cell2mat(cvcorr{i}.cvcorr')',2)==56
%     arrayfun(@(x)plot(repmat(x,1,2),[0.3,1],'k--'),[12.5,16.5,24.5,28.5,48.5,52.5])
    cellfun(@(x)fill([x,flip(x)],[0.3,0.3,1,1],[0.8500 0.3250 0.0980],'EdgeColor','none'),{[12.5,16.5],[48.5,52.5]})
    cellfun(@(x)fill([x,flip(x)],[0.3,0.3,1,1],[0.4660 0.6740 0.1880],'EdgeColor','none'),{[24.5,28.5],[36.5,40.5]})
    alpha(0.1)
elseif size(cell2mat(cvcorr{i}.cvcorr')',2)==140
%     arrayfun(@(x)plot(repmat(x,1,2),[0.3,1],'k--'),[29.5,40.5,59.5,70.5,89.5,100.5,109.5,120.5])
    cellfun(@(x)fill([x,flip(x)],[0.3,0.3,1,1],[0.8500 0.553333250 0.0980],'EdgeColor','none'),{[29.5,40.5],[109.5,120.5]})
    cellfun(@(x)fill([x,flip(x)],[0.3,0.3,1,1],'EdgeColor','none'),{[59.5,70.5],[89.5,100.5]})
    alpha(0.1)
end
saveas(gcf,fullfile(HomePath,'decoding',sprintf('%s_decoding',prefix)),'png')


colors={'r','b','k'};
fh=figure('Color','w','Position',[100,100,235,235]);
hold on;
for mi=1:3
    cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('SU%d',x)){mi}.cvcorr{1},50:50:300,'UniformOutput',false));
    mm=mean(cvmat).*100;
    cim=bootci(500,@(x) mean(x).*100,cvmat);
    ph(mi)=plot(50:50:300,mm,strjoin({'-',colors{mi}},''));
    fill([50:50:300,fliplr(50:50:300)],[cim(1,:),fliplr(cim(2,:))],colors{mi},'EdgeColor','none','FaceAlpha',0.1);
end

cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('SU%d',x)){mi}.errcorr{1},50:50:300,'UniformOutput',false));
mm=mean(cvmat).*100;
cim=bootci(500,@(x) mean(x).*100,cvmat);
ph(4)=plot(50:50:300,mm,'--r');
fill([50:50:300,fliplr(50:50:300)],[cim(1,:),fliplr(cim(2,:))],'r','EdgeColor','none','FaceAlpha',0.1);

ylabel('WM classification accuracy');
xlabel('Number of FC pairs');
xlim([50,300])
legend(ph,{'Same memory','Diff. memory','Nonmemory','(Error trial)'},Location,"northoutside");

exportgraphics(fh,'FC_coding_SVM.pdf')

%%
function [out,NeuronSet]=NeuronSetForDecoding(HomePath,id,opt)
arguments  
    HomePath (1,:) char
    id (:,1) double
    opt.trial (1,:) char 
    opt.rps (1,:) double = 1
    opt.NeuronNumConstrcit (1,:) logical = true
    opt.trialnum (1,1) double
    opt.neuronnum (1,1) double = 10000 
    opt.FRFile(1,:) char ='FR_All_1000ms.hdf5'
end

cluster_id=h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/cluster_id');
path=regexp(h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/path'),'(\w|\\|-)*','match','once');
NeuronSet.id=id;
NeuronSet.path=arrayfun(@(x)path{ismember(cluster_id,x),1},NeuronSet.id,'UniformOutput',false);

out.s1=[];
out.s2=[];
for idx=1:length(NeuronSet.id)    
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
   
    temp.FR_All=h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},opt.FRFile),'/FR_All');    
    temp.trials=util.markLPerf(h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},opt.FRFile),'/Trials'),0.7,0.8,120);
    temp.SU_id=h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},opt.FRFile),'/SU_id');    
   
    [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',opt.trial); 
    if isempty(sel_S1)||isempty(sel_S2)||nnz(sel_S1)<5||nnz(sel_S2)<5
        idx=idx+ncell;
        clear sel_S1 sel_S2 FR temp s1 s2
        continue
    end
    temp.FR=temp.FR_All(:,:,temp.SU_id==rem(NeuronSet.id(idx),10000));
    %     FR=(temp.FR-mean(temp.FR,1))./repmat(std(temp.FR),size(temp.FR,1),1,1);
 
    temp.FR_base=reshape(temp.FR([sel_S1;sel_S2],1:2,:),[],1);        
  
    
    FR=(temp.FR-mean(temp.FR_base,1))./std(temp.FR_base,0,1);

%     FR=(temp.FR-mean(temp.FR([sel_S1;sel_S2],:,:),1))./std(temp.FR([sel_S1;sel_S2],:,:),0,1);
%     
%     
    temp.s1=reshape(arrayfun(@(z)FR(z,:,:),cell2mat(arrayfun(@(y)sel_S1(y),cell2mat(arrayfun(@(x)randperm(length(sel_S1)...
        ,opt.trialnum),1:opt.rps,'UniformOutput',false)')','UniformOutput',false)),'UniformOutput',false),[],1);
    temp.s2=reshape(arrayfun(@(z)FR(z,:,:),cell2mat(arrayfun(@(y)sel_S2(y),cell2mat(arrayfun(@(x)randperm(length(sel_S2)...
        ,opt.trialnum),1:opt.rps,'UniformOutput',false)')','UniformOutput',false)),'UniformOutput',false),[],1);
    s1=permute(reshape(cat(1,temp.s1{:}),opt.trialnum,opt.rps,size(FR,2),[]),[4,1,3,2]);
    s2=permute(reshape(cat(1,temp.s2{:}),opt.trialnum,opt.rps,size(FR,2),[]),[4,1,3,2]);
    out.s1=cat(1,out.s1,s1(:,:,12,:));
    out.s2=cat(1,out.s2,s2(:,:,12,:));
    
    
    clear sel_S1 sel_S2 FR temp s1 s2
end

end
function [out,NeuronSet]=NeuronSetForDecoding1(HomePath,id,opt)
arguments  
    HomePath (1,:) char
    id (:,1) double
    opt.trial (1,:) char 
    opt.rps (1,:) double = 1
    opt.NeuronNumConstrcit (1,:) logical = true
    opt.trialnum (1,1) double
    opt.neuronnum (1,1) double = 10000 
    opt.FRFile(1,:) char ='FR_All_1000ms.hdf5'
end

cluster_id=h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/cluster_id');
path=regexp(h5read(fullfile(HomePath,'Selectivity_1129.hdf5'),'/path'),'(\w|\\|-)*','match','once');
NeuronSet.id=id;
NeuronSet.path=arrayfun(@(x)path{ismember(cluster_id,x),1},NeuronSet.id,'UniformOutput',false);

for idx=1:length(NeuronSet.id)    
    if ~ispc
        NeuronSet.path{idx}=replace(NeuronSet.path{idx},'\','/');
    end
   
    temp.FR_All=h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},opt.FRFile),'/FR_All');    
    temp.trials=util.markLPerf(h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},opt.FRFile),'/Trials'),0.7,0.8,120);
    temp.SU_id=h5read(fullfile(HomePath,'DataSum',NeuronSet.path{idx},opt.FRFile),'/SU_id');    
   
    [sel_S1,sel_S2]=util.ExtractTrial(temp.trials,'task','dualtask','trials',opt.trial); 
    
    temp.FR=temp.FR_All(:,:,temp.SU_id==rem(NeuronSet.id(idx),10000)); 
    temp.FR_base=reshape(temp.FR([sel_S1;sel_S2],1:2,:),[],1);    
    FR=(temp.FR-mean(temp.FR_base,1))./std(temp.FR_base,0,1);
    
    out{1}=FR(sel_S1,12);
    out{2}=FR(sel_S2,12);
end
end

function plotLine(Data,c,l)
m=mean(Data,1);
plot(1:size(m,2),m,c,'LineStyle',l);
ci=bootci(1000,@(x) mean(x),Data);
fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.2);
end