clear
clc
homedir='F:/pixel-dualtask';
reg=["within_region","cross_region"];
mem=["congru","incongru","non_mem","mem_nonmem","nonmem_mem"];
trial_type=["distractorNo-correct","distractorNoGo-correct","distractorGo-correct"...
        ,"distractorNo-error","distractorNoGo-error","distractorGo-error"];
trial_color=["k","b","r","k","b","r"];
trial_linestyle=["-","-","-","--","--","--"];
phase=["Early delay","Middle delay","Late delay"];
%% extract data
f=dir(fullfile(homedir,'xcorr','coding','fc_coding_f*.mat'));
[sums,conn]=deal([]);
for i=1:size(f,1)
    fstr=load(fullfile(f(i).folder,f(i).name));
    sums=cat(1,sums,fstr.sums1);
end
sums(:,5)=mat2cell(util.mem_type(homedir,cell2mat(sums(:,3))),ones(1,size(sums,1)),2);
sums(:,6)=mat2cell(util.IDtoReg(homedir,cell2mat(sums(:,3))),ones(1,size(sums,1)),2);
sums_mem=cell2mat(sums(:,5));
for regidx=reg
    switch(regidx)         
        case 'within_region'
            conn_reg=cellfun(@(x)strcmp(x{1}{7},x{2}{7}),sums(:,6))&cellfun(@(x)strcmp(x{1}{4},'CTX')&strcmp(x{2}{4},'CTX'),sums(:,6))...
                &~cellfun(@(x)isempty(x{1}{7})|isempty(x{2}{7}),sums(:,6));
        case 'cross_region'
            conn_reg=~cellfun(@(x)strcmp(x{1}{7},x{2}{7}),sums(:,6))&cellfun(@(x)strcmp(x{1}{4},'CTX')&strcmp(x{2}{4},'CTX'),sums(:,6))...
                &~cellfun(@(x)isempty(x{1}{7})|isempty(x{2}{7}),sums(:,6));
    end
    for midx=mem        
        switch(midx)
            case 'congru'
                conn{end+1}=sums((all(ismember(sums_mem,1:3),2)|all(ismember(sums_mem,4:6),2))&conn_reg,:);
            case 'incongru'
                conn{end+1}=sums((any(ismember(sums_mem,1:3),2)&any(ismember(sums_mem,4:6),2))&conn_reg,:);
            case 'non_mem'
                conn{end+1}=sums((all(ismember(sums_mem,-1),2))&conn_reg,:);
            case 'mem_nonmem'
                conn{end+1}=sums((any(ismember(sums_mem(:,1),1:6),2)&any(ismember(sums_mem(:,2),-1),2))&conn_reg,:);
            case 'nonmem_mem'
                conn{end+1}=sums((any(ismember(sums_mem(:,1),-1),2)&any(ismember(sums_mem(:,2),1:6),2))&conn_reg,:);
        end
    end
end
%% FC selectivity index (only for congruent pair)
f=figure('Color','w','Position',[100,100,400,250]);
for regidx=1:2
    FCSP_sel=[];
    fc_congru=conn{5*(regidx-1)+1};
    s=2*all(ismember(cell2mat(fc_congru(:,5)),1:3),2)-1;
    for t=1:length(trial_type)
        FCSP_sel{t,1}=s.*cellfun(@(x)mean((x(5:6,1,t)-x(5:6,4,t))./(x(5:6,1,t)+x(5:6,4,t))),fc_congru(:,4));
        FCSP_sel{t,2}=s.*cellfun(@(x)mean((x(8:9,1,t)-x(8:9,4,t))./(x(8:9,1,t)+x(8:9,4,t))),fc_congru(:,4));
        FCSP_sel{t,3}=s.*cellfun(@(x)(x(12,1,t)-x(12,4,t))./(x(12,1,t)+x(12,4,t)),fc_congru(:,4));
    end
    avail=all(cell2mat(reshape(cellfun(@(x)~isnan(x),FCSP_sel,'UniformOutput',false),1,[])),2);
    FCSP_sel=cellfun(@(x)x(avail,1),FCSP_sel,'UniformOutput',false);
    FCSP_mean=cellfun(@(x)mean(x),FCSP_sel);
    FCSP_std=cellfun(@(x)std(x)/sqrt(length(x)),FCSP_sel);    
    
    for i=1:3
        subplot(2,3,i+3*(regidx-1))
        hold on
        for t=1:length(trial_type)
            bar(t,FCSP_mean(t,i),'EdgeColor',trial_color(t),'FaceColor','none','LineStyle',trial_linestyle(t))
            errorbar(t,FCSP_mean(t,i),FCSP_std(t,i),'Color',trial_color(t),'LineStyle','none')
        end
        set(gca,'XTick',[],'ylim',[0 0.25])
        xlabel(phase(i))
        ylabel('FCSP sel. indx')
        g1=cell2mat(arrayfun(@(x)x*ones(1,nnz(avail)),[1,1,1,2,2,2],'UniformOutput',false));
        g2=cell2mat(arrayfun(@(x)x*ones(1,nnz(avail)),[1,2,3,1,2,3],'UniformOutput',false));
        g3=cell2mat(arrayfun(@(x)x*ones(1,nnz(avail)),[1,2,3],'UniformOutput',false));
        [p1(:,i),~,~]=anovan(cell2mat(FCSP_sel(:,i)),{g1,g2},'varnames',{'Performance','distractor'},'display','off');
        [p2(:,i),~,~]=anovan(cell2mat(FCSP_sel(1:end/2,i)),{g3},'varnames',{'distractor'},'display','off');
        [p2(:,i+3),~,~]=anovan(cell2mat(FCSP_sel(end/2+1:end,i)),{g3},'varnames',{'distractor'},'display','off');
    end
end
exportgraphics(f,fullfile(homedir,'xcorr','FCSP_SelIndex.pdf'))

%% FCSP rate

for regidx=1:numel(reg)
    f=figure('Color','w','Position',[100,100,400,600]);
    for midx=1:numel(mem)
        FCSP_sel=[];
        for t=1:length(trial_type)
            FCSP_sel{t,1}=cellfun(@(x)mean([x(5:6,1,t);x(5:6,4,t)]),conn{5*(regidx-1)+midx}(:,4));
            FCSP_sel{t,2}=cellfun(@(x)mean([x(8:9,1,t);x(8:9,4,t)]),conn{5*(regidx-1)+midx}(:,4));
            FCSP_sel{t,3}=cellfun(@(x)mean([x(12,1,t);x(12,4,t)]),conn{5*(regidx-1)+midx}(:,4));
        end
        avail=all(cell2mat(reshape(cellfun(@(x)~isnan(x),FCSP_sel,'UniformOutput',false),1,[])),2);
        FCSP_sel=cellfun(@(x)x(avail,1),FCSP_sel,'UniformOutput',false);
        FCSP_mean=cellfun(@(x)mean(x),FCSP_sel);
        FCSP_std=cellfun(@(x)std(x)/sqrt(length(x)),FCSP_sel);
        
       
        for i=1:3
            subplot(5,3,i+3*(midx-1))
            hold on
            for t=1:length(trial_type)
                bar(t,FCSP_mean(t,i),'EdgeColor',trial_color(t),'FaceColor','none','LineStyle',trial_linestyle(t))
                errorbar(t,FCSP_mean(t,i),FCSP_std(t,i),'Color',trial_color(t),'LineStyle','none')
            end
            set(gca,'Xtick',[],'ylim',[0 3])
            xlabel(phase(i))
            ylabel('FCSP rate (HZ)')
            g1=cell2mat(arrayfun(@(x)x*ones(1,nnz(avail)),[1,1,1,2,2,2],'UniformOutput',false));
            g2=cell2mat(arrayfun(@(x)x*ones(1,nnz(avail)),[1,2,3,1,2,3],'UniformOutput',false));
            g3=cell2mat(arrayfun(@(x)x*ones(1,nnz(avail)),[1,2,3],'UniformOutput',false));
            [p1(:,i),~,~]=anovan(cell2mat(FCSP_sel(:,i)),{g1,g2},'varnames',{'Performance','distractor'},'display','off');
            [p2(:,i),~,~]=anovan(cell2mat(FCSP_sel(1:end/2,i)),{g3},'varnames',{'distractor'},'display','off');
            [p2(:,i+3),~,~]=anovan(cell2mat(FCSP_sel(end/2+1:end,i)),{g3},'varnames',{'distractor'},'display','off');
        end
    end
    exportgraphics(f,fullfile(homedir,'xcorr',sprintf('FCSP_rate_%s.pdf',reg(regidx))))
end



%%
for t=1:length(trial_type)
    fc_sel{t,1}=sumssssss(cellfun(@(x)any(x(5:6,19,t)<0.05/8),sumssssss(:,4)),:);
    fc_sel{t,2}=sumssssss(cellfun(@(x)any(x(8:9,19,t)<0.05/8),sumssssss(:,4)),:);
    fc_sel{t,3}=sumssssss(cellfun(@(x)any(x(12,19,t)<0.05/8),sumssssss(:,4)),:);
end
sel_frac=cellfun(@(x)length(x),fc_sel)/length(sumssssss);


%%
fc_type=["congru","incongru","non-mem"];
c={'k','r','b','k','r','b'};
f=figure('Color','w','Position',[100,100,400,250]);
hold on
for i=1:3
    subplot(1,3,i)
    hold on
    for t=1:length(trial_type)
        bar(t,sel_frac(t,i),0.18,'k','EdgeColor',c{t})
    end
    p(i,1)=statistics.chi2test(length(sumssssss)*[sel_frac(:,i),1-sel_frac(:,i)]);
    p(i,2)=statistics.chi2test(length(sumssssss)*[sel_frac(1:3,i),1-sel_frac(1:3,i)]);
    p(i,3)=statistics.chi2test(length(sumssssss)*[sel_frac(4:6,i),1-sel_frac(4:6,i)]);
end

%% function
function normalize(FR)

end