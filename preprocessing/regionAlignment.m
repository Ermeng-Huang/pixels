%??????sorting??????????
clear
clc
%%
if false
track_table=readtable('D:\pixel-optogenetic\Neuropixels_AIopto_track_20210307.xlsx','ReadVariableNames',true);
cd('D:\pixel-optogenetic')
fl=dir('D:\pixel-optogenetic\DataSum\*\*\M*meta');
su_reg=[];

for onefile=fl'    
    temp=reg(onefile,track_table);
    su_reg=[su_reg;temp];    
end
save('reg_path.mat','su_reg')
end
%% calcualte the count of single unit in per region

load('D:\pixel-optogenetic\reg_path.mat')
load('D:\pixel-optogenetic\Performance_allneuron_0308.mat')
reglist=unique(su_reg(:,2));
reg_none={'fiber tracts','VS','unlabeled','root'};
reglist(ismember(reglist,reg_none))=[];
for i=1:size(reglist,1)
   reglist{i,2}=nnz(strcmp(su_reg(:,2),reglist{i,1}));  
   reglist{i,3}=nnz(strcmp(su_reg(PerfList(:,2)==1,2),reglist{i,1}));  %learning
   reglist{i,4}=nnz(strcmp(su_reg(PerfList(:,3)==1,2),reglist{i,1})); 
   reg_logic(i,1)=reglist{i,3}>40 ;
   reg_logic(i,2)=reglist{i,4}>40 ;
   reg_logic(i,3)=reglist{i,2}>40 ;
   
end
writetable(cell2table(reglist),'reg_num_20210412.csv')
reg_L=reglist(reg_logic(:,1),1);
reg_W=reglist(reg_logic(:,2),1);
reg_keep=reglist(reg_logic(:,3),1);
save('reg_keep_40_0412.mat','reg_L','reg_W')
for i=1:length(reg_L)
    reg_L{i,2}=reg2ccfid(reg_L{i});
end
writetable(cell2table(reg_L),'reg_L_20210412.csv')

%% Function

function region=reg(onefile,track_table) 
n=0;
rootpath=onefile.folder;
filename=split(onefile.name,'_');
temp=track_table(strcmp(track_table.MouseID,filename{1,1})&track_table.Date==str2double(filename{2,1})&track_table.ImecID==str2double(filename{5,1}(1)),:);
probeID=temp.ProbeID;
h=readtable(['F:\pixel-optogenetic\Histology\' filename{1,1} '\results\Probe ' num2str(probeID) '.csv']);

load(['F:\pixel-optogenetic\Histology\' filename{1,1} '\results\HistLength.mat']);
load('D:\pixel-optogenetic\reg_ccfid_map','reg2tree')
ProbeLength=Probe_length(probeID);
% h.tip2upperBorder=max(h.lowerBorder)-h.upperBorder;
% h.tip2lowerBorder=max(h.lowerBorder)-h.lowerBorder;

clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
metaf=ls(fullfile(rootpath,'*.meta'));
fh=fopen(fullfile(rootpath,metaf));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/30000/2*1.0;
clusterInfo_good=clusterInfo(strcmp(clusterInfo.KSLabel,'good') & clusterInfo.n_spikes>spkNThresh,:);
disp(rootpath)    
for i=1:size(clusterInfo_good,1)
    depth=ProbeLength-clusterInfo_good.depth(i);
    su_id=clusterInfo_good.id(i);            
    try
        if isempty(cell2mat(h{h.upperBorder<=depth& h.lowerBorder>depth,3}))            
            region_su='unlabeled';
            n=n+1;
            disp(n)       
        else
            reg_class=reg2tree(cell2mat(h{h.upperBorder<=depth& h.lowerBorder>depth,3}));
            if length(reg_class)==1
                region_su=reg_class{1};
            elseif strcmp(reg_class{2},'grey')
                try
                    region_su=reg_class{7};
                catch
                    region_su=reg_class{end};
                end
            elseif strcmp(reg_class{2},'fiber tracts')
                region_su=reg_class{2};
            elseif strcmp(reg_class{2},'VS')
                region_su=reg_class{2};
            else
                region_su=reg_class{end};
            end
            clear reg_class
        end
    catch
        disp(reg_class{1})
    end
    region{i,1}=su_id;
    region{i,2}=region_su;
    region{i,3}=regexp(rootpath,'(?<=DataSum\\)(.*)','match','once');

end
su_id2reg=table(region);
writetable(su_id2reg,fullfile(rootpath,'su_id2reg_5.csv'))

end

% sus_trans=h5read('D:\code\transient_6_AIopto_learning.hdf5','/sus_trans');
% reg_list=h5read('D:\code\transient_6_AIopto_learning.hdf5','/reg');
% cid_list=h5read('D:\code\transient_6_AIopto_learning.hdf5','/cluster_id');
% path_list=h5read('D:\code\transient_6_AIopto_learning.hdf5','/path');
% 
% for i=1:size(reg_list)    
%     if ~isempty(regexp(reg_list{i,1},'CA[123]', 'once'))
%         rename=regexp(reg_list{i,1},'CA[123]','match','once');
%     elseif ~isempty(regexp(reg_list{i,1},'DG', 'once'))
%         rename=regexp(reg_list{i,1},'DG','match','once');            
%     else
%         rename=regexp(reg_list{i,1},'[A-Za-z-]{1,}','match','once');
%     end
%     if ~isempty(rename)
%         reg_list{i,1}=rename;
%     end
% end
% reg=unique(table2cell(su_reg(:,2)));
% for i=1:size(reg,1)
%    reg{i,2}=nnz(strcmp(reg_list,reg{i,1})); 
% end
% writetable(cell2table(reg),'reg_num.csv')
