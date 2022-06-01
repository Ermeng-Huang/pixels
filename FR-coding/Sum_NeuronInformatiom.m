%%% Follow code 'Selectivity.m' and 'regionAlignment.m'
clear
clc
homedir='F:\pixel-dualtask';
%%
load(fullfile(homedir,'sust_trans_1_206.mat'),'sust_trans','path')
sel=double(FindSel1(sust_trans)); % old criterion

fl=unique(path);
SU_id=[];reg=[];
for i=1:length(fl)
    [SU_id_temp,reg_temp]=plotOneTrack(fullfile(homedir,'DataSum',fl{i}));
    SU_id=[SU_id;SU_id_temp];
    reg=[reg;reg_temp];
end

if ~all(SU_id==rem(sust_trans(:,end),10000))
    disp('check!')
end 
hdf5write(fullfile(homedir,'Selectivity_0327.hdf5'), '/Pvalue_noPermutaion', sust_trans(:,1:end-1),'/cluster_id', sust_trans(:,end)...
    ,'/path', path,'/reg', reg,'/sust_trans_noPermutaion',sel);

%% Main function
function sel=FindSel(sust_trans)
% for two-way anova
sel_temp=abs(sust_trans(:,3:end-1))<0.05/(size(sust_trans(:,3:end-1),2)/2);
%find sustained and transient
sel(:,1)=all(sel_temp(:,1:end/2),2)+2*all(sel_temp(:,end/2+1:end),2);
sel(:,2)=any(sel_temp(:,1:end/2),2)+2*any(sel_temp(:,end/2+1:end),2);


pre1=arrayfun(@(x,y)x>0&y==1,sust_trans(:,3:end-1),sel_temp);
pre2=arrayfun(@(x,y)x<0&y==1,sust_trans(:,3:end-1),sel_temp)*2;
pre=pre1+pre2;
sel(any(pre==1,2)&any(pre==2,2),1:2)=0;

% sel(:,1)=(all(sel_temp(:,1:end/2),2)&(all(pre(:,1:end/2)==1,2)|all(pre(:,1:end/2)==2,2)))...
%     +2*(all(sel_temp(:,end/2+1:end),2)&(all(pre(:,end/2+1:end)==1,2)|all(pre(:,end/2+1:end)==2,2)));
% sel(:,2)=(any(sel_temp(:,1:end/2),2)& ((any(pre(:,1:end/2)==1,2)&~any(pre(:,1:end/2)==2,2))|(~any(pre(:,1:end/2)==1,2)&any(pre(:,1:end/2)==2,2))) )...
%     +2*(any(sel_temp(:,end/2+1:end),2)&((any(pre(:,end/2+1:end)==1,2)&~any(pre(:,end/2+1:end)==2,2))|(~any(pre(:,end/2+1:end)==1,2)&any(pre(:,end/2+1:end)==2,2))));
sel=double(sel);
sel(:,5:size(sust_trans(:,3:end-1),2)+4)=pre;

sel(:,3:4)=1*(sust_trans(:,1:2)<0.05&sust_trans(:,1:2)>0);
sel(:,3:4)=2*(sust_trans(:,1:2)>-0.05&sust_trans(:,1:2)<0);
end

function sel=FindSel1(sust_trans)
% for ranksum
sel_temp=abs(sust_trans(:,3:10))<(0.01);
%find sustained and transient
sel(:,1)=all(sel_temp,2);
sel(:,2)=any(sel_temp,2)&~all(sel_temp,2);
pre1=arrayfun(@(x,y)x>0&y==1,sust_trans(:,3:10),sel_temp);
pre2=arrayfun(@(x,y)x<0&y==1,sust_trans(:,3:10),sel_temp)*2;
pre=pre1+pre2;
sel(any(pre==1,2)&any(pre==2,2),1:2)=0;

sel(:,5:12)=sel_temp;
sel(:,3:4)=1*(sust_trans(:,1:2)<0.01&sust_trans(:,1:2)>0);
sel(:,3:4)=2*(sust_trans(:,1:2)>-0.01&sust_trans(:,1:2)<0);

end

function [SU_id,reg]=plotOneTrack(rootpath)
SU_id=plotOneDir(rootpath);
try
    reg=table2cell(readtable(fullfile(rootpath,'su_id2reg_allclass.csv'),'Format','%s%s%s%s%s%s%s%s%s%s%s')); %%%
catch
%     reg=cellstr(cell(length(SU_id),11));
    reg=[];
end
end

function [cluster_ids]=plotOneDir(rootpath)
sps=30000;
FR_Th=1.0;
metaf=ls(fullfile(rootpath,'*.meta'));
fh=fopen(fullfile(rootpath,metaf));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/sps/2*FR_Th;
clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
waveformGood=strcmp(clusterInfo{:,4},'good');
freqGood=clusterInfo{:,10}>spkNThresh;
cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
end