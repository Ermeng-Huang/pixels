clc
clear
close all

if ispc
    CodePath='D:\code';
    HomePath='D:\pixel-optogenetic';
else
    CodePath='/home/hem/Code';
    HomePath='/home/hem/datashare/AI-opto';
end
load(fullfile(HomePath,'session_list.mat'))
path=util.Path_default;

sus_trans=h5read(path.selectivity,'/sus_trans_noPermutaion')';
sel=h5read('D:\pixel-optogenetic\selectivity.hdf5','/selectivity');
sust=sel(sus_trans(:,1)==1,:);
tran_1=sel(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==1,4:10);
tran_2=sel(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==2,4:10);
tran_3=sel(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==3,4:10);
tran_4=sel(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==4,4:10);
tran_5=sel(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==5,4:10);

tran_1(sus_trans(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==1,6:12)==0)=0;
tran_2(sus_trans(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==2,6:12)==0)=0;
tran_3(sus_trans(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==3,6:12)==0)=0;
tran_4(sus_trans(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==4,6:12)==0)=0;
tran_5(sus_trans(sus_trans(:,2)==1&sum(sus_trans(:,7:12)>0,2)==5,6:12)==0)=0;

[S1_1,S2_1]=sortSequenece1(tran_1);
[S1_2,S2_2]=sortSequenece2(tran_2);
[S1_3,S2_3]=sortSequenece3(tran_3);
[S1_4,S2_4]=sortSequenece4(tran_4);
[S1_5,S2_5]=sortSequenece5(tran_5);
fh=figure('Color','w','Position',[100,100,250,400]);
colormap([[[0:0.05:0.9;0:0.05:0.9]' ones(19,1)];[1,1,1];[ones(19,1) [0.9:-0.05:0;0.9:-0.05:0]']]);
imagesc([S1_5;S1_4;S1_3;S1_2;S1_1;S2_1;S2_2;S2_3;S2_4;S2_5])
caxis([-1 1])
colorbar
exportgraphics(fh,fullfile('selectivity_sequence.pdf'),'ContentType','vector')

function [S1,S2]=sortSequenece1(tran_1)
[m1,n1]=sort(tran_1(:,2),'descend');
[m2,n2]=sort(tran_1(:,3),'descend');
[m3,n3]=sort(tran_1(:,4),'descend');
[m4,n4]=sort(tran_1(:,5),'descend');
[m5,n5]=sort(tran_1(:,6),'descend');
[m6,n6]=sort(tran_1(:,7),'descend');
S1=[tran_1(n1(m1>0),:);tran_1(n2(m2>0),:);tran_1(n3(m3>0),:);tran_1(n4(m4>0),:);tran_1(n5(m5>0),:);tran_1(n6(m6>0),:)];
S2=[tran_1(flipud(n1(m1<0)),:);tran_1(flipud(n2(m2<0)),:);tran_1(flipud(n3(m3<0)),:);tran_1(flipud(n4(m4<0)),:);tran_1(flipud(n5(m5<0)),:);tran_1(flipud(n6(m6<0)),:)];
end
function [S1,S2]=sortSequenece2(tran_1)
[m1,n1]=sort(sum(tran_1(:,2:3),2),'descend');
t2=tran_1(n1(m1==0),:);
[m2,n2]=sort(sum(t2(:,3:4),2),'descend');
t3=t2(n2(m2==0),:);
[m3,n3]=sort(sum(t3(:,4:5),2),'descend');
t4=t3(n3(m3==0),:);
[m4,n4]=sort(sum(t4(:,5:6),2),'descend');
t5=t4(n4(m4==0),:);
[m5,n5]=sort(sum(t5(:,6:7),2),'descend');
S1=[tran_1(n1(m1>0),:);t2(n2(m2>0),:);t3(n3(m3>0),:);t4(n4(m4>0),:);t5(n5(m5>0),:)];
S2=[tran_1(flipud(n1(m1<0)),:);t2(flipud(n2(m2<0)),:);t3(flipud(n3(m3<0)),:);t4(flipud(n4(m4<0)),:);t5(flipud(n5(m5<0)),:)];
end
function [S1,S2]=sortSequenece3(tran_1)
[m1,n1]=sort(sum(tran_1(:,2:4),2),'descend');
t2=tran_1(n1(m1==0),:);
[m2,n2]=sort(sum(t2(:,3:5),2),'descend');
t3=t2(n2(m2==0),:);
[m3,n3]=sort(sum(t3(:,4:6),2),'descend');
t4=t3(n3(m3==0),:);
[m4,n4]=sort(sum(t4(:,5:7),2),'descend');
S1=[tran_1(n1(m1>0),:);t2(n2(m2>0),:);t3(n3(m3>0),:);t4(n4(m4>0),:)];
S2=[tran_1(flipud(n1(m1<0)),:);t2(flipud(n2(m2<0)),:);t3(flipud(n3(m3<0)),:);t4(flipud(n4(m4<0)),:)];
end
function [S1,S2]=sortSequenece4(tran_1)
[m1,n1]=sort(sum(tran_1(:,2:5),2),'descend');
t2=tran_1(n1(m1==0),:);
[m2,n2]=sort(sum(t2(:,3:6),2),'descend');
t3=t2(n2(m2==0),:);
[m3,n3]=sort(sum(t3(:,4:7),2),'descend');
S1=[tran_1(n1(m1>0),:);t2(n2(m2>0),:);t3(n3(m3>0),:)];
S2=[tran_1(flipud(n1(m1<0)),:);t2(flipud(n2(m2<0)),:);t3(flipud(n3(m3<0)),:)];
end
function [S1,S2]=sortSequenece5(tran_1)
[m1,n1]=sort(sum(tran_1(:,2:6),2),'descend');
t2=tran_1(n1(m1==0),:);
[m2,n2]=sort(sum(t2(:,3:7),2),'descend');
S1=[tran_1(n1(m1>0),:);t2(n2(m2>0),:)];
S2=[tran_1(flipud(n1(m1<0)),:);t2(flipud(n2(m2<0)),:)];
end