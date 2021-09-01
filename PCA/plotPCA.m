close
clc
clear
path='D:\pixel-optogenetic\PCA\late3s\';
comp_sel=h5read(fullfile(path,'PCA_comp.hdf5'),'/comp_sel');
comp_sel_repeat=h5read(fullfile(path,'PCA_comp.hdf5'),'/comp_nonsel_repeat');

comp_sel=mean(comp_sel_repeat,4);

fh=figure('Color','w','Position',[50,50,800,600]);
view([14,7])
c={'b-','b--','k-','k--'};
hold on
% for j=9:52
for j=1:size(mean(comp_sel_repeat,4),2)/4-1
    plot3(comp_sel(1,j:j+1),comp_sel(2,j:j+1),comp_sel(3,j:j+1),'Color','b','LineStyle','--','Marker','.','LineWidth',2,'MarkerSize',5)
    plot3(comp_sel(1,size(comp_sel,2)/4+j:size(comp_sel,2)/4+j+1),comp_sel(2,size(comp_sel,2)/4+j:size(comp_sel,2)/4+j+1),comp_sel(3,size(comp_sel,2)/4+j:size(comp_sel,2)/4+j+1),'Color','b','LineStyle',':','Marker','.','LineWidth',2,'MarkerSize',5)
    plot3(comp_sel(1,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),comp_sel(2,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),comp_sel(3,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),'Color','k','LineStyle','--','Marker','.','LineWidth',2,'MarkerSize',5)
    plot3(comp_sel(1,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),comp_sel(2,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),comp_sel(3,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),'Color','k','LineStyle',':','Marker','.','LineWidth',2,'MarkerSize',5)    
end

legend('S1-laser on','S2-laser on','S1-laser off','S2-laser off')
marker{5}='^';marker{41}='+';marker{33}='*';marker{34}='*';
for j=[5] %laser phase
    plot3(comp_sel(1,j),comp_sel(2,j),comp_sel(3,j),'Color','r','LineStyle','none','Marker',marker{j},'LineWidth',2)
    plot3(comp_sel(1,size(comp_sel,2)/4+j),comp_sel(2,size(comp_sel,2)/4+j),comp_sel(3,size(comp_sel,2)/4+j),'Color','r','LineStyle','none','Marker',marker{j},'LineWidth',2)
    plot3(comp_sel(1,2*size(comp_sel,2)/4+j),comp_sel(2,2*size(comp_sel,2)/4+j),comp_sel(3,2*size(comp_sel,2)/4+j),'Color','r','LineStyle','none','Marker',marker{j},'LineWidth',2)
    plot3(comp_sel(1,3*size(comp_sel,2)/4+j),comp_sel(2,3*size(comp_sel,2)/4+j),comp_sel(3,3*size(comp_sel,2)/4+j),'Color','r','LineStyle','none','Marker',marker{j},'LineWidth',2)
end
% 
xlabel('PC1')
ylabel('PC2')
zlabel('PC3') 


exportgraphics(fh,fullfile(path,'PCA_nonsel.pdf'),'ContentType','vector')
% 
close

fh=figure('Color','w','Position',[50,50,600,600]);
view([-13,7])
xlabel('PC1')
ylabel('PC2')
zlabel('PC3') 
hold on 
v = VideoWriter(fullfile(path,'pca_nonsel.avi'),'Uncompressed AVI');
v.FrameRate=4;
open(v);

for j=2:size(mean(comp_sel_repeat,4),2)/4
  plot3(comp_sel(1,j:j+1),comp_sel(2,j:j+1),comp_sel(3,j:j+1),'Color','b','LineStyle','--','Marker','.','LineWidth',2,'MarkerSize',5)
    plot3(comp_sel(1,size(comp_sel,2)/4+j:size(comp_sel,2)/4+j+1),comp_sel(2,size(comp_sel,2)/4+j:size(comp_sel,2)/4+j+1),comp_sel(3,size(comp_sel,2)/4+j:size(comp_sel,2)/4+j+1),'Color','b','LineStyle',':','Marker','.','LineWidth',2,'MarkerSize',5)
    plot3(comp_sel(1,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),comp_sel(2,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),comp_sel(3,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),'Color','k','LineStyle','--','Marker','.','LineWidth',2,'MarkerSize',5)
    plot3(comp_sel(1,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),comp_sel(2,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),comp_sel(3,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),'Color','k','LineStyle',':','Marker','.','LineWidth',2,'MarkerSize',5)    
    if j==33 
       plot3(comp_sel(1,j:j+1),comp_sel(2,j:j+1),comp_sel(3,j:j+1),'Color','r','LineStyle','none','Marker','*','MarkerSize',15)
       plot3(comp_sel(1,1*size(comp_sel,2)/4+j:1*size(comp_sel,2)/4+j+1),comp_sel(2,1*size(comp_sel,2)/4+j:1*size(comp_sel,2)/4+j+1),comp_sel(3,1*size(comp_sel,2)/4+j:1*size(comp_sel,2)/4+j+1),'Color','r','LineStyle','none','Marker','*','MarkerSize',15)
       plot3(comp_sel(1,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),comp_sel(2,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),comp_sel(3,2*size(comp_sel,2)/4+j:2*size(comp_sel,2)/4+j+1),'Color','r','LineStyle','none','Marker','*','MarkerSize',15)
       plot3(comp_sel(1,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),comp_sel(2,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),comp_sel(3,3*size(comp_sel,2)/4+j:3*size(comp_sel,2)/4+j+1),'Color','r','LineStyle','none','Marker','*','MarkerSize',15)
    end
    
    if j<=12
        title('baseline')
        
    elseif j<=16
        title('Sample odor')
    elseif j<=40
        title('delay')
    elseif j<=44
        title('Test odor')
    elseif j<=48
        title('response window')
    else
        title('ITI')
    end
    frame=getframe(gcf);
    writeVideo(v,frame);
end

close(v)

%%
close
dist_sel_on=h5read(fullfile(path,'pca_dist_100.hdf5'),'/dist_sel_on')';
dist_sel_off=h5read(fullfile(path,'pca_dist_100.hdf5'),'/dist_sel_off')';
dist_nonsel_on=h5read(fullfile(path,'pca_dist_100.hdf5'),'/dist_nonsel_on')';
dist_nonsel_off=h5read(fullfile(path,'pca_dist_100.hdf5'),'/dist_nonsel_off')';

fh=figure('Color','w','Position',[50,50,400,200]);

subplot(1,2,1)
hold on
m=mean(dist_sel_on,1);
s=std(dist_sel_on)/sqrt(size(dist_sel_on,1));
plot(1:length(m),m,'b')
a = fill([1:length(m), fliplr(1:length(m))], [m+s, fliplr(m-s)],'b', 'edgecolor','none');
alpha(a,0.2);

m=mean(dist_sel_off,1);
s=std(dist_sel_off)/sqrt(size(dist_sel_off,1));
plot(1:length(m),m,'k')
a = fill([1:length(m), fliplr(1:length(m))], [m+s, fliplr(m-s)],'k', 'edgecolor','none');
alpha(a,0.2);

for i=1:size(dist_sel_on,2)
if ranksum(dist_sel_on(:,i),dist_sel_off(:,i))<0.05/(size(dist_sel_off,1)*length(m))
    plot(i,50,'k*')
end
end
plot([4.5 4.5],[0 50],'k--')
plot([6.5 6.5],[0 50],'k--')
set(gca,'XTick',0.5:4:12,'XTickLabel',{'-1','0','1','','','','5','',' ',' ','',''})
xlabel('time')
ylabel('S1-S2 PC distance')
ylim([0 50])
% plot([13 13],[0 50],'k--')
% plot([16 16],[0 50],'k--')
% plot([33 33],[0 50],'k--')
% plot([34 34],[0 50],'k--')
% plot([40 40],[0 50],'k--')
% plot([44 44],[0 50],'k--')

% set(gca,'XTick',9:4:48,'XTickLabel',{' ','0',' ','','','','5','',' ',' ','',''})
% xlim([9 52])
% legend('laser on','laser on','laser off','laser off')
box off

subplot(1,2,2)
hold on
m=mean(dist_nonsel_on,1);
s=std(dist_nonsel_on)/sqrt(size(dist_nonsel_on,1));
plot(1:length(m),m,'b')
a = fill([1:length(m), fliplr(1:length(m))], [m+s, fliplr(m-s)],'b', 'edgecolor','none');
alpha(a,0.2);

m=mean(dist_nonsel_off,1);
s=std(dist_nonsel_off)/sqrt(size(dist_nonsel_off,1));
plot(1:length(m),m,'k')
a = fill([1:length(m), fliplr(1:length(m))], [m+s, fliplr(m-s)],'k', 'edgecolor','none');
alpha(a,0.2);

for i=1:size(dist_sel_on,2)
if signrank(dist_nonsel_on(:,i),dist_nonsel_off(:,i))<0.05/(size(dist_sel_off,1)*length(m))
    plot(i,50,'k*')
end
end
plot([4.5 4.5],[0 50],'k--')
plot([6.5 6.5],[0 50],'k--')
set(gca,'XTick',0.5:4:12,'XTickLabel',{'-1','0','1','','','','5','',' ',' ','',''})
xlabel('time')
ylabel('S1-S2 PC distance')

% plot([13 13],[0 120],'k--')
% plot([16 16],[0 120],'k--')
% plot([33 33],[0 120],'k--')
% plot([34 34],[0 120],'k--')
% plot([40 40],[0 120],'k--')
% plot([44 44],[0 120],'k--')

% set(gca,'XTick',9:4:48,'XTickLabel',{' ','0',' ','','','','5','',' ',' ','',''})
% xlim([9 52])

exportgraphics(fh,fullfile(path,'PCA_dist.pdf'),'ContentType','vector')

