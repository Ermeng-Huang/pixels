function plotLine(Data,c,l)
    m=mean(Data,1);
    plot(1:size(m,2),m,c,'LineStyle',l);
%     s=std(Data)/sqrt(size(Data,1));
%     fill([1:size(m,2),fliplr(1:size(m,2))],[m+s,fliplr(m-s)],c,'EdgeColor','none','FaceAlpha',0.2);
    ci=bootci(1000,@(x) mean(x),Data);
    fill([1:size(m,2),fliplr(1:size(m,2))],[ci(1,:),fliplr(ci(2,:))],c,'EdgeColor','none','FaceAlpha',0.2);
end