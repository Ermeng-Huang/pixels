%%% Performance
clear
clc
close all

homedir='X:\';
trial_min=120;
crierian=0.7;
lick_crierian=0.8;
load(fullfile(homedir,'session_list.mat'))
%% Reaction Time
file_id=unique(floor(h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/cluster_id')/100000));
s=session(file_id);
for i=1:size(s,1)    
    [event{i},lick{i}]=ExtractEvent(util.ser2mat2(fullfile(homedir,'behavior',sprintf('%s.ser',s{i}(1:12)))),crierian,lick_crierian,trial_min);
    rt{i}=RT(event{i},lick{i});
end

%%

for i=1:size(s,1)
    % hit
    rt_m(1,i)=mean(rt{i}(event{i}(:,6)==7&event{i}(:,7)==-1&event{i}(:,14)==1),1)/1000;
    rt_m(2,i)=mean(rt{i}(event{i}(:,6)==7&event{i}(:,7)~=-1&event{i}(:,14)==1),1)/1000;
    % false
    rt_m(3,i)=mean(rt{i}(event{i}(:,6)==4&event{i}(:,7)==-1),1)/1000;
    rt_m(4,i)=mean(rt{i}(event{i}(:,6)==4&event{i}(:,7)~=-1),1)/1000;
end

fh=figure('Color','w','Position',[100,100,215,250]);
subplot(2,1,1)
hold on
bar(1:2,mean(rt_m(1:2,:),2),'w')
errorbar(1:2,mean(rt_m(1:2,:),2),std(rt_m(1:2,:),1,2)/sqrt(i),'linestyle','none','color','k')
set(gca,'XTick',1.5,'XTickLabel',{'hit trials'})
ylabel('Reaction Time (s)')

subplot(2,1,2)
hold on
bar(1:2,mean(rt_m(3:4,:),2),'w')
errorbar(1:2,mean(rt_m(3:4,:),2),std(rt_m(3:4,:),1,2)/sqrt(i),'linestyle','none','color','k')
set(gca,'XTick',1.5,'XTickLabel',{'false trials'})
ylabel('Reaction Time (s)')

% ylabel('Trial ID #')
exportgraphics(fh,'E:\dualtask\ReactionTime.pdf','ContentType','vector');
ranksum(rt_m(1,:),rt_m(2,:))
ranksum(rt_m(3,:),rt_m(4,:))

return
%% lick raster
trials_choose=cell2mat(cellfun(@(x)x(x(:,14)==1,:),event,'UniformOutput',false)');
lick_choose=cellfun(@(x,y)y(x(:,14)==1)',event,lick,'UniformOutput',false);
lick_choose=[lick_choose{:}]';

fh=figure('Color','w','Position',[100,100,320,250]);
hold on
b=80;
num=24;
for i=1:num
    if trials_choose(i+b,7)==-1 % blank-DPA
        fill([-300,1400,1400,-300],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'r','FaceAlpha',0.2);
    elseif trials_choose(i+b,7)~=-1 % dualtask-DPA
        fill([-300,1400,1400,-300],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'b','FaceAlpha',0.2);
    end
    try
        plot([lick_choose{i+b}/10 lick_choose{i+b}/10], [num+1-i-0.5 num+1-i+0.5],'k');

    catch
    end    
    
    % DR response
    if trials_choose(i+b,12)==7 %hit
       fill([1400,1450,1450,1400],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'r');
    elseif trials_choose(i+b,12)==5 %CR
       fill([1400,1450,1450,1400],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'k');  
    elseif trials_choose(i+b,12)==4 %F
       fill([1400,1450,1450,1400],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'b');
    end

%     % DPA sample
%     if trials_choose(i+b,2)==9
%        fill([1400,1450,1450,1400],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'r');
%     elseif trials_choose(i+b,2)==10
%        fill([1400,1450,1450,1400],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'b');
%     end

    % DPA response
    if trials_choose(i+b,6)==7 %hit
       fill([1500,1550,1550,1500],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'r');
    elseif trials_choose(i+b,6)==5 %CR
       fill([1500,1550,1550,1500],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'k');  
    elseif trials_choose(i+b,6)==4 %F
       fill([1500,1550,1550,1500],[num+1-i-0.5,num+1-i-0.5,num+1-i+0.5,num+1-i+0.5], 'b');
    end
    

   
end


plot([0 0],[0 i+0.5],'g--')
plot([300 300],[0 i+0.5],'g--')
plot([600 600],[0 i+0.5],'g--')
plot([900 900],[0 i+0.5],'g--')

xlim([-100 1550])
ylim([0+0.5 i+0.5])
set(gca,'XTick',0:500:1400,'XTickLabel',0:5:10)
xlabel('Time (s)')
ylabel('Trial ID #')
exportgraphics(fh,'E:\dualtask\lick_raster_24trials.pdf','ContentType','vector');

%% lick psth of different trial type
f=figure('Color','w','Position',[100,100,300,650]);
%%%Only use hit trials
% DistactorNo
subplot(3,1,1)
hold on
util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-1000:250:13000)/4,z,'UniformOutput',false))...
    ,cellfun(@(x,y)y(x(:,6)==7&x(:,8)==-1&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','-');
% util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-3000:250:11000)/4,z,'UniformOutput',false))...
%     ,cellfun(@(x,y)y(x(:,6)==4&x(:,8)==-1&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','-');
% DistactorGo
arrayfun(@(x)plot([x,x],[0,1],'k--'),[4.5,8.5,40.5,44.5])
arrayfun(@(x)plot([x,x],[0,1],'k:'),[16.5,20.5,28.5,32.5])
set(gca,'xlim',[0,56],'ylim',[0,1],'xTick',4.5:20:56,'xTickLabel',{'0','5','10'})
ylabel('Lick Rate (Hz)')
title('distractorNo')

subplot(3,1,2)
hold on
util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-1000:250:13000)/4,z,'UniformOutput',false))...
    ,cellfun(@(x,y)y(x(:,6)==7&x(:,8)==66&x(:,12)==7&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','-');
% util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-3000:250:11000)/4,z,'UniformOutput',false))...
%     ,cellfun(@(x,y)y(x(:,6)==4&x(:,8)==66&x(:,12)==7&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','-');
util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-1000:250:13000)/4,z,'UniformOutput',false))...
    ,cellfun(@(x,y)y(x(:,6)==7&x(:,8)==66&x(:,12)==6&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','--');
% util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-3000:250:11000)/4,z,'UniformOutput',false))...
%     ,cellfun(@(x,y)y(x(:,6)==4&x(:,8)==66&x(:,12)==6&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','--');
arrayfun(@(x)plot([x,x],[0,1],'k--'),[4.5,8.5,40.5,44.5])
arrayfun(@(x)plot([x,x],[0,1],'k:'),[16.5,20.5,28.5,32.5])
set(gca,'xlim',[0,56],'ylim',[0,1],'xTick',4.5:20:56,'xTickLabel',{'0','5','10'})
ylabel('Lick Rate (Hz)')
title('distractorGo')

subplot(3,1,3)
hold on
util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-1000:250:13000)/4,z,'UniformOutput',false))...
    ,cellfun(@(x,y)y(x(:,6)==7&x(:,8)==64&x(:,12)==5&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','-');
% util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-3000:250:11000)/4,z,'UniformOutput',false))...
%     ,cellfun(@(x,y)y(x(:,6)==4&x(:,8)==64&x(:,12)==5&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','-');
util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-1000:250:13000)/4,z,'UniformOutput',false))...
    ,cellfun(@(x,y)y(x(:,6)==7&x(:,8)==64&x(:,12)==4&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','--');
% util.plotLine(cell2mat(cellfun(@(z)cell2mat(cellfun(@(s)histcounts(double(s),-3000:250:11000)/4,z,'UniformOutput',false))...
%     ,cellfun(@(x,y)y(x(:,6)==4&x(:,8)==64&x(:,12)==4&x(:,14)==1),event,lick,'UniformOutput',false),'UniformOutput',false)'),'k','--');
arrayfun(@(x)plot([x,x],[0,1],'k--'),[4.5,8.5,40.5,44.5])
arrayfun(@(x)plot([x,x],[0,1],'k:'),[16.5,20.5,28.5,32.5])
set(gca,'xlim',[0,56],'ylim',[0,1],'xTick',4.5:20:56,'xTickLabel',{'0','5','10'})
ylabel('Lick Rate (Hz)')
title('distractorNoGo')
% exportgraphics(f,fullfile(homedir,'Performance','LickInTreeCondition.png'))
% exportgraphics(f,fullfile(homedir,'Performance','LickInTreeCondition.pdf'))
function [out]=markLPerf(facSeq,crierian,lick_crierian,trial_min)
i=trial_min;
facSeq_WT=zeros(length(facSeq),1);
while i<=length(facSeq)
    goodOff=nnz(xor(facSeq(i-trial_min+1:i,5)==facSeq(i-trial_min+1:i,6),facSeq(i-trial_min+1:i,7)>0)& facSeq(i-trial_min+1:i,9)==-1);
%     lickoff=nnz(facSeq(i-trial_min+1:i,5)~=facSeq(i-trial_min+1:i,6)& facSeq(i-trial_min+1:i,7)>0); %% exclude unlicked in distractor No trial
    lickoff=nnz(facSeq(i-trial_min+1:i,5)~=facSeq(i-trial_min+1:i,6)& facSeq(i-trial_min+1:i,7)>0 & facSeq(i-trial_min+1:i,9)==-1);
    if goodOff>=crierian*trial_min/3 && lickoff > (0.5*lick_crierian*trial_min)/3 %.trial_min correct rate
        facSeq_WT(i-trial_min+1:i,1)=1;
    end
    i=i+1;
end
out=[facSeq,facSeq_WT];
end

function [event,lick_trial]=ExtractEvent(in,crierian,lick_crierian,trial_min)
% DPA
Sample=in((in(:,3)==10&in(:,4)==0)|(in(:,3)==10&in(:,4)==100)|(in(:,3)==9&in(:,4)==1),[1,3]);
Test=in((in(:,3)==10&in(:,4)==2)|(in(:,3)==9&in(:,4)==3),[1,3]);
Sample(~arrayfun(@(x)any(all((Test(:,1)-9500)<x & (Test(:,1)-8500)>x,2),1),Sample(:,1)),:)=[];
Response=in((in(:,3)==4&in(:,4)==2)|(in(:,3)==5&in(:,4)==2)|(in(:,3)==6&in(:,4)==2)|(in(:,3)==7&in(:,4)==2)...
    |(in(:,3)==4&in(:,4)==1)|(in(:,3)==5&in(:,4)==1)|(in(:,3)==6&in(:,4)==1)|(in(:,3)==7&in(:,4)==1),[1,3]);
Response(~arrayfun(@(x)any(all((Test(:,1)-500)<x & (Test(:,1)+2500)>x,2),1),Response(:,1)),:)=[];
%%% unknown question num(Sample)~= num(Test)and num(Response) file 35 39
if min([size(Sample,1),size(Test,1),size(Response,1)]) ~= max([size(Sample,1),size(Test,1),size(Response,1)])
    
    if min([size(Sample,1),size(Test,1),size(Response,1)])==size(Sample,1)
        miss.id=find(~arrayfun(@(x)any(all((Sample(:,1)+9500)>x & (Sample(:,1)+8500)<x,2),1),Test(:,1)));
        Test(miss.id,:)=[];
        Response(miss.id,:)=[];
    elseif min([size(Sample,1),size(Test,1),size(Response,1)])==size(Test,1)
        miss.id=find(~arrayfun(@(x)any(all((Test(:,1)-8500)>x & (Test(:,1)-9500)<x,2),1),Sample(:,1)));
        Sample(miss.id,:)=[];
        Response(miss.id,:)=[];
    else
        miss.id=find(~arrayfun(@(x)any(all((Response(:,1)-2500)<x & (Response(:,1)-500)>x,2),1),Test(:,1)));
        Sample(miss.id,:)=[];
        Test(miss.id,:)=[];        
    end
end

% DR
DRSample=in((in(:,3)==66&in(:,4)==4)|(in(:,3)==64&in(:,4)==5),[1,3]);
DRTest=in((in(:,3)==64&in(:,4)==6),[1,3]);
DRSample(~arrayfun(@(x)any(all((DRTest(:,1)-3500)<x & (DRTest(:,1)-2500)>x,2),1),DRSample(:,1)),:)=[];
DRResponse=in((in(:,3)==4&in(:,4)==3)|(in(:,3)==5&in(:,4)==3)|(in(:,3)==6&in(:,4)==3)|(in(:,3)==7&in(:,4)==3),[1,3]);
DRResponse(~arrayfun(@(x)any(all((DRTest(:,1)-500)<x & (DRTest(:,1)+2500)>x,2),1),DRResponse(:,1)),:)=[];
%%% unknown question num(DRResponse)~= num(DRSample)
if size(DRSample,1)~=size(DRResponse,1)
    if size(DRSample,1)<size(DRResponse,1)
        arrayfun(@(x)any(all((DRSample(:,1)+3500)<x & (DRSample(:,1)+5500)>x,2),1),DRResponse(:,1))
    else
       miss.t=DRSample(~arrayfun(@(x)any(all((DRResponse(:,1)-3500)>x & (DRResponse(:,1)-5500)<x,2),1),DRSample(:,1)),1); 
       miss.id=find(~arrayfun(@(x)any(all((DRResponse(:,1)-3500)>x & (DRResponse(:,1)-5500)<x,2),1),DRSample(:,1)));
       if DRSample(~arrayfun(@(x)any(all((DRResponse(:,1)-3500)>x & (DRResponse(:,1)-5500)<x,2),1),DRSample(:,1)),2)==64 %no-go           
           if any(in(in(:,1)>(miss.t+3100)&in(:,1)<(miss.t+5200),3)==0) % no lick
               DRResponse=[DRResponse(1:miss.id,:);[-1000000,5];DRResponse(1:miss.id,:)];
           else
               DRResponse=[DRResponse(1:miss.id-1,:);[-1000000,4];DRResponse(miss.id:end,:)];           
           end
       else %go
           if any(in(in(:,1)>(miss.t+3100)&in(:,1)<(miss.t+5200),3)==0) % no lick
               DRResponse=[DRResponse(1:miss.id-1,:);[-1000000,6];DRResponse(miss.id:end,:)];
           else
               DRResponse=[DRResponse(1:miss.id-1,:);[-1000000,7];DRResponse(miss.id:end,:)];           
           end
       end
    end    
end

if size(DRSample,1)==size(Sample,1)
    event=[Sample,Test,Response,DRSample,DRTest,DRResponse];
else
    for i=1:size(Sample,1)
        if nnz((DRSample(:,1)-Sample(i,1))<3100&(DRSample(:,1)-Sample(i,1))>2900)==0 % blank trial
            DRSample_new(i,:)=int32(-ones(1,6));
        else
            trial_id=find((DRSample(:,1)-Sample(i,1))<3100&(DRSample(:,1)-Sample(i,1))>2900);
            DRSample_new(i,:)=[DRSample(trial_id,:),DRTest(trial_id,:),DRResponse(trial_id,:)];     
            
        end
    end
   
    event=[Sample,Test,Response,DRSample_new];
end
lick=in(in(:,3)==0&in(:,4)==1,1);
lick_trial=arrayfun(@(x)lick(lick>(x-3000)&lick<(x+11000),1)-x,Sample(:,1),'UniformOutput',false);

event(:,13)=0;
event(event(:,6)==5|event(:,6)==7,13)=1;
event(:,14)=0;
i=trial_min;
% only use blank trial
while i<=size(event,1)
    goodOff=nnz(event(i-trial_min+1:i,13)==1&event(i-trial_min+1:i,7)==-1); %7-hit,5-correct
    lickoff=nnz(event(i-trial_min+1:i,6)==7&event(i-trial_min+1:i,7)==-1);
    if goodOff>=crierian*trial_min/3 && lickoff > (0.5*lick_crierian*trial_min)/3 %.trial_min correct rate
        event(i-trial_min+1:i,14)=1;
    end
    i=i+1;
end



% if all(cell2mat(arrayfun(@(y)arrayfun(@(x)nnz(event(:,2)==y&event(:,8)==x&all(event(:,13:14),2))>5,[-1,66,64]),[9,10],'UniformOutput',false)))% trial number of each condition
%     availed=true;
% else
%     availed=false;
% end
end

function out=RT(event,lick)
test_onset=double(event(:,3)-event(:,1));
for i=1:size(event,1)
    try
        out(i,1)=double(lick{i}(find(lick{i}>test_onset(i)&lick{i}<test_onset(i)+2000,1)))-test_onset(i);
    catch
        out(i,1)=10000;
    end
end
end
