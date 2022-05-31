function [out]=markLPerf(facSeq,crierian,lick_crierian,trial_min,task)
switch(task)
    case 'dualtask'
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
    case 'DPA'
        i=trial_min;
        facSeq_WT=zeros(length(facSeq),1);
        facSeq(:,end+1)=xor(facSeq(:,3)==facSeq(:,4),facSeq(:,5)>0);
        while i<=length(facSeq)
            goodOff=nnz(xor(facSeq(i-trial_min+1:i,3)==facSeq(i-trial_min+1:i,4),facSeq(i-trial_min+1:i,5)>0));
            lickoff=nnz(facSeq(i-trial_min+1:i,3)~=facSeq(i-trial_min+1:i,4)& facSeq(i-trial_min+1:i,5)>0 );
            if goodOff>=crierian*trial_min && lickoff >= (0.5*lick_crierian*trial_min) %.trial_min correct rate
                facSeq_WT(i-trial_min+1:i,1)=1;
            end
            i=i+1;
        end
        out=[facSeq,facSeq_WT];
end

% 
% function [out]=markLPerf(facSeq,trial_min,lick_criteria,Perf_criteria)
% i=trial_min;
% facSeq_WT=zeros(length(facSeq),1);
% while i<=length(facSeq)
%     goodOff=nnz(xor(facSeq(i-trial_min+1:i,5)==facSeq(i-trial_min+1:i,6) , facSeq(i-trial_min+1:i,7)>0));
%     lickOff=nnz(facSeq(i-trial_min+1:i,5)~=facSeq(i-trial_min+1:i,6)&facSeq(i-trial_min+1:i,7)>0);
%     if goodOff>=Perf_criteria*trial_min && lickOff>=0.5*lick_criteria*trial_min %.60 correct rate
%         facSeq_WT(i-trial_min+1:i,1)=1;
%     end
%     i=i+1;
% end
% out=[facSeq,facSeq_WT];
% end
