clc;
clear;
homedir='D:\pixel-optogenetic\';
%% Find ring according to functional coupling
if ~isfile(fullfile(homedir,'xcorr','rings_bz.mat'))
    load(fullfile(homedir,'states_conn_bz_0314.mat'))
    for sess=1:size(sig_con,1)
        for ring_size=3:5
            rings{sess,ring_size-2}=unique(flexsort(find_rings_bz(sig_con{sess,1},ring_size)),'row');
        end
        % shuffled
        
    end
    save(fullfile(homedir,'xcorr','rings_bz.mat'),'rings');
else
    load(fullfile(homedir,'xcorr','rings_bz.mat'));
end

%% shuffled
rings=[];
load(fullfile(homedir,'states_conn_bz_0314.mat'))
for rpt=1:100    
    for sess=1:size(sig_con,1)
        shufs(sess,rpt)=shuffle_conn_onerpt(sig_con(sess,:));
        for ring_size=3:5
            rings{sess,ring_size-2,rpt}=unique(flexsort(find_rings_bz(shufs{sess,rpt},ring_size)),'row');            
        end  
    end
end


save(fullfile(homedir,'xcorr','rings_bz_shufs.mat'),'shufs','rings')


function shuf=shuffle_conn_onerpt(sig_con)
for si=1:size(sig_con,1)
    curr_cnt=length(sig_con{si, 1});
    if curr_cnt==0
        continue
    end    
    shuf{si,1}=sig_con{si, 2}.pair_comb(randsample(1:length(sig_con{si, 2}.pair_comb),curr_cnt),:);
end
end
%% function
function out=find_rings_bz(in,msize)
arguments
    in (:,2) double {mustBeNonempty(in)} %sig or pair?
    msize (1,1) double {mustBeMember(msize,3:5)} % size of ring
end
out=[];
pre_unit_set=unique(in(:,1));
for a_first=pre_unit_set(:)'
    seconds=in(in(:,1)==a_first,2);
    for a_second=seconds(:)'
        thirds=in(in(:,1)==a_second,2);
        thirds(thirds==a_first)=[];
        if msize==3
            ring=in(ismember(in(:,1),thirds) & in(:,2)==a_first,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([a_first a_second],size(subgrp,1),1);
                out=[out;subgrp];
            end
        else
            for a_third=thirds(:)'
                fourths=in(in(:,1)==a_third,2);
                fourths(fourths==a_second)=[];
                if msize==4
                    ring=in(ismember(in(:,1),fourths) & in(:,2)==a_first,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring];
                        subgrp(:,1:3)=repmat([a_first a_second a_third],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                elseif msize==5
                    fourths(fourths==a_first)=[];
                    for a_fourth=fourths(:)'
                        fifths=in(in(:,1)==a_fourth, 2);
                        fifths(ismember(fifths,[a_first a_second]))=[];
                        ring=in(ismember(in(:,1),fifths) & in(:,2)==a_first,1);
                        if ~isempty(ring)
                            subgrp=[ring,ring,ring,ring,ring];
                            subgrp(:,1:4)=repmat([a_first a_second a_third a_fourth],size(subgrp,1),1);
                            out=[out;subgrp];
                        end
                    end
                end
            end
        end
    end
end
end

function out=flexsort(in)
arguments
    in double %list of rings
end
out=in;
for i=1:size(out,1)
    [~,sess]=min(out(i,:));
    while sess>1
        out(i,:)=circshift(out(i,:),1);
        [~,sess]=min(out(i,:));
    end
end
end

