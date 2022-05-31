function out = mem_type(homedir,in)
cluster_id=h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/cluster_id');
sus_trans=h5read(fullfile(homedir,'Selectivity_1129.hdf5'),'/sust_trans_noPermutaion');
if iscell(in)
    all_su=unique(cell2mat(in(:)));
else
    all_su=unique(in(:));
end
for i=1:numel(all_su)
    if ~any(sus_trans(cluster_id==all_su(i),1:2)) %non-mem
        mem(i)=-1;
    else
        if sus_trans(cluster_id==all_su(i),1)==1  %sustained
            prefer=unique(sus_trans(cluster_id==all_su(i),5:12));
            if prefer(prefer~=0)==1  %s1
                mem(i)=1;
            else
                mem(i)=4;
            end
        elseif sus_trans(cluster_id==all_su(i),2)==1 || sus_trans(cluster_id==all_su(i),2)==3
            prefer=unique(sus_trans(cluster_id==all_su(i),5:12));
            if prefer(prefer~=0)==1 %s1
                mem(i)=2;
            else
                mem(i)=5;
            end
        elseif sus_trans(cluster_id==all_su(i),2)==2
            prefer=unique(sus_trans(cluster_id==all_su(i),13:20));
            if prefer(prefer~=0)==1 %s1
                mem(i)=3;
            else
                mem(i)=6;
            end
        end
    end
end
out=arrayfun(@(x)mem(all_su==x),in);
end
