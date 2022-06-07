clearvars -Except av st
close all

% directory of reference atlas files
annotation_volume_location = 'D:\histology\SmartTrack-Histology-master\annotation_volume_10um_by_index.npy';
structure_tree_location = 'D:\histology\SmartTrack-Histology-master\structure_tree_safe_2017.csv';
% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

cid_list=rem(h5read('Selectivity_1129.hdf5','/cluster_id'),10000);
path_list=regexp(h5read('Selectivity_1129.hdf5','/path'),'(\w|\\|-)*','match','once');
track_list=unique(path_list);
load('DepthToTip.mat')
coord=[];
track_table=readtable('D:\histology\dualtask\Neuropixels-dualtask_210906.xlsx','ReadVariableNames',true);

for i=1:length(track_list)
    mouseID=regexp(track_list{i},'(?=M)[A-Za-z0-9]*','match','once');
    date=regexp(track_list{i},'(?=20)[0-9]*','match','once');
    imec=regexp(track_list{i},'(?<=imec)(\d)','match','once');
    temp=track_table(strcmp(track_table.MouseID,mouseID)&track_table.Date==str2double(date)&track_table.ImecID==str2double(imec),:);
    probeID=temp.ProbeID;
    if  isempty(probeID)||isnan(probeID)
        coord = cat(1,coord,nan(nnz(ismember(path_list,track_list{i})),3));
        continue
    end
        
    load(fullfile('D:\histology\dualtask',sprintf('%s',mouseID),'processed/probe_pointselectrode_1.mat'));
    load(fullfile('D:\histology\dualtask',sprintf('%s',mouseID),'results/HistLength.mat'));
    
    ProbeLength=Probe_length(probeID);
    depth_temp=split(temp.Site{1},',');
    depth_record=round(abs(str2num(depth_temp{end})*1000)); % depth by recording
    depth_site=max(DepthToTip(ismember(path_list,track_list{i}))); % max depth by site    
    if depth_record>ProbeLength&&depth_site>depth_record
        ProbeLength=depth_site;
    elseif depth_site<depth_record && abs(depth_record-ProbeLength)>300
        ProbeLength=depth_record;
    end

    realDepth = ProbeLength - DepthToTip(ismember(path_list,track_list{i}));    
    curr_probePoints=pointList.pointList{probeID,1}(:, [3 2 1]);
    % m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
    [m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));

    % ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
    if p(2)<0
        p = -p;
    end
    % determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
    ann = 10;
    out_of_brain = false;
    while ~(ann==1 && out_of_brain) % && distance_stepped > .5*active_probe_length)
        m = m-p; % step 10um, backwards up the track
        ann = av(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
        if strcmp(st.safe_name(ann), 'root')
            % make sure this isn't just a 'root' area within the brain
            m_further_up = m - p*20; % is there more brain 200 microns up along the track?
            ann_further_up = av(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
            if strcmp(st.safe_name(ann_further_up), 'root')
                out_of_brain = true;
            end
        end
    end
    coord = cat(1,coord,[m(1)+p(1)*realDepth/10,m(2)+p(2)*realDepth/10,m(3)+p(3)*realDepth/10]);
end
save('D:\histology\sucoords.mat','coord','-v7.3');       

 
