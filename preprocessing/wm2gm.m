clearvars -Except av st
% white matter to grey matter
load('D:\histology\dualtask\sucoords.mat');
annotation_volume_location = 'D:\histology\SmartTrack-Histology-master\annotation_volume_10um_by_index.npy';
structure_tree_location = 'D:\histology\SmartTrack-Histology-master\structure_tree_safe_2017.csv';
% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end
regs = deblank(h5read('Selectivity_0602.hdf5','/reg'));
acronym_list = regs(:,7);
w2g=zeros(size(regs,1),1);
% distance queried for confidence metric 
probe_radius = 50; % 500um
tic
for i = 1:length(coord)
    if all(isnan(coord(i,:))) % track is losting
        continue;
    end
    ann = av(round(coord(i,1)),round(coord(i,2)),round(coord(i,3)));
    if strcmp(regs(i,2),'VS') || strcmp(regs(i,2),'fiber tracts') 
%     if strcmp(regs(:,2),'grey') || isAreaOrContains(st.acronym(ann), 'fiber tracts', st)
        radius = 1;    
        while radius <= probe_radius
            [x,y,z] = ndgrid(round(coord(i,1))-radius:round(coord(i,1))+radius,...
                round(coord(i,2))-radius:round(coord(i,2))+radius,round(coord(i,3))-radius:round(coord(i,3))+radius);
            mask = arrayfun(@(x,y,z) norm([x,y,z]-coord(i,:))<=radius,x(:),y(:),z(:));
            anns = arrayfun(@(x,y,z) av(x,y,z),x(mask),y(mask),z(mask));
            acronyms = st.acronym(anns);
            isA = isAreaOrContains(acronyms, 'grey', st);
            if any(isA)
                tbl = tabulate(acronyms(isA));
                [~,I] = max(cell2mat(tbl(:,3)));
                acronym_list{i} = tbl{I,1};
                w2g(i)=1;
                break;
            end
            radius = radius + 1;
        end
    elseif strcmp(regs(i,6),'OLF') && isempty(acronym_list{i})
        radius = 1;        
        while radius <= probe_radius
            [x,y,z] = ndgrid(round(coord(i,1))-radius:round(coord(i,1))+radius,...
                round(coord(i,2))-radius:round(coord(i,2))+radius,round(coord(i,3))-radius:round(coord(i,3))+radius);
            mask = arrayfun(@(x,y,z) norm([x,y,z]-coord(i,:))<=radius,x(:),y(:),z(:));
            anns = arrayfun(@(x,y,z) av(x,y,z),x(mask),y(mask),z(mask));
            acronyms = st.acronym(anns);
            isA = isAreaOrContains(acronyms, 'grey', st) & isAreaOrContains(acronyms, 'OLF', st) & ~strcmp(acronyms,'OLF');
            if any(isA)
                tbl = tabulate(acronyms(isA));
                [~,I] = max(cell2mat(tbl(:,3)));
                acronym_list{i} = tbl{I,1};
                w2g(i)=2;
                break;
            end
            radius = radius + 1;
        end       
    elseif strcmp(regs(i,5),'STR') && isempty(acronym_list{i})
        radius = 1;       
        while radius <= probe_radius
            [x,y,z] = ndgrid(round(coord(i,1))-radius:round(coord(i,1))+radius,...
                round(coord(i,2))-radius:round(coord(i,2))+radius,round(coord(i,3))-radius:round(coord(i,3))+radius);
            mask = arrayfun(@(x,y,z) norm([x,y,z]-coord(i,:))<=radius,x(:),y(:),z(:));
            anns = arrayfun(@(x,y,z) av(x,y,z),x(mask),y(mask),z(mask));
            acronyms = st.acronym(anns);
            isA = isAreaOrContains(acronyms, 'grey', st) & isAreaOrContains(acronyms, 'STR', st) & ~strcmp(acronyms,'STR');
            if any(isA)
                tbl = tabulate(acronyms(isA));
                [~,I] = max(cell2mat(tbl(:,3)));
                acronym_list{i} = tbl{I,1};
                w2g(i)=3;
                break;
            end
            radius = radius + 1;
        end         
%     elseif isAreaOrContains(st.acronym(ann), 'TH', st) && isempty(acronym_list{i})
%         radius = 1;
%         tic
%         while radius <= probe_radius
%             [x,y,z] = ndgrid(round(coord(i,1))-radius:round(coord(i,1))+radius,...
%                 round(coord(i,2))-radius:round(coord(i,2))+radius,round(coord(i,3))-radius:round(coord(i,3))+radius);
%             mask = arrayfun(@(x,y,z) norm([x,y,z]-coord(i,:))<=radius,x(:),y(:),z(:));
%             anns = arrayfun(@(x,y,z) av(x,y,z),x(mask),y(mask),z(mask));
%             acronyms = st.acronym(anns);
%             isA = isAreaOrContains(acronyms, 'grey', st) & isAreaOrContains(acronyms, 'TH', st) & ~strcmp(acronyms,'TH');
%             if any(isA)
%                 tbl = tabulate(acronyms(isA));
%                 [~,I] = max(cell2mat(tbl(:,3)));
%                 acronym_list{i} = tbl{I,1};
%                 w2g(i)=4;
%                 break;
%             end
%             radius = radius + 1;
%         end
%         toc     
    end
end
toc
load('reg_ccfid_map.mat','reg2tree')
reg_new=regs;
for i = find(w2g)'
    reg_new(i,1:length(reg2tree(acronym_list{i})))=reg2tree(acronym_list{i});
end

save('regs_w2g.mat','reg_new','w2g','regs','acronym_list');
