# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:31:41 2021

@author: hem
"""
# Calculate projection matrix by AllenInsititude Data

import numpy as np
import matplotlib.pyplot as plt
import warnings
import pandas as pd
import csv
import h5py
# warnings.filterwarnings('ignore')

import sys
sys.path.insert(1,'D:\code\AllenSDK-master')
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

# tell the cache class what resolution (in microns) of data you want to download
mcc = MouseConnectivityCache(resolution=25)

# use the structure tree class to get information about the isocortex structure
structure_tree = mcc.get_structure_tree()
Id = structure_tree.get_structures_by_name(['Agranular insular area'])[0]['id']

# a list of dictionaries containing metadata for non-Cre experiments
experiments = mcc.get_experiments(file_name='non_cre.json',
                                  injection_structure_ids=[Id])
experiment_ids = [ e['id'] for e in experiments ]

with open('reg_L_20210412.csv', 'r') as f:    
    l = list(csv.reader(f))[1:]
    cortex_phy = [list(i) for i in zip(*l)]
      
# projection_id=structure_tree.get_structures_by_name(cortex_phy[0])[0]['id']
projection_structure_ids=cortex_phy[1]     
# for i in range(len(cortex_phy[0])): 
#     cortex[i] = structure_tree.get_structures_by_name(cortex_phy[0][i])[0] 
#     # ctx_children = structure_tree.child_ids( [isocortex['id']] )[0]

# # download the projection density volume for one of the experiments
# pd = mcc.get_projection_density(experiments[0]['id'])

# download the projection density volume for all experiments
pm = mcc.get_projection_matrix(experiment_ids = experiment_ids, 
                               projection_structure_ids = projection_structure_ids,
                               hemisphere_ids= [2], # right hemisphere, ipsilateral
                               parameter = 'projection_density')

row_labels = pm['rows'] # these are just experiment ids
column_labels = [ c['label'] for c in pm['columns'] ] 
matrix = pm['matrix']

with h5py.File('AllenData_Projection_20210413.hdf5', "w") as fw:
        fw.create_dataset('projection_density', data=matrix.astype('double'))   
        fw.create_dataset('projection_structure',data=np.array(column_labels).astype('S10'))    
        fw.create_dataset('experiment_id', data=np.array(row_labels).astype('S9'))
        
        
fig, ax = plt.subplots(figsize=(15,15))
heatmap = ax.pcolor(matrix, cmap=plt.cm.afmhot)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(matrix.shape[1])+0.5, minor=False)
ax.set_yticks(np.arange(matrix.shape[0])+0.5, minor=False)

ax.set_xlim([0, matrix.shape[1]])
ax.set_ylim([0, matrix.shape[0]])          

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(column_labels,ratation=45,frontsize='x-large', minor=False)
ax.set_yticklabels(row_labels, minor=False)
plt.show()
fig.savefig('ProjectionByAllen' + ".png", dpi=300, bbox_inches="tight")
# summary_structures = structure_tree.get_structures_by_set_id([167587189])
# summary=pd.DataFrame(summary_structures)