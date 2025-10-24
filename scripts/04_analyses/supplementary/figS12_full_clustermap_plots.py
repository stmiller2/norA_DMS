#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 17:40:39 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
from plot_utils import spec_cmap, spec_vmin, spec_vmax

#Load data
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='WT')
ligcols = data.columns[:8]

#Get cluster assignments
cluster_assignments = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv', index_col = 'mut')['cluster']
data['cluster'] = cluster_assignments

#Organize data by cluster
min_cluster_size = 3
mask = data['cluster'].value_counts() >= min_cluster_size
clusters = data['cluster'].value_counts()[mask].index
data_by_cluster = {cluster: data[data['cluster'] == cluster] for cluster in clusters}
keys = list(data_by_cluster.keys())
keys.sort()
data_by_cluster = {i: data_by_cluster[i] for i in keys}

#Plot
fig, axes = plt.subplots(nrows=len(clusters), figsize = (4, 12))
cbar_ax = fig.add_axes([0.9, 0.05, .05, .2])
for idx, (cluster, data_cluster) in enumerate(data_by_cluster.items()):
    sns.heatmap(data_cluster[ligcols], cmap=spec_cmap, center = 0, 
                ax=axes[idx], vmin = spec_vmin, vmax = spec_vmax,
                cbar = (idx == 0), cbar_ax = None if idx else cbar_ax,
                cbar_kws = None if idx else {'label': 'Standardized\nfunctional score'})    
    axes[idx].set_yticks([])
    axes[idx].set_ylabel(f'c{int(cluster)} (N = {len(data_cluster)})', 
                         rotation = 0, 
                         fontsize = 12, 
                         ha = 'right')
    axes[idx].yaxis.set_label_coords(-0.01,0.05)
    if idx == 0:
        axes[idx].set_title(f'Euclidean average, 54 clusters\nMinimum cluster size: {min_cluster_size}')
    if idx == len(data_by_cluster.items())-1:
        axes[idx].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
        axes[idx].set_xticklabels(ligcols, rotation = 30, ha = 'right',
                                  fontsize = 12)
    else:
        axes[idx].set_xticks([])
plt.subplots_adjust(top=0.93,
                    bottom=0.05,
                    left=0.08,
                    right=0.84,
                    hspace=0.2,
                    wspace=0.2)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS12_full_clustermap/clustermap.png',
            bbox_inches='tight', dpi=300, transparent = True)
