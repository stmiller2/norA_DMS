#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:24:04 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import plot_utils as pl


""" Load data """
data = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'WT')
data.columns = ['acriflavine', 'ethidium', 'norfloxacin', 'ofloxacin', 'pentamidine', 'pipemidic acid', 'puromycin', 'TPP']
ligcols = data.columns

#get cluster assignments
cluster_assignments = pd.read_csv('/Volumes/sraman4/General/Silas/norA_DMS/writing/revision/figures/clustering/cluster_assignments.csv',
                                  index_col = 'mut')
data['cluster'] = cluster_assignments['cluster']

#get custom colormap & color palletes
from plot_utils import sub_colors, spec_cmap, spec_vmin, spec_vmax

#%%
""" FIG 1B COMPONENT: SINGLE-CLUSTER INSET HEATMAP """
fig, ax = plt.subplots(figsize = (3,.75))
sns.heatmap(data[data['cluster'] == 37][ligcols],
            cmap = spec_cmap,
            center = 0,
            vmin = spec_vmin, vmax = spec_vmax,
            cbar = False)
ax.set_xticks([])
ax.set_yticks([])
for x in np.arange(1, 8, 1):
    ax.axvline(x, color = 'black', linewidth = 0.1)
ax.set_ylabel('')
ax.set_xlabel('')
for spine in ['top', 'bottom', 'left', 'right']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1)
plt.show()
fig.savefig(f'{base_dir}/results/fig1/charge_cluster_inset.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 1D: HISTOGRAM COMPARISON """
fig, axes = plt.subplots(nrows = 8, figsize = (3,4.5), sharey = True)
for i, (lig, ax) in enumerate(zip(ligcols, axes)):
    cleaned = data[data.index != 'WT'][lig].replace(-3, np.nan).dropna()
    ax.hist(cleaned, bins = np.arange(-3, 1, .05), color=sub_colors[i], label = lig,
            edgecolor = 'black', histtype = 'stepfilled', linewidth = .75)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(-3.25, 1.25)
    ax.text(0, 0.075, lig if lig == 'TPP' else lig.capitalize(), 
            ha = 'left', transform = ax.transAxes,
            bbox = dict(facecolor = 'white', edgecolor = 'white', alpha = 0.5, pad=0.05))
    ax.patch.set_alpha(0)
    for spine in ['left','top','right']:
        ax.spines[spine].set_visible(False)
    if i == len(ligcols)//2:
        ax.set_ylabel('Frequency')
        ax.yaxis.set_label_coords(-0.05, .4)
    if i == len(ligcols)-1:
        ax.set_xticks([-3, -2, -1, 0, 1])
        ax.set_xlabel('Standardized functional score')
plt.subplots_adjust(hspace = -0.6)
plt.show()
fig.savefig(f'{base_dir}/results/fig1/histograms.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 1E: CORRELATION HEATMAP """
correlation_matrix = data[ligcols].corr(method = 'spearman') 
#correlation_matrix = data[~data['cluster'].isin([7 , 26])][ligcols].corr(method = 'spearman') #Remove univ. permitted & disabling mutations (focus only on specificity-driving muts)

fig = sns.clustermap(correlation_matrix,
                     cmap = sns.light_palette("#E44B3A", as_cmap=True),
                     vmin = correlation_matrix.min().min(), vmax = 1,
                     annot = True,
                     figsize = (4,4),
                     dendrogram_ratio = .1,
                     cbar_pos = None,
                     tree_kws = {'linewidth': 1},
                     cbar_kws={'orientation': 'horizontal', 'label': 'Pearson correlation'},
                     annot_kws={'size':8})
labs = [label.get_text().capitalize() for label in fig.ax_heatmap.get_xticklabels()]
labs = ["TPP" if lab == "Tpp"
        else "Pipemidic\nacid" if lab == "Pipemidic acid"
        else lab for lab in labs]
plt.xticks(np.arange(0.5,8.5,1), labs, ha = 'right', rotation_mode = 'anchor', fontsize = 8)
plt.yticks(np.arange(0.5,8.5,1), labs, ha = 'left', fontsize = 8)
hm = fig.ax_heatmap
hm.set_ylabel(None)
hm.set_aspect(1./hm.get_data_ratio())
plt.setp(hm.xaxis.get_majorticklabels(), rotation = 45, ha = 'right')
plt.show()
fig.savefig(f'{base_dir}/results/fig1/correlation.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

