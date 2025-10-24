#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:54:15 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import matplotlib.pyplot as plt
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
                                wt_policy = 'exclude')
data.columns = ['acriflavine', 'ethidium', 'norfloxacin', 'ofloxacin', 'pentamidine', 'pipemidic acid', 'puromycin', 'TPP']
ligcols = data.columns

#get cluster assignments
cluster_assignments = pd.read_csv('/Volumes/sraman4/General/Silas/norA_DMS/writing/revision/figures/clustering/cluster_assignments.csv',
                                  index_col = 'mut')
data['cluster'] = cluster_assignments['cluster']

#set up minimum activity mask
active_mask = (data[ligcols] > -1).sum(axis = 1)>0 #Only variants with 1+ F>-1

#get custom colormap
from plot_utils import aa_color_dict, spec_cmap, spec_vmin, spec_vmax

def print_pymol_commands(data, cluster, aa_color_dict):
    cluster_data = data[data['cluster'] == cluster]
    positions = '+'.join(cluster_data.index.str[1:-1].unique().tolist())
    print(f'sele resi {positions}')
    print(f'set_name sele, cluster{cluster}_positions')
    seen = set()
    for mut in cluster_data.index:
        wt = mut[0]
        pos = ''.join([c for c in mut if c.isdigit()])
        if pos in seen:
            continue
        seen.add(pos)
        color = '0x'+aa_color_dict.get(wt, "white")[1:]
        print(f"show sticks, resi {pos}")
        print(f"color {color}, resi {pos}")


#%%
""" Cluster 37 - charge """
c = 37
clusterdata = data[data['cluster'] == c].loc[active_mask.reindex(data[data['cluster'] == c].index, fill_value=False), ligcols] #remove variants which may be misfolded (below activity cutoff)
hm = pl.cluster_heatmap_small(clusterdata, cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax)
plt.show()
sp = pl.stack_plot(clusterdata, 2, 6, text_x = -0.21)
plt.show()
hm.savefig(f'{base_dir}/results/fig3/charge_hm.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
sp.savefig(f'{base_dir}/results/fig3/charge_stack.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
print_pymol_commands(data, c, aa_color_dict) #Use 7LO8

#%%
""" Cluster 31 - puromycin tolerant """
c = 31
clusterdata = data[data['cluster'] == c].loc[active_mask.reindex(data[data['cluster'] == c].index, fill_value=False), ligcols] #remove variants which may be misfolded (below activity cutoff)
hm = pl.cluster_heatmap_small(clusterdata, cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax)
plt.show()
sp = pl.stack_plot(clusterdata, 5, 4, text_x = -0.21)
plt.show()
hm.savefig(f'{base_dir}/results/fig3/pur_hm.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
sp.savefig(f'{base_dir}/results/fig3/pur_stack.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
print_pymol_commands(data, c, aa_color_dict) #Use 9B3M

#%%
""" Cluster 32 - TPSA """
c = 32
clusterdata = data[data['cluster'] == c].loc[active_mask.reindex(data[data['cluster'] == c].index, fill_value=False), ligcols] #remove variants which may be misfolded (below activity cutoff)
hm = pl.cluster_heatmap_small(clusterdata, cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax)
plt.show()
sp = pl.stack_plot(clusterdata, 2, 5, text_x = -0.21)
plt.show()
hm.savefig(f'{base_dir}/results/fig3/tpsa_hm.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
sp.savefig(f'{base_dir}/results/fig3/tpsa_stack.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
print_pymol_commands(data, c, aa_color_dict) #Use 9B3M

#%%
""" Cluster 12 - TPP """
c = 12
clusterdata = data[data['cluster'] == c].loc[active_mask.reindex(data[data['cluster'] == c].index, fill_value=False), ligcols] #remove variants which may be misfolded (below activity cutoff)
hm = pl.cluster_heatmap_small(clusterdata, cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax)
plt.show()
sp = pl.stack_plot(clusterdata, 4, 6)
plt.show()
hm.savefig(f'{base_dir}/results/fig3/tpp_hm.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
sp.savefig(f'{base_dir}/results/fig3/tpp_stack.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
print_pymol_commands(data, 12, aa_color_dict) #Use 7LO8

#%%
""" Cluster 25 - pentamidine """
c = 25
clusterdata = data[data['cluster'] == c].loc[active_mask.reindex(data[data['cluster'] == c].index, fill_value=False), ligcols] #remove variants which may be misfolded (below activity cutoff)
hm = pl.cluster_heatmap_small(clusterdata, cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax)
plt.show()
sp = pl.stack_plot(clusterdata, 2, 3)
plt.show()
hm.savefig(f'{base_dir}/results/fig3/pen_hm.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
sp.savefig(f'{base_dir}/results/fig3/pen_stack.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
print_pymol_commands(data, 25, aa_color_dict) #Use 9B3M

#%%
""" Cluster 2 - ethidium """
c = 2
clusterdata = data[data['cluster'] == c].loc[active_mask.reindex(data[data['cluster'] == c].index, fill_value=False), ligcols] #remove variants which may be misfolded (below activity cutoff)
hm = pl.cluster_heatmap_small(clusterdata, cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax)
plt.show()
sp = pl.stack_plot(clusterdata, 5, 2)
plt.show()
hm.savefig(f'{base_dir}/results/fig3/eth_hm.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
sp.savefig(f'{base_dir}/results/fig3/eth_stack.png',
           bbox_inches = 'tight', dpi = 300, transparent = True)
print_pymol_commands(data, 2, aa_color_dict) #Use 9B3M
