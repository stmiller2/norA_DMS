#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:33:57 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import plot_utils as pl

#Load data
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='exclude')
data.columns = ['acriflavine', 'ethidium', 'norfloxacin', 'ofloxacin', 'pentamidine', 'pipemidic acid', 'puromycin', 'TPP']
ligcols = data.columns

#Get cluster assignments
cluster_assignments = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv', index_col = 'mut')['cluster']
data['cluster'] = cluster_assignments

#Set up minimum activity mask
active_mask = (data[ligcols] > -1).sum(axis = 1)>0 #Only variants with 1+ F>-1

#Get residue distances
dists = pd.read_csv(f'{base_dir}/results/general/residue_distances.csv',
                    index_col = 'mut')
data['bindingsite_dist'] = dists[['9b3m_bindingsite', '7lo8_bindingsite']].min(axis = 1) #minimum distance to the binding site in either conformation
data.loc[data['bindingsite_dist'].isna(), 'bindingsite_dist'] = dists.loc[data['bindingsite_dist'].isna(), 'alphafold_bindingsite'] #use alphafold for site not resolved in cryoem structures
data['coupling_dist'] = dists[['9b3m_coupling', '7lo8_coupling']].min(axis = 1) #minimum distance to the binding site in either conformation
data.loc[data['coupling_dist'].isna(), 'coupling_dist'] = dists.loc[data['coupling_dist'].isna(), 'alphafold_coupling'] #use alphafold for site not resolved in cryoem structures

#Get shell assigments
def assign_shell(d):
    if d == 0:    #directly binding-site exposed
        return "Shell 1"
    elif d <= 5:  #1–5 Å
        return "Shell 2"
    elif d <= 12: #5–12 Å
        return "Shell 3"
    else:
        return np.nan
data['shell'] = data['bindingsite_dist'].apply(assign_shell)

#Get SI scores
data_transformed = data[ligcols].replace(-3.0, np.nan).clip(upper = 0)
data['SI'] = data_transformed.max(axis = 1) - data_transformed.min(axis = 1)
data['SI'] = (data['SI'] - data['SI'].median())/data['SI'].std()
spec_cutoff = data['SI'].quantile(.95) #Look for mutations with SI scores above 95th percentile (~2.18)
spec = data[(data['SI'] >= spec_cutoff)]
spec = spec[active_mask.reindex(spec.index, fill_value=False)]

#%%
""" FIG S14A: SHELL 2 STACK PLOT """

spec_shell2 = spec[spec['shell'] == 'Shell 2']
fig = pl.stack_plot(spec_shell2, 6, 6)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS14_shell_highlights/shell2_stack.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S14D: SHELL 3 STACK PLOT """

spec_shell3 = spec[spec['shell'] == 'Shell 3']
fig = pl.stack_plot(spec_shell3, 6, 6)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS14_shell_highlights/shell3_stack.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S14B, C, E, F: SINGLE-RESIDUE HEATMAPS FOR SELECT SHELL HIGHLIGHTS """

for res in [220, 247, 35, 177]:
    fig = pl.residue_heatmap(res, cbar = False, return_fig = True)
    plt.show()
    
    fig.savefig(f'{base_dir}/results/supplementary/figS14_shell_highlights/{fig.axes[0].get_title()}.png', 
                bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" S14G-LL: SINGLE-RESIDUE HEATMAPS FOR SELECT CLUSTER HIGHLIGHTS """

for res in [244, 111, 310, 82, 15, 144]:
    fig = pl.residue_heatmap(res, cbar = False, return_fig = True)
    plt.show()
    
    fig.savefig(f'{base_dir}/results/supplementary/figS14_shell_highlights/{fig.axes[0].get_title()}.png', 
                bbox_inches = 'tight', dpi = 300, transparent = True)
    
