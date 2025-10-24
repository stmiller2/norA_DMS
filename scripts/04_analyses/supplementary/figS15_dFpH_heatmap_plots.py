#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 13:14:26 2024

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import matplotlib.pyplot as plt
import os
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import efficiency_utils as ef
import plot_utils as pl
from plot_utils import nor_effc_vmin, nor_effc_vmax, nor_effc_cmap

#Load data
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'mut_names')
effc.loc[effc.index.str[0] == effc.index.str[-1], 'nor_dFpH'] = 0 #Fill WT rows with 0
effc['position'] = effc.index.str[1:-1].astype(int)
effc['mutation'] = effc.index.str[-1]
effc.loc[effc['nor_dFpH'].isna(), 'nor_SE'] = np.nan #Don't show error bars if there's no score

spec = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'mut_names')
ligcols = spec.columns

#set up minimum activity mask
active_mask = (spec[ligcols] > -1).sum(axis = 1)>0 #Only variants with F>-1 in at least one substrate
active_mask = active_mask | (effc[['nor_pH60', 'nor_pH70']] > -1).sum(axis = 1) > 0 #OR F>-1 in at least one pH condition

#%%
""" FIG S15A: HEATMAP OF PH SENSITIVITY -- NORFLOXACIN """

fig = pl.dms_heatmap(effc.loc[active_mask], vmin = nor_effc_vmin, vmax = nor_effc_vmax, 
                      valcol = 'nor_dFpH', secol = 'nor_SE', 
                      title = 'ΔF$_{pH}$ of norfloxacin transport', 
                      cmap = nor_effc_cmap, error_limit = 0.5)
plt.show()

fig.savefig(f'{base_dir}/results/supplementary/figS15_dFpH_heatmap/dFpH_heatmap.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S15B: VOLCANO PLOT -- NORFLOXACIN """

#Get FDR-adjusted p values
effc['z'] = effc['nor_dFpH'] / effc['nor_SE']
effc['p'] = 2 * stats.norm.sf(abs(effc['z']))
mask = np.isfinite(effc['z'])
effc['fdr'] = np.nan
effc.loc[mask, 'fdr'] = multipletests(effc.loc[mask, 'p'], method = 'fdr_bh')[1]
effc['fdr'] = effc['fdr'].replace(0, effc[effc['fdr'] > 0]['fdr'].min()) #Fill values that round to zero with the next lowest p value

#Set plot parameters
sigma = effc['nor_dFpH'].std()
delta_cutoff = (-2*sigma, 2*sigma) #Significance cutoff: 2 standard deviations
sig = ((effc['fdr'] <= 0.05) & ((effc['nor_dFpH'] <= delta_cutoff[0]) | (effc['nor_dFpH'] >= delta_cutoff[1]))).tolist()
colors = ['#573280' if f < delta_cutoff[0] else '#F58300' if f > delta_cutoff[1] else 'tab:gray' for f in effc['nor_dFpH']]
colors = [colors[i] if sig[i] else 'tab:gray' for i in range(len(sig))]

#Plot
fig, ax = plt.subplots(figsize = (3, 3))
ax.scatter(effc['nor_dFpH'], -1*np.log10(effc['fdr']), s = 8, c = colors, alpha = 0.3)
ax.axhline(-np.log10(0.05), color = 'dimgray', linestyle = '--')
ax.axvline(delta_cutoff[0], color = 'dimgray', linestyle = '--')
ax.axvline(delta_cutoff[1], color = 'dimgray', linestyle = '--')
ax.set_ylabel('-log$_{10}$(FDR)')
ax.set_xlabel('F$_{pH7.0}$ − F$_{pH6.0}$', fontsize = 12)
ax.set_title('Norfloxacin', fontsize = 10)
lim = 1.1*abs(effc['nor_dFpH']).max()
ax.set_xlim(-lim, lim)
ax.set_xticks([-4, -2, 0, 2, 4])
ax.set_aspect(1/ax.get_data_ratio())
ax.text(delta_cutoff[1]*1.1, ax.get_ylim()[1]*.9, '2σ', fontsize = 9)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS15_dFpH_heatmap/volcano_nor.png', 
            dpi = 300, bbox_inches='tight', transparent = True)

#%%
""" FIG S15C: VOLCANO PLOT -- ACRIFLAVINE """

#Get FDR-adjusted p values
effc['z'] = effc['acr_dFpH'] / effc['acr_SE']
effc['p'] = 2 * stats.norm.sf(abs(effc['z']))
mask = np.isfinite(effc['z'])
effc['fdr'] = np.nan
effc.loc[mask, 'fdr'] = multipletests(effc.loc[mask, 'p'], method = 'fdr_bh')[1]
effc['fdr'] = effc['fdr'].replace(0, effc[effc['fdr'] > 0]['fdr'].min()) #Fill values that round to zero with the next lowest p value

#Set plot parameters
sigma = effc['acr_dFpH'].std()
delta_cutoff = (-2*sigma, 2*sigma) #Significance cutoff: 2 standard deviations
sig = ((effc['fdr'] <= 0.05) & ((effc['acr_dFpH'] <= delta_cutoff[0]) | (effc['acr_dFpH'] >= delta_cutoff[1]))).tolist()
colors = ['#573280' if f < delta_cutoff[0] else '#F58300' if f > delta_cutoff[1] else 'tab:gray' for f in effc['acr_dFpH']]
colors = [colors[i] if sig[i] else 'tab:gray' for i in range(len(sig))]

#Plot
fig, ax = plt.subplots(figsize = (3, 3))
ax.scatter(effc['acr_dFpH'], -1*np.log10(effc['fdr']), s = 8, c = colors, alpha = 0.3)
ax.axhline(-np.log10(0.05), color = 'dimgray', linestyle = '--')
ax.axvline(delta_cutoff[0], color = 'dimgray', linestyle = '--')
ax.axvline(delta_cutoff[1], color = 'dimgray', linestyle = '--')
ax.set_ylabel('-log$_{10}$(FDR)')
ax.set_xlabel('F$_{pH7.0}$ − F$_{pH6.0}$', fontsize = 12)
ax.set_title('Acriflavine', fontsize = 10)
lim = 1.1*abs(effc['acr_dFpH']).max()
ax.set_xlim(-lim, lim)
ax.set_xticks([-4, -2, 0, 2, 4])
ax.set_aspect(1/ax.get_data_ratio())
ax.text(delta_cutoff[1]*1.1, ax.get_ylim()[1]*.9, '2σ', fontsize = 9)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS15_dFpH_heatmap/volcano_acr.png', 
            dpi = 300, bbox_inches='tight', transparent = True)

#%%
""" FIG S15D: CORRELATION OF PROMISCUITY WITH EFFICIENCY AS MEASURED BY ACRIFLAVINE PH SENSITIVITY """

#Get data
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'exclude')
spec = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')
ligcols = spec.columns
spec['breadth'] = spec[ligcols].mean(axis = 1)

#Set up minimum activity mask
active_mask = (spec[ligcols] > -1).sum(axis = 1)>0 #Only variants with F>-1

#Plot
fig, ax = plt.subplots(figsize = (5,5), dpi = 300)
ax.set_title('variants with F > -1 for at least one substrate', fontsize = 10)
x = spec.loc[active_mask, 'breadth'] #Remove variants that may be misfolded (below activity cutoff) for functional analysis
y = effc.loc[active_mask, 'acr_dFpH'] #Remove variants that may be misfolded (below activity cutoff) for functional analysis
ax.scatter(x, y, alpha = 0.3, color = 'tab:gray')
ax.set_xlabel('Promiscuity (F$_{avg}$)', fontsize = 12)
ax.set_ylabel('Efficiency (F$_{pH7.0}$ − F$_{pH6.0}$)', fontsize = 12)

# Histogram on the x-axis
ax_histx = ax.inset_axes([0, 1, 1, 0.1])
ax_histx.hist(x, bins=100, color='lightgray', edgecolor = 'black', linewidth = .5, histtype = 'stepfilled')
ax_histx.set_xlim(ax.get_xlim())
ax_histx.set_xticks([])
ax_histx.set_yticks([])

# Histogram on the y-axis
ax_histy = ax.inset_axes([1, 0, 0.1, 1])
ax_histy.hist(y, bins=100, orientation='horizontal', color='lightgray', edgecolor = 'black', linewidth = .5, histtype = 'stepfilled')
ax_histy.set_ylim(ax.get_ylim())
ax_histy.set_xticks([])
ax_histy.set_yticks([])

# Add zeros markers
plt.axvline(0, color = 'gray', linestyle = 'dashed')
plt.axhline(0, color = 'gray', linestyle = 'dashed')

# Display Spearman correlation
valid_mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isinf(x) & ~np.isinf(y)
rho, pval = stats.spearmanr(x[valid_mask], y[valid_mask])
text_str = f"Spearman R = {rho:.2f}\np < $10^{-16}$" #pval rounds to zero
bbox_props = dict(boxstyle="round,pad=0.5", edgecolor="black", linewidth = .75, facecolor=(1,1,1,.2))
ax.text(0.43, 0.92, text_str, transform=ax.transAxes, fontsize=10, verticalalignment='top', ha = 'right', bbox=bbox_props)

for spine in ['left', 'right', 'top', 'bottom']:
    ax_histx.spines[spine].set_visible(False)
    ax_histy.spines[spine].set_visible(False)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS15_dFpH_heatmap/acr_scatter.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)    
    