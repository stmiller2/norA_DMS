#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 13:14:26 2024

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from joypy import joyplot
import matplotlib as mpl
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import efficiency_utils as ef

""" Load data """
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'exclude')
effc['position'] = effc.index.str[1:-1].astype(int)
effc['mutation'] = effc.index.str[-1]

spec = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')
ligcols = spec.columns

#set up minimum activity mask
active_mask = (spec[ligcols] > -1).sum(axis = 1)>0 #Only variants with 1+ F>-1

spec['breadth'] = spec[ligcols].mean(axis = 1)
spec['num_active'] = (spec[ligcols] > -1).sum(axis = 1)

#get custom colormap
from plot_utils import nor_effc_cmap, nor_effc_vmin, nor_effc_vmax

#%%
""" FIG 4A: LINES SHOWING CHANGE ACROSS PH CONDITIONS -- NORFLOXACIN """

a = effc['nor_pH60']
b = effc['nor_pH70']
slopes = effc['nor_dFpH']

limit = max(abs(nor_effc_vmin*1.1), abs(nor_effc_vmax*1.1))
norm = mpl.colors.Normalize(-limit, limit)
sorted_indices = slopes.sort_values(key = abs).index
sigma = slopes.std()

# Draw lines between corresponding points in Acid and Base
fig, ax = plt.subplots(figsize=(2, 8))
s = 8; color = 'gray'
ax.scatter(['pH 6.0'] * len(a), a, marker='o', s = s, color = color, edgecolor = 'black', linewidth = 0.5, zorder = 3)
ax.scatter(['pH 7.0'] * len(b), b, marker='o', s = s, color = color, edgecolor = 'black', linewidth = 0.5, zorder = 3)
for i in sorted_indices:
    color = nor_effc_cmap(norm(slopes[i]))
    linewidth = 2
    ax.plot(['pH 6.0', 'pH 7.0'], [a[i], b[i]], color=color, linewidth=linewidth, zorder = 1)
    
#Show WT
ax.plot(['pH 6.0', 'pH 7.0'], [0, 0], color = 'gray', linewidth = 2, zorder = 1)
ax.text(0.75,0.05,'WT', fontsize = 8, weight = 'bold')

#Show select variants
highlight_color = '#090909'
bbox_props = dict(boxstyle="round,pad=0.2", edgecolor="black", linewidth = .5, facecolor=(1,1,1,.5))
ax.plot(['pH 6.0', 'pH 7.0'], [a['E222D'], b['E222D']], color = highlight_color, linewidth = 1, zorder = 1)
ax.text(0.15,-1.6,'E222D', fontsize = 7, weight = 'bold', bbox = bbox_props)
ax.plot(['pH 6.0', 'pH 7.0'], [a['D307E'], b['D307E']], color = highlight_color, linewidth = 1, zorder = 1)
ax.text(0.6,-1.2,'D307E', fontsize = 7, weight = 'bold', bbox = bbox_props)
ax.plot(['pH 6.0', 'pH 7.0'], [a['R98K'], b['R98K']], color = highlight_color, linewidth = 1, zorder = 1)
ax.text(0.5,-2.4,'R98K', fontsize = 7, weight = 'bold', bbox = bbox_props)

#Set labels
ax.set_ylabel('Functional score – Norfloxacin')
ax.margins(x=0.25)

plt.show()
fig.savefig(f'{base_dir}/results/fig4/slopes.png', 
            dpi = 300, bbox_inches='tight', transparent = True)

#%%
""" FIG 4A: HISTOGRAM / COLORBAR -- NORFLOXACIN """

fig, ax = plt.subplots(figsize = (.75,3))

#Calculate histogram data
bins = np.linspace(min(slopes), max(slopes), 40)
hist, bin_edges = np.histogram(slopes, bins=bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
bin_width = bin_edges[1] - bin_edges[0]

#Plot histogram
ax.barh(bin_centers, hist, height=bin_width, color=nor_effc_cmap(norm(bin_centers)), edgecolor='black', linewidth = 0.8)
xlims = ax.get_xlim()
ax.set_xlim(xlims[0], xlims[1]*1.5)

#Plot colorbar
cbar = mpl.colorbar.ColorbarBase(ax.inset_axes([-.2, 0, .2, 1]), cmap=nor_effc_cmap, norm=norm, orientation='vertical')
ax.set_title('ΔF$_{pH}$', y = 1.07, x = -0.15)
ax.text(-0.1, 1.03, '(F$_{pH7.0}$ − F$_{pH6.0}$)', transform=ax.transAxes, ha = 'center', fontsize = 9)
ymin, ymax = ax.get_ylim()
cbar.ax.set_ylim(ymin, ymax-0.06)
cbar.ax.yaxis.set_ticks_position('left')  # Move ticks to the left
cbar.outline.set_linewidth(0.5)

#Add sigma lines and label
ax.axhline(-1*sigma, linestyle='--', color='gray', linewidth=0.8)
ax.axhline(sigma, linestyle='--', color='gray', linewidth=0.8, label='σ')  # Adjusting for horizontal lines
ax.text(.8*hist.max(), sigma*1.1, s =' σ', fontsize = 7)
ax.text(.8*hist.max(), sigma*-1.25, s ='-σ', fontsize = 7)

#Adjust style
for spine in ['top', 'bottom', 'right']:
    ax.spines[spine].set_visible(False) # Remove border
ax.set_yticks([])  # Remove x-axis ticks
ax.set_xticks([]) # Remove y-axis ticks

plt.show()
fig.savefig(f'{base_dir}/results/fig4/hist_nor.png', 
            dpi =300, bbox_inches='tight', transparent = True)

#%%
""" FIG 4C: KDE POTS COMPARING SPECIFICITY CLUSTERS """

#Get cluster assignments
clust = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv',
                    index_col = 'mut')['cluster']
clust = clust[clust.index != 'WT']

clusters = {31: 'Puromycin-specific\ntolerance',
            12: 'TPP-specific\ntolerance',
            32: 'High TPSA-specific\ntolerance',
            37: 'Cation-specific\nimpairment',
            25: 'Pentamidine-specific\nimpairment',
            2: 'Ethidium-specific\ntolerance',
            14: 'Universally\nenriched'}

colors = ['#4f8da7', '#86B26E', '#e1b23a', '#ee842b', '#e45c3a', '#C66C79', '#A87BB7', '#DADCDF']
plt.rc("font", size=12) 
df_list = []

for cluster, label in clusters.items():
    subset = effc.loc[active_mask & (clust == cluster), 'nor_dFpH'] #Remove variants that may be misfolded (below activity cutoff) for functional analysis
    df_list.append(pd.DataFrame({'nor_dFpH': subset, 'Cluster': label}))
df_list.append(pd.DataFrame({'nor_dFpH': effc['nor_dFpH'], 'Cluster': 'All Data'}))
plot_data = pd.concat(df_list)
fig, axes = joyplot(
            data=[plot_data[plot_data['Cluster'] == label]['nor_dFpH'] for label in clusters.values()] + [plot_data[plot_data['Cluster'] == 'All Data']['nor_dFpH']],
            labels=list(clusters.values()) + ['All Data'],
            kind="kde",
            overlap=1,
            fade=False,
            figsize = (6, 6),
            color = colors,
            x_range = (-3.9, 2.5)
            )

plt.axvline(0, ymax = .9, color = 'gray', linestyle = '--')
plt.xlabel('$ΔF_{pH}$')
plt.show()

fig.savefig(f'{base_dir}/results/fig4/spec_kde.png',
            dpi = 300, bbox_inches = 'tight', transparent = True)

#%%
""" FIG 4B: EFFICIENCY SCORE VS. DISTANCE FROM COUPLING RESIDUES """

#Get coupling distances
dists = pd.read_csv(f'{base_dir}/results/general/residue_distances.csv', index_col = 'mut')
spec['dist_to_coupling'] = dists[['9b3m_coupling', '7lo8_coupling']].min(axis = 1) #minimum distance to the binding site in either conformation
spec.loc[spec['dist_to_coupling'].isna(), 'dist_to_coupling'] = dists.loc[spec['dist_to_coupling'].isna(), 'alphafold_coupling'] #use alphafold for site not resolved in cryoem structures

#Custom colormap
colors = ['#381462', '#573280', '#C66C79',  '#ee842b',  '#FAD16B']
cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_plasma", colors, N=256)

#Mask NaN values
mask = ~(effc['nor_dFpH'].isna() | spec['dist_to_coupling'].isna())

#Get density
xy = np.vstack([effc.loc[mask, 'nor_dFpH'], spec.loc[mask, 'dist_to_coupling']])
density = stats.gaussian_kde(xy)(xy)
idx = density.argsort()
x = effc.loc[mask, 'nor_dFpH'].iloc[idx]
y = spec.loc[mask, 'dist_to_coupling'].iloc[idx]
density = density[idx]

#Plot
fig, ax = plt.subplots(figsize = (4.5,4))
sc = ax.scatter(x, y, c=density, cmap=cmap, alpha = 0.7, s = 15, edgecolors=None)
cbar = plt.colorbar(sc, ax=ax, fraction = 0.025, pad = 0.04)
cbar.set_label('Density')
cbar.ax.set_yticks([])

ax.set_title('Mutations near coupling residues are more\nlikely to impair efficiency', fontsize = 12)
ax.set_xlabel('Energy efficiency (ΔF$_{pH})$')
ax.set_ylabel('Distance to nearest coupling residue (Å)')
plt.show()

fig.savefig(f'{base_dir}/results/fig4/distance_scatter.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 4D: SPECIFICITY-EFFICIENCY CORRELATION SCATTER """

#Plot
fig, ax = plt.subplots(figsize = (5,5), dpi = 300)
ax.set_title('variants with F > -1 for at least one substrate', fontsize = 10)
x = spec.loc[active_mask, 'breadth'] #Remove variants that may be misfolded (below activity cutoff) for functional analysis
y = effc.loc[active_mask, 'nor_dFpH'] #Remove variants that may be misfolded (below activity cutoff) for functional analysis
ax.scatter(x, y, alpha = 0.3, color = 'tab:gray')
ax.set_xlabel('Promiscuity (F$_{avg}$)', fontsize = 12)
ax.set_ylabel('Efficiency (F$_{pH7.0}$ − F$_{pH6.0}$)', fontsize = 12)

#Histogram on the x-axis
ax_histx = ax.inset_axes([0, 1, 1, 0.1])
ax_histx.hist(x, bins=100, color='lightgray', edgecolor = 'black', linewidth = .5, histtype = 'stepfilled')
ax_histx.set_xlim(ax.get_xlim())
ax_histx.set_xticks([])
ax_histx.set_yticks([])

#Histogram on the y-axis
ax_histy = ax.inset_axes([1, 0, 0.1, 1])
ax_histy.hist(y, bins=100, orientation='horizontal', color='lightgray', edgecolor = 'black', linewidth = .5, histtype = 'stepfilled')
ax_histy.set_ylim(ax.get_ylim())
ax_histy.set_xticks([])
ax_histy.set_yticks([])

#Add zeros markers
plt.axvline(0, color = 'gray', linestyle = 'dashed')
plt.axhline(0, color = 'gray', linestyle = 'dashed')

#Display Spearman correlation
valid_mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isinf(x) & ~np.isinf(y)
rho, pval = stats.spearmanr(x[valid_mask], y[valid_mask])
text_str = f"Spearman R = {rho:.2f}\np < $10^{-16}$"
bbox_props = dict(boxstyle="round,pad=0.5", edgecolor="black", linewidth = .75, facecolor=(1,1,1,.2))
ax.text(0.43, 0.92, text_str, transform=ax.transAxes, fontsize=10, verticalalignment='top', ha = 'right', bbox=bbox_props)

for spine in ['left', 'right', 'top', 'bottom']:
    ax_histx.spines[spine].set_visible(False)
    ax_histy.spines[spine].set_visible(False)

plt.show()
fig.savefig(f'{base_dir}/results/fig4/scatter.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)    


#%%
""" FIG 4E: VIOLIN PLOTS BY NUM. TRANSPORTED SUBSTRATES """

#Mask NaN values
valid_mask = effc['nor_dFpH'].notna()
rnge = [1,2,3,4,5,6,7,8]
grouped_data = [effc[valid_mask].loc[spec['num_active'] == i, 'nor_dFpH'] for i in rnge]

#Plot
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
for i, data in enumerate(grouped_data):
    jitter = np.random.normal(0, 0.07, size=len(data))
    ax.scatter(np.full(len(data), i+1) + jitter, data, c='tab:gray', alpha = 0.5,
               s=10, edgecolor='none')
violin_parts = ax.violinplot(grouped_data, positions=rnge, widths=0.5, showmeans=False, showmedians=True, showextrema=False)
for pc in violin_parts['bodies']:
    pc.set_edgecolor('black')
    pc.set_facecolor('white')
    pc.set_alpha(.75)
violin_parts['cmedians'].set_color('black')

#Label and annotate
ax.set_xlabel('Promiscuity (Num. substrates with $F$ > -1)', fontsize = 12)
ax.set_ylabel('Efficiency (F$_{pH7.0}$ − F$_{pH6.0}$)', fontsize = 12)
ax.text(0, ax.get_ylim()[1], 'N =', ha='center', va='bottom', fontsize=10)
for i, data in enumerate(grouped_data):
    ax.text(i+1, ax.get_ylim()[1], f'{len(data)}', ha='center', va='bottom', fontsize=10)    
ax.set_xticks([1,2,3,4,5,6,7,8], [1,2,3,4,5,6,7,8])

plt.show()
fig.savefig(f'{base_dir}/results/fig4/violin.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

