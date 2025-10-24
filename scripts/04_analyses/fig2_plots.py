#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 09:46:39 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import scipy.cluster.hierarchy as sch
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
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
                                wt_policy = 'exclude')
data.columns = ['acriflavine', 'ethidium', 'norfloxacin', 'ofloxacin', 'pentamidine', 'pipemidic acid', 'puromycin', 'TPP']
ligcols = data.columns

#get cluster assignments
cluster_assignments = pd.read_csv('/Volumes/sraman4/General/Silas/norA_DMS/writing/revision/figures/clustering/cluster_assignments.csv',
                                  index_col = 'mut')
data['cluster'] = cluster_assignments['cluster']

#set up minimum activity mask (here, cluster 7 is manually included so it can be plotted in Fig. 2A)
active_mask = (data[ligcols] > -1).sum(axis = 1)>0 #Only variants with 1+ F>-1
active_mask.loc[data['cluster'] == 7] = True #Keep cluster 7 for plotting in split heatmap

#get custom colormap
from plot_utils import spec_cmap, spec_vmin, spec_vmax

#%%
""" FIG 2A: FULL SPLIT CLUSTERMAP """
min_cluster_size = 3
pctile_cutoff = 25

#Get cluster rank order
clusters = pd.DataFrame()
clusters.index = data['cluster'].value_counts()[data['cluster'].value_counts() > 1].index
clusters['size'] = data['cluster'].value_counts()[data['cluster'].value_counts() > 1].values
dists = []

for cluster in clusters.index:
    linkage_matrix = sch.linkage(data[data['cluster'] == cluster][ligcols].T,
                                                     metric = 'euclidean', method = 'average')
    # This is the euclidean distance between the final two clusters of ligands to be merged,
    # normalized to the number of datapoints in the dataset. It is a measure of how well the
    # mutations in this cluster separate different ligands - i.e., the magnitude of 
    # specificity-driving-ness for the given cluster.
    final_merge_distance_normalized = linkage_matrix[-1,2]/len(data[data['cluster'] == cluster])
    dists.append(final_merge_distance_normalized)
    
clusters['fmd'] = dists


clusters = clusters.sort_values('fmd')
spec_clusters = clusters[clusters['fmd'] > np.percentile(clusters['fmd'], pctile_cutoff)]
spec_clusters = spec_clusters[spec_clusters['size'] >= min_cluster_size]
data_by_cluster = {cluster: data[data['cluster'] == cluster] for cluster in [26, 7, 14] + spec_clusters.index.tolist()}
# --- ^ Biologically relevant clusters to display at the top of the plot regardless of fmd:
# ---   26: Universally permitted; 7: Universally disabling; 14: Universally enriched

n_clusters = sum([len(c[active_mask.reindex(c.index, fill_value=False)]) >= min_cluster_size
                  for c in data_by_cluster.values()])  # number of clusters meeting min cluster size after filtering for activity

mosaic = []
height_ratios = []

for i in range(n_clusters):
    mosaic.append([str(i)])
    
    #heights of the cluster rows
    if i in (0, 1):
        height_ratios.append(2)   #first two clusters tall
    else:
        height_ratios.append(1)   #rest uniform
        
    if i == 0:  #separator after first cluster
        mosaic.append(["."])
        height_ratios.append(0.75)
    elif i == 1 and n_clusters > 2:  #separator after second cluster
        mosaic.append(["."])
        height_ratios.append(0.75)

fig = plt.figure(figsize=(3, 8.5))
axes = fig.subplot_mosaic(mosaic, height_ratios=height_ratios)

cbar_ax = fig.add_axes([-0.03, 0.05, .05, .18])

plot_idx = 0
for cluster, data_cluster in data_by_cluster.items():
    cluster_active_mask = active_mask.reindex(data_cluster.index) #reindex activity mask for variants in this cluster
    data_cluster = data_cluster[cluster_active_mask] #remove variants which may be misfolded (below activity cutoff)

    if len(data_cluster) < min_cluster_size:
        continue

    ax = axes[str(plot_idx)]
    sns.heatmap(
        data_cluster.iloc[:, :8], cmap=spec_cmap, center=0,
        ax=ax, vmin=spec_vmin, vmax=spec_vmax,
        cbar=(plot_idx == 0), cbar_ax=None if plot_idx else cbar_ax
    )
    ax.set_yticks([])
    ax.set_ylabel(len(data_cluster), rotation=0, fontsize=8, ha='right')
    ax.yaxis.set_label_coords(-0.02, 0.3)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(.5)
    for x in np.arange(1, 8, 1):
        ax.axvline(x, color='black', linewidth=0.1)
    if plot_idx == n_clusters - 1:  # careful: n_clusters should now mean "plotted clusters"
        ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
        ax.set_xticklabels(ligcols, rotation=30, ha='right', rotation_mode = 'anchor', fontsize=10)
    else:
        ax.set_xticks([])

    plot_idx += 1

cbar_ax.yaxis.tick_left()
cbar_ax.yaxis.set_label_position('left')
cbar_ax.set_yticks([2, 1, 0, -1, -2, -3])
cbar_ax.set_yticklabels([2, 1, 0, -1, -2, -3], fontsize = 11)
cbar_ax.yaxis.set_label_coords(-1.8, .5)
cbar_ax.set_ylabel('Standardized\nfunctional score', fontsize = 11)

# Titles
axes['0'].set_title('Universally permitted mutations', weight='bold', fontsize=10, y=0.9)
axes['1'].set_title('Universally disabling mutations', weight='bold', fontsize=10, y=0.9)
axes['2'].set_title('Specificity-driving mutations', weight='bold', fontsize=10, y=.8)

plt.subplots_adjust(top=0.93, bottom=0.05, left=0.08, right=0.84, hspace=0, wspace=0.2)
fig.text(0.065, .915, 'Cluster size', ha='right', weight='bold', fontsize=9)
plt.show()

fig.savefig(f'{base_dir}/results/fig2/split_clustermap.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 2E: POSITIONAL DEPENDENCE LINE PLOT """
data['position'] = data.index.str[1:-1].astype(int)
neu = data[data['cluster'] == 26]; imp = data[data['cluster'] == 7]
window = 5
helices = [(3,33),(37,65),(67,84),(92,119),(125,149),(157,176),(203,232),(240,264),(268,288),(292,320),(325,353),(357,382)]
colors = ['#95110F','#2F6035','#365659','#EB3633','#4A9653','#55878B','#EE5A58','#69B572','#81AEB1','#F49190','#92C998','#9ABEC1']

fig = plt.figure(figsize = (8, 1.5))
axes = fig.subplot_mosaic([['a'],['b']],
                          height_ratios = [10, 1.3])
neu_hist = np.histogram(neu['position'], bins = 388)[0]
neu_ma = np.convolve(neu_hist/21, np.ones(window)/window, mode = 'same')
imp_hist = np.histogram(imp['position'], bins = 388)[0]
imp_ma = np.convolve(imp_hist/21, np.ones(window)/window, mode = 'same')
axes['a'].plot(np.arange(1, 389, 1), neu_ma, color='tab:gray', linewidth=1, label = 'Universally permitted')
axes['a'].fill_between(np.arange(1, 389, 1), neu_ma, color = '#d7dbe2', alpha = 0.4)
axes['a'].plot(np.arange(1, 389, 1), imp_ma, color='#1E86B3', linewidth=1, label = 'Universally disabling')
axes['a'].fill_between(np.arange(1, 389, 1), imp_ma, color = '#1E86B3', alpha = 0.4)
axes['a'].set_xticks([])
axes['a'].set_yticks([0, .5, 1])
axes['a'].set_ylabel('Mutation frequency')
axes['a'].spines['top'].set_visible(False)
axes['a'].spines['right'].set_visible(False)
axes['a'].legend(bbox_to_anchor=(0.03, 1.4), loc = 'upper left')
#TM domain graphic
axes['b'].plot(np.arange(1, 390, 1), 0.5*np.ones(389), color = 'black', linewidth = 2)
for i, helix in enumerate(helices):
    start, end = helix
    axes['b'].fill_between(np.arange(1, 389, 1)[start:end], 0, 1, color = colors[i], zorder = 3)
    axes['b'].text((start + (end - start) / 2)+0.5, .4, f'{i + 1}', ha='center', va='center', color = 'black',  fontsize = 8)
axes['b'].text(-10, 0, 'N', weight = 'bold', fontsize = 12)
axes['b'].text(393, 0, 'C', weight = 'bold', fontsize = 12)
axes['b'].set_yticks([])  
axes['b'].set_xticks([])  
axes['b'].set_xlabel('Amino acid position')
axes['b'].spines['top'].set_visible(False)
axes['b'].spines['right'].set_visible(False)
axes['b'].spines['bottom'].set_visible(False)
axes['b'].spines['left'].set_visible(False)
plt.subplots_adjust(hspace = 0.1, wspace = 0.1)
plt.show()
fig.savefig(f'{base_dir}/results/fig2/position_frequency.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 2D: IMPAIRED MUTATION FREQUENCY ON STRUCTURE - PYMOL SELECTION LISTS AND COLORBAR """
imp_counts = data[data['cluster'] == 7]['position'].value_counts()
select_residues = lambda l, u: '+'.join(imp_counts[(l <= imp_counts) & (imp_counts <= u)].index.astype(str))
print('color white')
for color, rng in zip(['0x6EA5C7', '0x215278', '0x00273D'], [(5, 9), (10, 14), (15, 20)]):
    print(f'\ncolor {color}, resi {select_residues(*rng)}')

#create colorbar
cmap_struct = mpl.colors.ListedColormap(['white', '#6EA5C7', '#215278', '#00273D'])
norm = mpl.colors.BoundaryNorm([0, 5, 10, 15, 21], cmap_struct.N)
fig, ax = plt.subplots(figsize=(0.2, 3))
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap_struct, norm=norm), 
                    cax=ax, orientation='vertical', spacing='proportional', 
                    extend='neither')
cbar.set_ticks([2.5, 7.5, 12.5, 18])
cbar.set_ticklabels(['0-4', '5-9', '10-14', '15-20'])
cbar.set_label('Num. disabling mutations', rotation=270, labelpad=15)
cbar.outline.set_linewidth(.5)
plt.show()

fig.savefig(f'{base_dir}/results/fig2/struct_cbar.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 2B: PERMITTED / DISABLING STACK PLOTS """
#Universally permitted cluster
fig = pl.stack_plot(data[data['cluster'] == 26], 5, 5, text_x = -0.21)
plt.show()
fig.savefig(f'{base_dir}/results/fig2/permitted_stack.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#Universally disabling cluster
fig = pl.stack_plot(data[data['cluster'] == 7], 6, 5, text_x = -0.21)
plt.show()
fig.savefig(f'{base_dir}/results/fig2/disabling_stack.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 2G, H: PROXIMAL / DISTAL SPECIFICITY DRIVER STACK PLOTS """
#Get binding site distances
dists = pd.read_csv(f'{base_dir}/results/general/residue_distances.csv', index_col = 'mut')
data['bindingsite_dist'] = dists[['9b3m_bindingsite', '7lo8_bindingsite']].min(axis = 1) #minimum distance to the binding site in either conformation
data.loc[data['bindingsite_dist'].isna(), 'bindingsite_dist'] = dists.loc[data['bindingsite_dist'].isna(), 'alphafold_bindingsite'] #use alphafold for residues not resolved in cryoEM structures

#Compute SI scores
data_transformed = data[ligcols].replace(-3.0, np.nan).clip(upper = 0)
data['SI'] = data_transformed.max(axis = 1) - data_transformed.min(axis = 1)
data['SI'] = (data['SI'] - data['SI'].median())/data['SI'].std()
spec_cutoff = data['SI'].quantile(.75) + (1.5*(data['SI'].quantile(.75) - data['SI'].quantile(.25))) #"Specificity-driving" means high outilers (above upper Tukey fence, ~2.77)
spec = data[data['SI'] >= spec_cutoff]
spec_bs = spec[spec['bindingsite_dist'] == 0]
spec_bs = spec_bs[active_mask.reindex(spec_bs.index, fill_value=False)] #remove variants which may be misfolded (below activity cutoff)
spec_nbs = spec[spec['bindingsite_dist'] != 0]
spec_nbs = spec_nbs[active_mask.reindex(spec_nbs.index, fill_value=False)] #remove variants which may be misfolded (below activity cutoff)

#Specificity-driving mutations within the binding site - stack
fig = pl.stack_plot(spec_bs, 6, 4, text_x = -0.21)
plt.show()
fig.savefig(f'{base_dir}/results/fig2/spec_bs_stack.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#Specificity-driving mutations outside the binding site - stack
fig = pl.stack_plot(spec_nbs, 6, 4, text_x = -0.21)
plt.show()
fig.savefig(f'{base_dir}/results/fig2/spec_nbs_stack.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 2C: BLOSUM90 SCORE DISTANCE DISTRIBUTIONS """
data['blosum90'] = data.index.to_series().apply(sp.physicochemical_distance)

fig, ax = plt.subplots(figsize = (7, 2))
sns.kdeplot(data[data['cluster'] == 7]['blosum90'], 
            fill = True, bw_adjust = 1.5, 
            label = 'Univ. disabling', alpha = .6,
            color = '#4F8DA7', edgecolor = 'black', linewidth = .75)
sns.kdeplot(data[data['cluster'] == 26]['blosum90'], 
            fill = True, bw_adjust = 1.5, 
            label = 'Univ. permitted', alpha = .4,
            color = '#dadcdf', edgecolor = 'black', linewidth = .75)
sns.kdeplot(data['blosum90'], 
            fill = False, bw_adjust = 1.8, 
            label = 'All data', color = 'gray', 
            linestyle = '--', linewidth = 1)

ks, p = stats.ks_2samp(data[data['cluster'] == 7]['blosum90'].dropna(), 
                       data[data['cluster'] == 26]['blosum90'].dropna())
bbox_props = dict(boxstyle="round,pad=0.5", edgecolor="black", linewidth = .75, facecolor=(1,1,1,.8))
ax.text(.31, .7, f'KS statistic = {ks:.2f}\np = {p:.2e}', ha = 'right', transform=plt.gca().transAxes, bbox = bbox_props)
ax.set_xlabel('BLOSUM90 mutation distance')

plt.legend()
plt.show()
fig.savefig(f'{base_dir}/results/fig2/blosum90_dist.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 2F: SPECIFICITY IMPORTANCE VS. DISTANCE FROM BINDING SITE """
bin_edges = np.concatenate(([int(data['bindingsite_dist'].max() + 1)], np.arange(20, -1, -1)))

#Mechanism-based bins ("Shell 1, Shell 2, Shell 3, beyond")
colors = ['#4f8da7', '#4f8da7', '#4f8da7', '#4f8da7', '#4f8da7', '#4f8da7', '#4f8da7', '#4f8da7', '#4f8da7', '#86B26E', '#86B26E', '#86B26E', '#86B26E', '#86B26E', '#86B26E', '#86B26E', '#e1b23a', '#e1b23a', '#e1b23a', '#e1b23a', '#EB4343']
binned_scores = []

#Print pymol commands for coloring structure figure
for i in range(len(bin_edges)-1):
    bin_mask = (data['bindingsite_dist'] < bin_edges[i]) & (data['bindingsite_dist'] >= bin_edges[i + 1]) & active_mask #remove variants which may be misfolded (below activity cutoff)
    bin_function_scores = np.array(data['SI'])[bin_mask]
    binned_scores.append(bin_function_scores[~np.isnan(bin_function_scores)])
    print(f'color 0x{colors[i][1:]}, resi {"+".join(data[bin_mask]["position"].unique().astype(str).tolist())}')

#Create boxplots
fig, ax = plt.subplots(figsize=(4, 7))
ticklabels = ['>20']+[f"{bin_edges[i+1]}-{bin_edges[i]}" for i in range(len(bin_edges) - 1)][1:-1] + ['0 ']
bp = ax.boxplot(binned_scores, tick_labels=ticklabels, showfliers=True, patch_artist=True,
                flierprops={'marker': 'o', 'markerfacecolor': '#4F8DA7', 'markersize': 5,
                            'alpha': 0.3}, vert=False)
ax.set_ylabel("Distance to nearest binding site residue (Å)", labelpad = 15)  
ax.text(1.9, 23.4, "Specificity importance", ha = 'center', fontsize = 11)
ax.text(1.9, 22.75, "(Deviations from median)", ha = 'center', fontsize = 9)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_yticklabels(ticklabels, rotation=0)
ax.set_xlim(-1.6, 5.5)

#Add shell labels
s1y = .975
s2y = .8575
s3y = .6
ax.annotate('Shell 1', xy=(1.025, s1y), xytext=(1.15, s1y), xycoords='axes fraction',
            fontsize=8, ha='center', va='center',
            arrowprops=dict(arrowstyle='-[, widthB=0.9, lengthB=.5', lw=.5, color='k'))
ax.annotate('Shell 2', xy=(1.025, s2y), xytext=(1.15, s2y), xycoords='axes fraction',
            fontsize=8, ha='center', va='center',
            arrowprops=dict(arrowstyle='-[, widthB=4.0, lengthB=.5', lw=.5, color='k'))
ax.annotate('Shell 3', xy=(1.025, s3y), xytext=(1.15, s3y), xycoords='axes fraction',
            fontsize=8, ha='center', va='center',
            arrowprops=dict(arrowstyle='-[, widthB=7.45, lengthB=.5', lw=.5, color='k'))

#Adjust style
for median in bp['medians']:
    median.set_color('black')
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
for marker, color in zip(bp['fliers'], colors):
    marker.set_markerfacecolor(color)
for spine in ['bottom', 'right']:
    ax.spines[spine].set_visible(False)

plt.show()
fig.savefig(f'{base_dir}/results/fig2/boxplot_dist.png',
            bbox_inches='tight', dpi=300, transparent = True)
