#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:57:34 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.cluster import hierarchy
from itertools import product
import os

os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import clustering_utils as cl
import plot_utils as pl

### Load data
data = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = False,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')

data = data.dropna()
ligcols = ['acr', 'eth', 'nor', 'ofl', 'pen', 'ppa', 'pur', 'tpp']

distance_metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 
                    'cosine', 'dice', 'euclidean', 'mahalanobis', 'minkowski', 
                    'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 
                    'sokalsneath', 'sqeuclidean', 'yule'] 
linkage_methods = ['single', 'average', 'complete', 'median', 'ward', 'weighted', 
                   'centroid']
combinations = list(product(distance_metrics, linkage_methods))

### Check cophenetic correlations of metric/method combinations
ccs = pd.DataFrame()
ccs.index = [f'{distance_metric}_{linkage_method}' for distance_metric, linkage_method in combinations]
for distance_metric, linkage_method in combinations:
    try:
        cc = cl.cophenetic_coefficient(data[ligcols], distance_metric, linkage_method)
        ccs.loc[f'{distance_metric}_{linkage_method}', 'f_cophenet_coeff'] = cc
        print(f'{distance_metric} {linkage_method} complete')
    except:
        pass
        print(f'{distance_metric} {linkage_method} skipped')

fig, ax = plt.subplots()
ccs = ccs.sort_values('f_cophenet_coeff', ascending = False)
mask = ccs['f_cophenet_coeff'] > 0.7
ax.bar(ccs[mask].index, ccs.loc[mask, 'f_cophenet_coeff'])
ax.set_xticklabels(ccs[mask].index, rotation = 45, fontsize = 8, ha = 'right')
ax.set_ylim(0.65, 1)
ax.set_ylabel('Cophenetic correlation coefficient')
ax.set_xlabel('Clustering method')
ax.set_title('Cophenetic correlations > 0.70')
plt.show()
fig.savefig(f'{base_dir}/results/clustering/cophenetic_correlations.png',
            bbox_inches = 'tight', dpi = 300)

### Check gap statistic for various nclust and plot
metric = 'euclidean'
method = 'average'

X = data[ligcols].values
linkage_matrix = cl.precompute_linkage(tuple(map(tuple, X)), metric, method)

def agglom(X, k):
    return cl.agglom_prelinked(X, k, linkage_matrix)

k_values = np.arange(2, 101)
gaps, sk = cl.gap_statistic(data[ligcols].values, agglom, k_values, n_refs=20, random_state=42)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(k_values, gaps, linewidth=1)
ax.scatter(k_values, gaps, color='black', s=10)
ax.axvline(54, color = 'red', linestyle = '--') #annotate nclust chosen based on gap stat & domain knowledge
ax.text(55, 1.2, 'nclust = 54', fontsize = 10)
ax.set_ylabel('Gap statistic')
ax.set_xlabel('Number of clusters')
ax.set_xticks(np.arange(2, 101, 3))
ax.set_xticklabels(np.arange(2, 101, 3), rotation = 45, fontsize = 10)
ax.set_title(f'Gap statistic analysis\nClustering metric/method: {metric} {method}')
plt.show()
fig.savefig(f'{base_dir}/results/clustering/gap_statistic.png',
            bbox_inches = 'tight', dpi = 300)

### View clustering with chosen clustering parameters 
metric = 'euclidean'; method = 'average'; nclust = 54

fig = cl.split_snsclustermap(data = data, 
                             datacols = ligcols, 
                             metric = metric, 
                             method = method,
                             split = False,
                             tree_linewidth = 1,
                             figsize = (5, 10),
                             title = 'Clustering on Euclidean distance, average linkage')
plt.show()
fig.savefig(f'{base_dir}/results/clustering/clustermap_uncut.png',
            bbox_inches = 'tight', dpi = 300)

### Plot uncut dendrogram
dendrogram = hierarchy.linkage(data[ligcols], metric = metric, method = method, optimal_ordering = False)
fig, ax = plt.subplots(figsize=(7, 10))
hierarchy.dendrogram(dendrogram)
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_xlabel('Variants')
ax.set_ylabel('Distance')
plt.show()
fig.savefig(f'{base_dir}/results/clustering/dendrogram_uncut.png',
            bbox_inches = 'tight', dpi = 300)

### Plot cut dendrogram 
fig, ax = plt.subplots(figsize=(7, 10))
hierarchy.dendrogram(dendrogram, color_threshold=dendrogram[-nclust, 2], above_threshold_color='k')
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_xlabel('Variants')
ax.set_ylabel('Distance')
plt.axhline(y=dendrogram[-nclust, 2], color='r', linestyle='--', label=f'Cluster Threshold ({nclust} clusters)')
plt.show()
fig.savefig(f'{base_dir}/results/clustering/dendrogram_cut.png',
            bbox_inches = 'tight', dpi = 300)

### Save cluster assignments
metric = 'euclidean'; method = 'average'; nclust = 54
data = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = False,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'WT')

cluster_assignments = cl.assign_clusters(dataframe = data, 
                                         distance_metric = metric, 
                                         linkage_method = method, 
                                         nclust = nclust)['cluster']
cluster_assignments.to_csv(f'{base_dir}/results/clustering/cluster_assignments.csv')

### View all clusters with at least 3 members
min_cluster_size = 3
mask = data['cluster'].value_counts() >= min_cluster_size
clusters = data['cluster'].value_counts()[mask].index
data_by_cluster = {cluster: data[data['cluster'] == cluster] for cluster in clusters}
keys = list(data_by_cluster.keys())
keys.sort()
data_by_cluster = {i: data_by_cluster[i] for i in keys} #sort by cluster id number

#plot
cmap, vmin, vmax = pl.custom_colormap(data[ligcols])
fig, axes = plt.subplots(nrows=len(clusters), figsize = (4, 12))
cbar_ax = fig.add_axes([0.9, 0.05, .05, .2])
for idx, (cluster, data_cluster) in enumerate(data_by_cluster.items()):
    sns.heatmap(data_cluster[ligcols], cmap=cmap, center = 0, 
                ax=axes[idx], vmin = vmin, vmax = vmax,
                cbar = (idx == 0), cbar_ax = None if idx else cbar_ax,
                cbar_kws = None if idx else {'label': 'Functional score'})    
    axes[idx].set_yticks([])
    axes[idx].set_ylabel(f'c{cluster} (N = {len(data_cluster)})', 
                         rotation = 0, 
                         fontsize = 12, 
                         ha = 'right')
    axes[idx].yaxis.set_label_coords(-0.01,0.05)
    if idx == 0:
        axes[idx].set_title(f'{metric.capitalize()} {method}, {nclust} clusters\nMinimum cluster size: {min_cluster_size}')
    
    if idx == len(data_by_cluster.items())-1:
        axes[idx].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
        axes[idx].set_xticklabels(ligcols, rotation = 30, ha = 'right',
                                  rotation_mode = 'anchor', fontsize = 12)
    else:
        axes[idx].set_xticks([])

plt.subplots_adjust(top=0.93,
                    bottom=0.05,
                    left=0.08,
                    right=0.84,
                    hspace=0.2,
                    wspace=0.2)
plt.show()
fig.savefig(f'{base_dir}/results/clustering/all_clusters.png',
            bbox_inches = 'tight', dpi = 300)

