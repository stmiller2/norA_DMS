#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 09:44:20 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
from sklearn.neighbors import NearestCentroid
import seaborn as sns
from matplotlib import pyplot as plt

def cophenetic_coefficient(data, distance_metric, linkage_method):
    pairwise_distances = pdist(data, metric=distance_metric)
    dendrogram = hierarchy.linkage(data, 
                                   metric=distance_metric, 
                                   method=linkage_method)
    cophenetic_distances = hierarchy.cophenet(dendrogram)
    correlation = np.corrcoef(pairwise_distances, cophenetic_distances)[0, 1]
    return correlation

def gap_statistic(X, clusterer, k_values, n_refs=20, random_state=None):
    rng = np.random.default_rng(random_state)
    shape = X.shape
    tops = X.max(axis=0)
    bottoms = X.min(axis=0)
    dists = np.zeros((len(k_values), n_refs + 1))

    for i, k in enumerate(k_values):
        if k%10 == 0:
            print(f'K value {k} complete')
        #Fit clustering on actual data
        _, labels = clusterer(X, k)
        Wk = _within_cluster_dispersion(X, labels)
        dists[i, 0] = np.log(Wk)

        #Fit clustering on reference datasets
        for j in range(n_refs):
            X_ref = rng.uniform(bottoms, tops, size=shape)
            _, labels_ref = clusterer(X_ref, k)
            Wk_ref = _within_cluster_dispersion(X_ref, labels_ref)
            dists[i, j + 1] = np.log(Wk_ref)

    #Compute gap statistic
    gaps = dists[:, 1:].mean(axis=1) - dists[:, 0]
    sk = dists[:, 1:].std(axis=1) * np.sqrt(1 + 1 / n_refs)
    return gaps, sk

from functools import lru_cache

@lru_cache(maxsize=None)
def precompute_linkage(X_tuple, metric, method):
    X = np.array(X_tuple)
    return hierarchy.linkage(X, metric=metric, method=method, optimal_ordering=False)

def agglom_prelinked(X, k, linkage_matrix):
    y_predict = hierarchy.fcluster(linkage_matrix, t=k, criterion='maxclust') - 1
    clf = NearestCentroid()
    clf.fit(X, y_predict)
    return clf.centroids_, y_predict

def _within_cluster_dispersion(X, labels):
    #vectorized version
    unique, inverse = np.unique(labels, return_inverse=True)
    counts = np.bincount(inverse)
    sums = np.zeros((len(unique), X.shape[1]))
    np.add.at(sums, inverse, X)
    centroids = sums / counts[:, None]
    diffs = X - centroids[inverse]
    return np.sum(np.linalg.norm(diffs, axis=1))

def split_snsclustermap(data, 
                     datacols, 
                     metric, 
                     method, 
                     nclusters = None,
                     cmap = 'coolwarm', 
                     center = 0,
                     figsize = (5,15),
                     cbar_label = 'Functional score',
                     tree_linewidth = 1.5,
                     split = False,
                     split_linewidth = 2,
                     tree_ratio = 0.2,
                     title = None,
                     vmin = None,
                     vmax = None,
                     col_cluster = True,
                     show_cbar = True):
    """
    Generate seaborn clustermap with colored clusters and (optional) lines 
    between clusters

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe containing the data.
    datacols : list
        List of columns in the dataframe that you want to plot.
    metric : string
        Distance metric for clustering.
    method : string
        Linkage method for clustering.
    nclusters : int, optional
        Number of clusters to split into, if you want to color clusters. The 
        default is None.
    cmap : string, optional
        Colormap for the heatmap. The default is 'coolwarm'.
    center : float, optional
        Where to center the heatmap color scale. The default is 0.
    figsize : tuple, optional
        (x, y) dimensions of the figure. The default is (5,15).
    cbar_label : string, optional
        Label for the colorbar. The default is 'Functional score'.
    tree_linewidth : float, optional
        Thickness of the dendrogram lines. The default is 1.5.
    split : bool, optional
        If true, lines will be drawn between clusters to further emphasize 
        their separation. Note that these lines are drawn on top of the heatmap
        and therefore will cover up a couple rows, depending on split_linewidth.
        The default is False.
    split_linewidth : float, optional
        Thickness of the cluster split lines. The default is 2.
    title : string, optional
        Title for the figure. The default is None.

    Returns
    -------
    fig : matplotlib.pyplot.figure
        Figure containing clustermap with your specifications.

    """
    if nclusters:
        dendrogram = hierarchy.linkage(data[datacols].dropna(), 
                                       metric = metric, 
                                       method = method, 
                                       optimal_ordering = False)
        data.loc[data[datacols].dropna().index, 'cluster'] = hierarchy.fcluster(dendrogram, 
                                                                                t = nclusters, 
                                                                                criterion = 'maxclust')
        cluster_sizes = [len(data[data['cluster'] == cluster]) for cluster in range(1,nclusters+1)]
    
    figx, figy = figsize
    fig = sns.clustermap(data[datacols], 
                         metric = metric,
                         method = method,
                         cmap = cmap,
                         center = center,
                         vmin = vmin, vmax = vmax,
                         col_cluster = col_cluster,
                         row_colors = plt.get_cmap('tab10')(data['cluster']%10) if nclusters else None,
                         figsize = figsize,
                         dendrogram_ratio = (tree_ratio, (tree_ratio / (figy / figx))),
                         cbar_pos = (0.85, 0.2, .05, .5) if show_cbar else None,
                         cbar_kws = {'label': cbar_label},
                         tree_kws = {'linewidth': tree_linewidth})
    hm = fig.ax_heatmap
    hm.set_ylabel(None)
    hm.set_yticks([])
    hm.set_yticklabels([])
    plt.setp(hm.xaxis.get_majorticklabels(), rotation = 30, ha = 'right')
    
    if title:
        plt.suptitle(title, y = 1.01)
    else:
        if nclusters:
            plt.suptitle(f'{metric} {method}, nclusters = {nclusters}', y = 1.01)
        else:
            plt.suptitle(f'{metric} {method}', y = 1.01)
    
    if split:
        for idx in [sum(cluster_sizes[:i]) for i in range(1, len(cluster_sizes)+1)]:
            hm.axhline(y=idx, color='white', linewidth=split_linewidth)
    
    return fig

def assign_clusters(dataframe, 
                    distance_metric, 
                    linkage_method, 
                    nclust, 
                    ligcols = ['acr', 'eth', 'nor', 'ofl', 'pen', 'ppa', 'pur', 'tpp']):
    #Compute dendrogram
    dendrogram = hierarchy.linkage(dataframe[ligcols].dropna(), metric = distance_metric, method = linkage_method, optimal_ordering = False)
    #Save clusters to dataframe
    dataframe.loc[dataframe[ligcols].dropna().index, 'cluster'] = hierarchy.fcluster(dendrogram, t = nclust, criterion = 'maxclust')
    return dataframe