#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 28 16:33:03 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import adjusted_rand_score
from scipy.cluster.hierarchy import linkage, fcluster
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import efficiency_utils as ef
import scoring_utils as sc #for correlation_scatter function

#Load data
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='exclude')
ligcols = data.columns
data['max'] = data[ligcols].max(axis = 1)

#Get cluster assignments
cluster_assignments = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv', index_col = 'mut')['cluster']
data['cluster'] = cluster_assignments

#%%
""" FIG S18A: HISTOGRAM BY CLUSTER """

#Plot
fig = plt.figure(figsize = (7,3))
bins = np.linspace(data['max'].min(), data['max'].max(), 100)
plt.hist(data['max'], bins=bins, alpha=1, histtype = 'stepfilled',
    edgecolor = 'black', label='All variants', color='white')
plt.hist(data.loc[data['cluster'] == 26, 'max'], bins=bins, alpha=1, histtype = 'stepfilled',
    edgecolor = 'black', linewidth = .5, label='Univ. permitted', color='#dadcdf')
plt.hist(data.loc[data['cluster'] == 7, 'max'], bins=bins, alpha=1, histtype = 'stepfilled',
    edgecolor = 'black', linewidth = .5, label='Univ. disabling', color='#4F8DA7')

#Annotate chosen minimal activity cutoff (-1)
plt.axvline(-1, color = 'red', linestyle = '--')
plt.xlabel('Standardized functional score (best-performing substrate)')
plt.ylabel('# of variants')
plt.legend()

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS18_sensitivity_analysis/histogram.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S18B: HIBIT ABUNDANCE VS MAXIMUM ACTIVITY """

#Load data
hibit = pd.read_csv(f'{base_dir}/results/general/hibit_sum.csv', index_col = 0)

#Plot
fig = sc.correlation_scatter(data['max'], hibit['mean'], 
                             xlab = 'Max F score', ylab = 'Protein abundance\n(HiBiT, WT-normalized)', 
                             yscale = 'log', corr = None, show = False, size = (4.5,3))
plt.axvline(-1, color = 'gray', linestyle = '--')
plt.xlim(-1.1)
plt.text(data.loc['E356P', 'max']*.9, hibit.loc['E356P', 'mean'], 'E356P', fontsize = 8)
plt.text(data.loc['M150E', 'max']*.6, hibit.loc['M150E', 'mean'], 'M150E', fontsize = 8)
plt.text(data.loc['F159K', 'max']*.8, hibit.loc['F159K', 'mean'], 'F159K', fontsize = 8)
plt.fill_between([-1.1, 2.4], 0.5, 2, alpha = 0.2, zorder = 0, color = 'gray')
plt.text(.97, .79, 'Expression within 2-fold of WT', fontsize = 8, ha = 'right', transform = plt.gca().transAxes)
plt.fill_between([-1.1, 2.4], 0.33, 3, alpha = 0.2, zorder = 0, color = 'gray')
plt.text(.97, .91, 'Expression within 3-fold of WT', fontsize = 8, ha = 'right', transform = plt.gca().transAxes)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS18_sensitivity_analysis/score_hibit_scatter.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S18C: OVERALL DATA RETENTION ACROSS DIFFERENT CUTOFFS """

#Compute data
cutoffs = np.arange(-2, 0.08, 0.08) 
frac_remaining = []
for cutoff in cutoffs:
    filtered = data[(data[ligcols].max(axis=1) > cutoff)].dropna()
    frac_remaining.append(len(filtered)/len(data.dropna()))

#Plot
fig, ax = plt.subplots(figsize=(7,3))
ax.plot(cutoffs, frac_remaining, linewidth = 2.5)
ax.set_ylabel('Proportion of data remaining')
ax.set_ylim(0,1)
ax.set_xlabel('Activity cutoff')
plt.title('Data remaining across activity cutoff stringencies')
plt.grid(True)
plt.show()

fig.savefig(f'{base_dir}/results/supplementary/figS18_sensitivity_analysis/data_remaining.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S18D: CLUSTER STABILITY ACROSS DIFFERENT CUTOFFS """

#Load data
cluster_assignments = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv', index_col = 'mut')
labels_orig = np.array([i[0] for i in cluster_assignments.dropna().reset_index(drop = True).values])
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=False, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='WT')
data = data.loc[cluster_assignments.dropna().index].reset_index(drop = True)

#Compute data
cutoffs = np.arange(-2, 0.08, 0.08) 
ari_scores = []
for cutoff in cutoffs:
    X_filtered = data[(data.max(axis=1) > cutoff)].dropna() 
    Z_filtered = linkage(X_filtered, method='average', metric='euclidean')
    n_clusters_filtered = min(round(np.sqrt(len(X_filtered)/2)), X_filtered.shape[0])
    labels_filtered = fcluster(Z_filtered, t=n_clusters_filtered, criterion='maxclust')
    overlap_idx = X_filtered.index
    ari = adjusted_rand_score(labels_orig[overlap_idx], labels_filtered)
    ari_scores.append(ari)

#Plot
fig, ax = plt.subplots(figsize=(7,3))
ax.plot(cutoffs, ari_scores, linewidth = 2.5)
ax.set_ylabel('Adjusted Rand Index\n(Similarity to unfiltered clustering)')
ax.set_ylim(0,1.02)  # ARI ranges 0-1
ax.set_xlabel('Activity cutoff')
plt.title('Cluster stability across activity cutoff stringencies')
plt.grid(True)
plt.show()

fig.savefig(f'{base_dir}/results/supplementary/figS18_sensitivity_analysis/cluster_ARI.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S18E: DISTAL SPECIFICITY-DRIVERS ACROSS DIFFERENT CUTOFFS """

# Load data
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='exclude')
dists = pd.read_csv(f'{base_dir}/results/general/residue_distances.csv', index_col = 'mut')
data['bindingsite_dist'] = dists[['9b3m_bindingsite', '7lo8_bindingsite']].min(axis = 1) #minimum distance to the binding site in either conformation
data.loc[data['bindingsite_dist'].isna(), 'bindingsite_dist'] = dists.loc[data['bindingsite_dist'].isna(), 'alphafold_bindingsite'] #use alphafold for site not resolved in cryoem structures

#Assign shells
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

#Compute data
cutoffs = np.arange(-2, 0.08, 0.08) 
spec_cutoff = data['SI'].quantile(.75) + (1.5*(data['SI'].quantile(.75) - data['SI'].quantile(.25))) #"Specificity-driving" means high outilers (above upper Tukey fence, ~2.77)
spec_mask = data['SI'] >= spec_cutoff

results = []
for cutoff in cutoffs:
    active_mask = (data[ligcols] > cutoff).sum(axis=1) > 0 #only variants with 1+ F>activity_cutoff
    subset = data[active_mask & spec_mask] #specificity-driving variants that meet the activity cutoff
    total_spec_variants = len(subset)
    # Proportion of specificity-driving variants within each shell
    for shell, group in subset.groupby('shell'):
        shell_variants = len(group)
        results.append((cutoff, shell, 100*shell_variants/total_spec_variants))
df_results = pd.DataFrame(results, columns=['cutoff', 'shell', 'prop'])

colors = ['#EB4343', '#e1b23a', '#86B26E']
labs = ['Shell 1 (0Å)', 'Shell 2 (1-5Å)', 'Shell 3 (5-12Å)']

#Plot
fig = plt.figure(figsize=(7,3))
for (_, grp), color, lab in zip(df_results.groupby('shell'), colors, labs):
    plt.plot(grp['cutoff'], grp['prop'], label=lab, color=color, linewidth = 2.5)
plt.ylim(0, 60)
plt.xlabel("Activity cutoff")
plt.ylabel("% of specificity-driving mutations")
plt.title("Frequency of distal specificity-driving mutations\nacross activity cutoff stringencies")
plt.legend(loc='upper left')
plt.grid(True)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS18_sensitivity_analysis/distal_shells.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" S18F; Correlation coefficient of specificity vs. efficiency """

#Load data
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'exclude')
spec = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='exclude')
spec['breadth'] = spec[ligcols].mean(axis = 1)

#Compute data
cutoffs = np.arange(-2, 0.08, 0.08) 
results_corr = []
for cutoff in cutoffs:
    active_mask = (spec[ligcols] > cutoff).sum(axis=1) > 0
    spec_subset = spec[active_mask]
    effc_subset = effc[active_mask]
    x = spec_subset['breadth']
    y = effc_subset['nor_dFpH']
    valid_mask = ~x.isna() & ~y.isna()
    rho, pval = stats.spearmanr(x[valid_mask], y[valid_mask])
    results_corr.append((cutoff, rho, pval))
df_corr = pd.DataFrame(results_corr, columns=['cutoff', 'rho', 'pval'])

#Plot
fig = plt.figure(figsize=(7,3))
plt.plot(df_corr['cutoff'], df_corr['rho'], linewidth=2.5)
plt.ylim(0,1)
plt.xlabel("Activity cutoff")
plt.ylabel("Spearman R")
plt.title("Efficiency (ΔF$_{pH}$) correlation with promiscuity (F$_{avg}$)\nacross activity cutoff stringencies")
plt.grid(True)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS18_sensitivity_analysis/promiscuity_efficiency_corr.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

