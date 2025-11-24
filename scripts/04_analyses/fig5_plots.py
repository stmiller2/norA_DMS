#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:03:18 2024

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import matplotlib as mpl
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import efficiency_utils as ef

""" Load data """

spec = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=False, 
                                return_format='dictionary', 
                                wt_policy='WT')

#get custom colormap
from plot_utils import spec_cmap, spec_vmin, spec_norm

#%%
""" FIG 5A: ETHIDIUM EFFLUX TIMECOURSE PLOTS """

#Get ethidium efflux platereader data
sheets = ['rep1', 'rep2', 'rep3'] #one replicate on each sheet
reps = []
for sh in sheets:
    df = pd.read_excel(f'{base_dir}/data/external/ethidium_efflux_data.xlsx', 
                       sheet_name=sh)
    for col in df.columns[1:]:
        df[col] = 100 * df[col] / df[col].iloc[:4].mean() #normalize to starting fluorescence
    df['Time'] = pd.to_timedelta(df['Time'].astype(str)).dt.total_seconds()
    reps.append(df.set_index('Time'))

#Combine reps & get error
rep_stack = np.stack([r.values for r in reps], axis=-1)
mean_vals = rep_stack.mean(axis=-1)
err_vals  = rep_stack.std(axis=-1, ddof=1) / np.sqrt(3) #SEM
mean_df = pd.DataFrame(mean_vals, index=reps[0].index, columns=reps[0].columns)
err_df  = pd.DataFrame(err_vals,  index=reps[0].index, columns=reps[0].columns)

#Get growth-based scores
variant_scores = spec['eth'].loc[mean_df.columns, ['f_hiq_norm', 'SE_norm']]

#Plot
fig, ax = plt.subplots(figsize=(5, 4))
for v in mean_df.columns:
    color = spec_cmap(spec_norm(variant_scores.loc[v, 'f_hiq_norm']))
    ax.plot(mean_df.index, mean_df[v], color=color, linewidth=2.5)
    ax.fill_between(mean_df.index,
                    mean_df[v] - err_df[v],
                    mean_df[v] + err_df[v],
                    color=color, alpha=0.2)
ax.annotate("0.4% glucose\nadded",
            xy=(136, mean_df.loc[136].mean()*.98),   # point on the curves
            xytext=(175, ax.get_ylim()[1]*0.75), # position of text
            arrowprops=dict(facecolor='black', arrowstyle="->"),
            fontsize=10, ha='center')
ax.set_xlabel("Time (s)", fontsize = 12)
ax.set_ylabel("Intracellular ethidium\n(% starting fluorescence)", fontsize = 12)
ax.set_title("Clonal ethidium efflux")
#Colorbar
sm = mpl.cm.ScalarMappable(cmap=spec_cmap, norm=spec_norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.ax.set_ylim(spec_vmin, 1)
cbar.set_label("Growth-based score ($F_{ethidium}$)", fontsize = 11)

plt.show()
fig.savefig(f'{base_dir}/results/fig5/ethidium_timecourse.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 5B: CORRELATION OF AREA UNDER EFFLUX CURVE WITH GROWTH BASED SCORES """

#Compute AUC
results = []
for i, v in enumerate(mean_df.columns):
    vals = rep_stack[:, i, :]
    aucs = []
    for r in range(vals.shape[1]):
        y = vals[:, r]; t = mean_df.index.values
        aucs.append(np.trapz(y, t))
    results.append({
        'variant'    : v,
        'score'      : variant_scores.loc[v, 'f_hiq_norm'],
        'score_err'  : variant_scores.loc[v, 'SE_norm'],
        'auc_mean'   : np.nanmean(aucs),
        'auc_err'    : np.nanstd(aucs, ddof=1)/np.sqrt(len(aucs)), #SEM
    })
res_df = pd.DataFrame(results)

#Plot correlation
fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.errorbar(res_df['score'], res_df['auc_mean'], 
            xerr=res_df['score_err'], yerr=res_df['auc_err'], 
            fmt='o', color='black', alpha = 0.8, ecolor='gray', capsize=3)

#Label variant names
labs = []
for xi, yi, name in zip(res_df['score'], res_df['auc_mean'], res_df['variant']):
    labs.append(ax.text(xi, yi, name, fontsize=8, ha = 'right'))
adjust_text(labs, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

#Display R value
r, pval = stats.spearmanr(res_df['score'], res_df['auc_mean'], nan_policy='omit')
props = dict(boxstyle="round,pad=0.5", edgecolor = 'black', linewidth = .75, facecolor=(1,1,1,.2))
ax.text(0.93, 0.93, f"Spearman R = {r:.2f}\np = {pval:.2e}",
        transform = ax.transAxes, fontsize = 8.5,
        verticalalignment = 'top', ha = 'right',
        bbox = props)

#Adjust style
ax.set_xlabel('Growth-based score ($F_{ethidium}$)', fontsize = 12)
ax.set_ylabel("Area under efflux curve", labelpad = 13, fontsize = 12)
ax.text(-0.125, 0.05, 'More\nefflux', ha = 'right', va = 'center',
        transform = ax.transAxes, fontsize = 8)
ax.text(-0.125, 0.95, 'Less\nefflux', ha = 'right', va = 'center',
        transform = ax.transAxes, fontsize = 8)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

plt.show()
fig.savefig(f'{base_dir}/results/fig5/ethidium_scatter.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG 5C: CORRELATION OF HIGH-THROUGHPUT AND CLONAL PH SENSITIVITIES """

clones = ['WT', 'T336E', 'H89L', 'F47D', 'I298D', 'N332Q', 'S133E', 'F159K', 'H354Q', 'I76M', 'P121R', 'F188C', 'A105E', 'M150E']

#Get data
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'WT')

clnl = pd.read_csv(f'{base_dir}/results/general/clonal_dFpH.csv',
                   index_col = 'strain')
clnl = clnl.loc[clones]
effc = effc.loc[clones]

r, p = stats.spearmanr(effc['nor_dFpH'], clnl['dFpH'])

#Plot
fig, ax = plt.subplots(figsize = (3.5,3.5))
ax.errorbar(effc['nor_dFpH'], clnl['dFpH'], 
            xerr=effc['nor_SE'], yerr=clnl['err_dFpH'], 
            fmt='o', color='black', alpha=0.8, ecolor='gray', capsize=3)
ax.set_xlabel('High-throughput ΔF$_{pH}$', fontsize = 12)
ax.set_ylabel('Clonal ΔF$_{pH}$', fontsize = 12)
ax.set_ylim(-.7, .35)
labs = []
for clone_name in clones:
    labs.append(ax.text(effc['nor_dFpH'].loc[clone_name], clnl['dFpH'].loc[clone_name], clone_name, fontsize=8, ha='right'))

props = dict(boxstyle="round,pad=0.5", edgecolor = 'black', linewidth = .75, facecolor=(1,1,1,.2))
ax.text(0.055, 0.95, f"Spearman R = {r:.2f}\np = {p:.3f}",
        transform = ax.transAxes, fontsize = 8.5,
        verticalalignment = 'top', ha = 'left',
        bbox = props)
adjust_text(labs, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))  
plt.show()

fig.savefig(f'{base_dir}/results/fig5/validation_scatter.png', 
            dpi = 300, bbox_inches = 'tight', transparent = True)

#%%
""" FIG 5D: HIBIT MEMBRANE-LOCALIZED ABUNDANCE BAR PLOT """

#Load data
df = pd.read_csv(f'{base_dir}/data/external/hibit_data.csv')

#HiBiT experiments were performed in batches across five days. To normalize 
#day-to-day variability, wild type was included in every batch, and values are
#normalized to the wild type value recorded on that particular day. 
batches = {1: ["P144_G163del", "I15N", "L30N", "P110S", "T211W", "Y278F"],
           2: ["F47D", "I76M", "H89L", "A105E", "G111V", "P121R", "S133E", "G139R"],
           3: ["E82Q", "F159K", "I177W", "F188C", "E222A", "I244F", "Q255W", "Y278I", "I298D", "F303Y", "F306G"],
           4: ["M150E", "R310N", "P311M", "N332Q", "T336E", "H354Q", "E356P"],
           5: ["P144A", "G147W", "G251D"]}

#Average technical replicates per biological replicate
bio_reps = ["rep1", "rep2", "rep3"]
bio_values = []
for i, rep in enumerate(bio_reps):
    tech_rows = df.iloc[i*3:(i+1)*3][df.columns[1:]] #first column titles replicates, can remove
    bio_values.append(tech_rows.mean(axis=0))
bio_df = pd.DataFrame(bio_values, index=bio_reps)

#Subtract background (no hibit) per replicate
bio_df = bio_df.subtract(bio_df["no_hibit"], axis=0)
bio_df = bio_df.drop(columns="no_hibit")

#Normalize each replicate by its own batch's WT
for batch_num, variant_cols in batches.items():
    variant_cols = [col for col in variant_cols if col in bio_df.columns]
    wt_col = f"WT_batch{batch_num}"
    bio_df[[wt_col] + variant_cols] = (bio_df[[wt_col] + variant_cols].div(bio_df[wt_col], axis=0))

#Keep only WT_batch1 as "WT" for plotting purposes
bio_df = bio_df.drop(columns=[c for c in bio_df.columns if c.startswith("WT_batch") and c != "WT_batch1"])
bio_df = bio_df.rename(columns={"WT_batch1": "WT"})

#Compute means and errors across normalized replicates
means = bio_df.mean(axis=0)
errs = bio_df.std(axis=0, ddof=1) / np.sqrt(3) #SEM

#Sort for plotting
means = means.sort_values()
errs = errs.loc[means.index]

labels = [r"$\bf{WT}$" if l == "WT" else
          r"$\bf{P144\_G163del}$" if l == "P144_G163del" 
          else l for l in 
          [label.replace("WT_batch1", "WT") for label in means.index]]

colors = ['#69B572' if l == r"$\bf{WT}$" else
          '#E1B23A' if l == r"$\bf{P144\_G163del}$" else 
          '#4F8DA7' for l in labels]

#Plot
fig =  plt.figure(figsize=(6, 3.5))
plt.bar(labels, means, yerr=errs, capsize=3, color=colors,
        error_kw=dict(elinewidth=1))
x_positions = np.arange(len(labels))
plt.axhline(1, linestyle='--', color='gray', linewidth=1)
plt.fill_between([-1,len(means)], 1/3, 3, color='gray', alpha=0.2, zorder=0)
plt.text(-.5, 2.25, '3-fold WT expression', fontsize = 9)
plt.fill_between([-1,len(means)], 1/2, 2, color='gray', alpha=0.2, zorder=0)
plt.text(-.5, 1.5, '2-fold WT expression', fontsize = 9)

#Plot replicate scatter points
for i, label in enumerate(means.index):
    if label not in bio_df.columns:
        continue
    for rep in bio_df.index:
        y = bio_df.loc[rep, label]
        plt.scatter(i + np.random.uniform(-0.15, 0.15), y,
                    color='black', alpha=.3, s=15)

#Adjust style
plt.xlim(-1, len(means))
plt.yscale('log')
plt.ylim(np.log(1.028), np.log(1000))
plt.xticks(rotation=45, ha='right', rotation_mode = 'anchor', fontsize = 8)
plt.yticks(fontsize = 10)
plt.ylabel("Normalized luminescence\n(log)", fontsize = 10, labelpad = -10)
plt.title("Relative protein abundance", fontsize = 11)
plt.tight_layout()

plt.show()
fig.savefig(f'{base_dir}/results/fig5/hibit_bar.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#Save summarized HiBiT data
hibit_sum = pd.DataFrame()
hibit_sum['mean'] = means
hibit_sum['err'] = errs
hibit_sum.to_csv(f'{base_dir}/results/general/hibit_sum.csv')

#%%
""" FIG 5E: CORRELATION HEATMAP OF SPECIFICITY, EFFICIENCY, ABUNDANCE """

# Get specificity data
spec = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')
spec.index = spec.index.str[1:]

#set up minimum activity mask
active_mask = (spec > -1).sum(axis = 1)>0 #Only variants with 1+ F>-1

spec['breadth'] = spec.mean(axis = 1) #get Favg

#Get efficiency data
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'exclude')
effc.index = effc.index.str[1:]

#Get Rosetta & ThermoMPNN predicted ddGs
computational_data = pd.read_csv(f'{base_dir}/data/external/computational_stability.csv',
                                 index_col = 0)

#HiBiT expression data
hibit_data = pd.read_csv(f'{base_dir}/results/general/hibit_sum.csv', index_col = 0)
hibit_data = hibit_data[~hibit_data.index.isin(['144_G163del', 'WT'])] #Remove non-library variants
hibit_data.index = hibit_data.index.str[1:]

#Assemble master dataframe
metrics = pd.DataFrame()
metrics['promiscuity'] = spec['breadth']
metrics['efficiency'] = effc['nor_dFpH']
metrics['rosetta'] = computational_data['rosetta_min']
metrics['thermo'] = computational_data['thermo_min']
metrics['hibit'] = hibit_data['mean']
corr_matrix = metrics.loc[active_mask].corr(method='spearman')


labs = ['Promiscuity\n($F_{avg}$)',
        'Efficiency\n($F_{pH7.0}$ - $F_{pH6.0}$)',
        'Rosetta\n(pred. ΔΔG)',
        'ThermoMPNN\n(pred. ΔΔG)',
        'Protein\nabundance\n(HiBiT lum.)']
labs_short = ['Promiscuity',
              'Efficiency',
              'Rosetta',
              'ThermoMPNN',
              'Protein\nabundance']

#Plot
fig, ax = plt.subplots(figsize = (4, 4))
heatmap = sns.heatmap(corr_matrix, annot=True, fmt='.2f', linewidths=0.25,
                      cmap='coolwarm', vmin=-1, vmax=1, xticklabels=labs_short, yticklabels=labs,
                      cbar = None, annot_kws={"size":12})
plt.title('Spearman correlation', fontsize = 14, y = 1.02)
plt.xticks(rotation=25, rotation_mode='anchor', ha = 'right', fontsize = 12)
plt.yticks(fontsize = 12)
ax.axhline(2, color='white', linewidth=3)
ax.axvline(2, color='white', linewidth=3)

plt.show()
fig.savefig(f'{base_dir}/results/fig5/abundance_correlations.png', 
            dpi = 300, bbox_inches = 'tight', transparent = True)
