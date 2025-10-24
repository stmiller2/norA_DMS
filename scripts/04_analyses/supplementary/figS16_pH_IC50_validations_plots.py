#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 10:36:53 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import efficiency_utils as ef

#%%
""" FIG S16A: IC50 REGRESSIONS FOR ALL VARIANTS CLONALLY TESTED """

COLORS = {6: 'tab:green', 7: 'tab:purple'}

#Load data
data = pd.read_csv(f'{base_dir}/data/external/pH_ic50_data.csv')
strains = data['strain'].unique()
n_strains = len(strains)

#Logistic model to fit ic50 curves
def logistic(x, mn, mx, IC50, hill):
    return mn + (mx - mn) / (1 + (x / IC50) ** hill)

#Get curvefit data for strain/pH combo
def fit_curves(data):
    results = {}
    for pH in [6, 7]:
        pH_data = data[data['pH'] == pH]
        fits, ic50s = [], []

        for rep in pH_data['replicate'].unique():
            subset = pH_data[pH_data['replicate'] == rep]
            x, y = subset['norfloxacin_conc_ng/mL'], subset['OD600']

            guess = [min(y), max(y), np.median(x), 3.5]
            bounds = (0, [np.inf, np.inf, np.inf, 10])
            popt, _ = curve_fit(logistic, x, y, p0=guess, bounds=bounds, maxfev=1000)

            x_fit = np.logspace(0, np.log10(x.max()), 100)
            fits.append(logistic(x_fit, *popt))
            ic50s.append(popt[2])

        results[pH] = {
            "x_fit": x_fit,
            "fits": np.vstack(fits),
            "ic50s": np.array(ic50s),
            "scatter": pH_data
        }
    return results

#Plot
fig, axes = plt.subplots(3, 5, figsize=(15, 7), sharex=True, sharey=True)
axes = axes.flatten()
for ax, strain in zip(axes, strains):
    strain_data = data[data['strain'] == strain]
    results = fit_curves(strain_data)
    ph_ic50s = {}
    for pH, res in results.items():
        x_fit = res["x_fit"]
        mean, std = res["fits"].mean(axis=0), res["fits"].std(axis=0)

        #Scatter replicate points
        for rep in res["scatter"]['replicate'].unique():
            subset = res["scatter"][res["scatter"]['replicate'] == rep]
            ax.scatter(subset['norfloxacin_conc_ng/mL'], subset['OD600'],
                       color=COLORS[pH], s=15)

        #Curvefit mean with error
        ax.plot(x_fit, mean, color=COLORS[pH], label=f'pH {pH}')
        ax.fill_between(x_fit, mean-std, mean+std, color=COLORS[pH], alpha=0.3)
        ax.axvline(res["ic50s"].mean(), 0, mean.max()/2.3,
                   ls='--', color=COLORS[pH])
        ph_ic50s[pH] = res["ic50s"].mean()
    ic50_foldchange = ph_ic50s[6] / ph_ic50s[7]
    ax.text(800, .95, f'{ic50_foldchange:.1f}-fold\ndecrease', 
            ha = 'right', va = 'bottom', fontsize = 9,
            bbox=dict(boxstyle="round,pad=0.3",
                      edgecolor="black",
                      facecolor="white",
                      alpha=0.6))
    ax.set_xscale('log')
    ax.set_ylim(0, 1.21)
    ax.set_xlim(10, 1000)
    ax.set_title(strain, fontsize=14)
for i, ax in enumerate(axes):
    if i == 5:
        ax.set_ylabel('OD600 at 16 hours')
    if i == 12:
        ax.set_xlabel('[norfloxacin] (ng/mL)')
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(1.07, .95), loc='upper right')
fig.tight_layout()

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS16_pH_IC50_validations/all_IC50s.png', 
            dpi=300, bbox_inches='tight', transparent = True)

#%%
""" FIG S16: HIBIT MEMBRANE-LOCALIZED ABUNDANCE FOR OF SELECT EFFICIENCY MUTANTS """

#Load data
df = pd.read_csv(f'{base_dir}/data/external/hibit_data.csv')
df = df[['E222A', 'F188C', 'T336E', 'M150E', 'P121R', 'I76M', 'H354Q', 'A105E', 'F47D', 'S133E', 'F159K', 'H89L', 'N332Q', 'I298D', #Efficiency mutants
         'WT_batch1', 'WT_batch2', 'WT_batch3', 'WT_batch4', #WT
         'no_hibit', 'P144_G163del']] #Controls

#HiBiT experiments were performed in batches across five days. To normalize 
#day-to-day variability, wild type was included in every batch, and values are
#normalized to the wild type value recorded on that particular day. 
batches = {1: ["P144_G163del"],
           2: ["F47D", "I76M", "H89L", "A105E", "P121R", "S133E"],
           3: ["F159K", "F188C", "E222A", "I298D"],
           4: ["M150E", "N332Q", "T336E", "H354Q"]}

#Average technical replicates per biological replicate
bio_reps = ["rep1", "rep2", "rep3"]
bio_values = []
for i, rep in enumerate(bio_reps):
    tech_rows = df.iloc[i*3:(i+1)*3]
    bio_values.append(tech_rows.mean(axis=0))
bio_df = pd.DataFrame(bio_values, index=bio_reps)

#Subtract background (no hibit) per replicate
bio_df = bio_df.subtract(bio_df["no_hibit"], axis=0)
bio_df = bio_df.drop(columns="no_hibit")

#Normalize each replicate by its own batch's WT
for batch_num, variant_cols in batches.items():
    variant_cols = [col for col in variant_cols if col in bio_df.columns]
    wt_col = f"WT_batch{batch_num}"
    bio_df[[wt_col] + variant_cols] = (
        bio_df[[wt_col] + variant_cols].div(bio_df[wt_col], axis=0)
    )

#Keep only WT_batch1 as "WT" for plotting purposes
bio_df = bio_df.drop(columns=[c for c in bio_df.columns if c.startswith("WT_batch") and c != "WT_batch1"])
bio_df = bio_df.rename(columns={"WT_batch1": "WT"})

#Compute means and errors across normalized replicates
means = bio_df.mean(axis=0)
errs = bio_df.std(axis=0, ddof=1)

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
fig =  plt.figure(figsize=(6, 4))
plt.bar(labels, means, yerr=errs, capsize=3, color=colors,
        error_kw=dict(elinewidth=1))
x_positions = np.arange(len(labels))
plt.axhline(1, linestyle='--', color='gray', linewidth=1)
plt.fill_between([-1,len(means)], 1/3, 3, color='gray', alpha=0.2, zorder=0)
plt.text(-.5, 2.25, '3-fold WT expression', fontsize = 9)
plt.fill_between([-1,len(means)], 1/2, 2, color='gray', alpha=0.2, zorder=0)
plt.text(-.5, 1.5, '2-fold WT expression', fontsize = 9)

#Plot replicate dots
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
plt.xticks(rotation=45, ha='right', rotation_mode = 'anchor', fontsize = 9)
plt.yticks(fontsize = 10)
plt.ylabel("Relative protein abundance (log)")
plt.tight_layout()

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS16_pH_IC50_validations/hibit_efficiency.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S16C: SCATTER HIBIT AND EFFICIENCY DATA """

#Load data
hibit_data = pd.read_csv(f'{base_dir}/results/general/hibit_sum.csv', index_col = 0)
hibit_data = hibit_data.loc[hibit_data.index.drop(['P144_G163del', 'WT'])] #Remove non-library variants
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

#Set up minimum activity mask
active_mask = (spec[ligcols] > -1).sum(axis = 1)>0 #Only variants with F>-1 

#Plot
fig, ax = plt.subplots(figsize = (5,5), dpi = 300)
x = hibit_data.loc[active_mask, 'mean']
y = effc.loc[x.index, 'nor_dFpH']
ax.scatter(x, y, alpha = 0.8, color = '#4F8DA7')
ax.set_xlabel('Relative protein abundance (log)', fontsize = 12)
ax.set_xscale('log')
ax.set_ylabel('Efficiency (F$_{pH7.0}$ − F$_{pH6.0}$)', fontsize = 12)

#Display Spearman correlation
rho, pval = stats.spearmanr(x, y)
text_str = f"Spearman R = {rho:.2f}\np = {pval:.2f}"
bbox_props = dict(boxstyle="round,pad=0.5", edgecolor="black", linewidth = .75, facecolor=(1,1,1,.2))
ax.text(.92, 0.92, text_str, transform=ax.transAxes, fontsize=10, verticalalignment='top', ha = 'right', bbox=bbox_props)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS16_pH_IC50_validations/hibit_efficiency_scatter.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)
