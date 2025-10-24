#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 08:29:09 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
os.chdir(f'{base_dir}/scripts/00_utils')
from plot_utils import sub_colors

#%%
""" FIG S1A: IC50 CURVES FOR ALL SUBSTRATES """

#load data
data = pd.read_csv(f'{base_dir}/data/external/substrates_ic50_data.csv')

#logistic model to fit ic50 curves
def logistic(x, mn, mx, IC50, hill):
    return mn + (mx - mn) / (1 + (x / IC50) ** hill)

#fit logistic curves per strain for each drug
def fit_curves(data):
    results = {}
    for strain in data['strain'].unique():
        strain_data = data[data['strain'] == strain]
        fits, ic50s = [], []

        for rep in strain_data['rep'].unique():
            subset = strain_data[strain_data['rep'] == rep]
            x, y = subset['conc'], subset['OD600_blank_subtracted']

            guess = [min(y), max(y), np.median(x), 3.5]
            bounds = (0, [np.inf, np.inf, np.inf, 10])
            popt, _ = curve_fit(logistic, x, y, p0=guess, bounds=bounds, maxfev=2000)

            x_fit = np.logspace(np.log10(subset.loc[subset['conc'] > 0, 'conc'].min()), np.log10(x.max()), 100)
            fits.append(logistic(x_fit, *popt))
            ic50s.append(popt[2])

        if fits:
            results[strain] = {
                "x_fit": x_fit,
                "fits": np.vstack(fits),
                "ic50s": np.array(ic50s),
                "scatter": strain_data
            }
    return results

#Plot
COLORS = {'WT': '#4F8DA7', 'E222A': '#E44B3A'}
drugs = data['drug'].unique()

fig, axes = plt.subplots(2, 4, figsize=(12, 6.5), sharex=False, sharey=True)
axes = axes.flatten()
foldchanges = {}
for ax, drug in zip(axes, drugs):
    drug_data = data[data['drug'] == drug]
    results = fit_curves(drug_data)
    ic50s = {}
    for strain, res in results.items():
        x_fit = res["x_fit"]
        mean, sem = res["fits"].mean(axis=0), res["fits"].std(axis=0)/np.sqrt(3) #SEM

        #Scatter replicate points
        for rep in res["scatter"]['rep'].unique():
            subset = res["scatter"][res["scatter"]['rep'] == rep]
            ax.scatter(subset['conc'], subset['OD600_blank_subtracted'], color=COLORS[strain], s=15)

        #Mean fit + error
        ax.plot(x_fit, mean, color=COLORS[strain], label=strain)
        ax.fill_between(x_fit, mean - sem, mean + sem, color=COLORS[strain], alpha=0.3)
        ax.axvline(res["ic50s"].mean(), 0, mean.max()/2.3, ls='--', color=COLORS[strain])
        ic50s[strain] = res["ic50s"].mean()

    #Record fold-changes for bar plot
    fold = ic50s['WT'] / ic50s['E222A']
    foldchanges[drug] = fold

    #Adjust style
    ax.set_xscale('log')
    ax.set_ylim(-0.02, 1.2)
    ax.set_xlim(min(drug_data['conc'][drug_data['conc']>0]), max(drug_data['conc']))
    title = drug if drug == 'TPP' else drug.capitalize()
    ax.set_title(title, fontsize=12)

#Shared labels
fig.text(0.5, 0.04, '[drug] (µg/mL)', ha='center')
fig.text(0.04, 0.5, 'OD600 at 16 hours', va='center', rotation='vertical')

#Legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(1.01, 0.95), loc='upper left')
fig.tight_layout(rect=[0.05, 0.05, 1, 1])

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS03_substrate_ic50s/curves.png',
            dpi = 300, bbox_inches = 'tight', transparent = True)

#%%
""" FIG S1B: QUANTIFIED FOLD-CHANGE IN IC50S FOR EACH SUBSTRATE """

#Plot
fig, ax = plt.subplots(figsize=(3, 5))
ax.bar(foldchanges.keys(), foldchanges.values(), color=sub_colors, edgecolor='black')
ax.set_ylabel('$IC_{50}$ fold change (WT / E222A)')
ax.set_title('Drug resistance conferred by NorA')
ax.set_ylim(0, max(foldchanges.values()) * 1.2)
plt.xticks(rotation=45, ha='right', rotation_mode = 'anchor')

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS03_substrate_ic50s/bar.png',
            dpi = 300, bbox_inches = 'tight', transparent = True)


