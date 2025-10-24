#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 15:10:42 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import ceil
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import barcode_mapping_utils as bcm

#%%
""" REPEAT BARCODING FOR REPRESENTATIVE SUBPOOL (2) """

# Load raw barcode mapping sequencing data (using tile 2 as representative)
seq_df = pd.read_csv(f'{base_dir}/data/processed/barcoding/good_reads/norA_T2_good_reads.csv', names = ['read'])
ref_lib = pd.read_csv(f'{base_dir}/data/processed/barcoding/ref_libs/norA_T2_ref_lib.csv', header = None)

# Define handles for grabbing barcode and ORF
prebc = 'GTAACCGCCACG'
postbc = 'TGAGATCCGGCT'
prerf = 'TTAATTATATGT'
postrf = 'CATCGTATGCCA'

### Extract barcodes
seq_df = bcm.extract_barcodes(seq_df, prebc, postbc)

### Translate variants
seq_df = bcm.translate_variants(seq_df, prerf, postrf)

### Group barcodes by levenshtein distance = 1
bc_groups = bcm.group_barcodes(seq_df)

### Collapse barcode groups
seq_df = bcm.collapse_groups(seq_df, bc_groups)

### Sample read cutoffs and filter
cutoff_sample = 25
pair_cutoff = 6

pass_pc = seq_df.copy()
pass_pc = pass_pc.assign(pair = pass_pc['bc'] + '-' + pass_pc['translation'])
counts = pass_pc['pair'].value_counts()
surviving = []
for cutoff in range(cutoff_sample):
    surviving.append(counts.gt(cutoff).sum())

#%%
""" FIG S1B: PAIR CUTOFF FILTER """
fig, ax = plt.subplots(figsize = (6, 4))

ax.plot(range(pair_cutoff + 1), surviving[:pair_cutoff + 1], color='gray')
ax.plot(range(pair_cutoff, cutoff_sample), surviving[pair_cutoff:], color='#4F8DA7')

ax.set_ylim(0, surviving[-1]*10)
ax.set_xlabel('Pair cutoff')
ax.set_ylabel('Unique barcode-variant pairs')
ax.set_title('Surviving barcode-variant pairs vs. pair cutoff')
ax.plot(pair_cutoff, surviving[pair_cutoff], 'o', color = '#4F8DA7')
ax.text(
    pair_cutoff-1,
    surviving[pair_cutoff-1] + (.6 * surviving[-1]),
    f'Chosen cutoff: ≥{pair_cutoff} observations\n{surviving[pair_cutoff-1]:,} remaining unique pairs',
    fontsize=11,
    bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='black', linewidth=1)
)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS01_barcoding/pair_cutoff.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%% Finish filtering by pair cutoff
counts = pass_pc['pair'].value_counts()
pairs_above_cutoff = counts[counts >= pair_cutoff].index
pass_pc = pass_pc[pass_pc['pair'].isin(pairs_above_cutoff)]
pass_pc = pass_pc.drop('pair', axis=1)
bc_counts = pass_pc.groupby('bc').size().reset_index(name = 'bc_freq')
pass_pc = pass_pc.merge(bc_counts, on = 'bc', how = 'left')

### Map barcodes
mapped_barcodes = bcm.map_barcodes(pass_pc)

### View chastity values and filter
chastity_cutoff = 0.975

chastity_values = pd.DataFrame()
chastity_values['bc'] = mapped_barcodes['bc'].unique()
chastity_values['chastity_value'] = mapped_barcodes.groupby('bc').apply(
    lambda mmdata: mmdata['map_freq'].nlargest(2).iloc[0] /
    mmdata['map_freq'].nlargest(2).sum()
    ).tolist()
chastity_values['fraction_of_total'] = mapped_barcodes.groupby('bc').apply(
    lambda mmdata: mmdata['map_freq'].max() / 
    mmdata['map_freq'].sum()
    ).tolist()
imperfect_chastity_values = chastity_values[
    chastity_values['chastity_value'] != 1.0
    ]['chastity_value']
mapped_barcodes = pd.merge(mapped_barcodes, chastity_values, on='bc', how='left')

#%%
""" FIG S1C: CHASTITY VALUE CUTOFF FILTER """

fig, ax = plt.subplots(figsize = (6,4))
y, x, _ = ax.hist(imperfect_chastity_values, bins = ceil(len(imperfect_chastity_values)/50), color = '#4F8DA7')
perfectmaps_uniqct = len(chastity_values)-len(imperfect_chastity_values)
ax.text(x[0]+(x[len(x)-1]-x[0])/20, y.max()*.75, 
         f'{perfectmaps_uniqct:,} barcodes map\nperfectly (not shown)\n\nChosen cutoff: ≥{chastity_cutoff} chastity',
         fontsize = 10,
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='black', linewidth=1))
ax.set_xlabel('Chastity value')
ax.set_ylabel('Number of barcodes')
ax.set_title('Distribution of chastity values')
ax.set_yticks([0, 100, 200, 300, 400])

ax.vlines(chastity_cutoff, 0, y.max(), linestyle='dashed', color = 'black')
# Shade gray to the left of the cutoff
for left, right, height in zip(x[:-1], x[1:], y):
    if left <= chastity_cutoff:
        ax.fill_between([left, right], 0, height, color='gray')

plt.show()

fig.savefig(f'{base_dir}/results/supplementary/figS01_barcoding/chastity_cutoff.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S1D: NORFLOXACIN BARCODE-VARIANT CORRELATION SCATTER """

#Get data
norvar = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/final_functional_scores/var_nor0064-sel.csv', index_col = 'mut')
norvar = norvar.drop(norvar[norvar['counts_selected'] == 0].index) #Remove complete dropouts
norbc = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/final_functional_scores/bc_nor0064-sel.csv', index_col = 'mut')
norbc = norbc.drop(norbc[norbc['counts_selected'] == 0].index) #Remove complete dropouts
norvar, norbc = norvar.align(norbc, join='inner', axis=0) #Only shared variants

#Plot
fig, ax = plt.subplots(figsize = (4,4))
ax.scatter(norbc['f_score'], norvar['f_score'], s=20, alpha=0.5, color = '#4F8DA7')
ax.set_xlabel('Enrichment score (barcode-derived)')
ax.set_ylabel('Enrichment score (directly measured)')
ax.set_title('0.064µg/mL Norfloxacin')
ticks = [-10, -8, -6, -4, -2, 0, 2]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlim(-11, 3.5)
ax.set_ylim(-11, 3.5)
rsq = np.corrcoef(norbc['f_score'], norvar['f_score'])[0,1]**2
ax.annotate(f'$R^2$ = {rsq:.2f}', xy=(0.85, 0.35), xycoords='axes fraction',
            ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='black', lw=1))

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS01_barcoding/norfloxacin_corr.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S1E: ACRIFLAVINE BARCODE-VARIANT CORRELATION SCATTER """

#Get data
acrvar = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/final_functional_scores/var_acr20-sel.csv', index_col = 'mut')
acrvar = acrvar.drop(acrvar[acrvar['counts_selected'] == 0].index) #Remove complete dropouts
acrbc = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/final_functional_scores/bc_acr20-sel.csv', index_col = 'mut')
acrbc = acrbc.drop(acrbc[acrbc['counts_selected'] == 0].index) #Remove complete dropouts
acrvar, acrbc = acrvar.align(acrbc, join='inner', axis=0) #Only shared variants

#Plot
fig, ax = plt.subplots(figsize = (4,4))
ax.scatter(acrbc['f_score'], acrvar['f_score'], s=20, alpha=0.5, color = '#4F8DA7')
ax.set_xlabel('Enrichment score (barcode-derived)')
ax.set_ylabel('Enrichment score (directly measured)')
ax.set_title('20 µg/mL Acriflavine')
ticks = [-10, -8, -6, -4, -2, 0, 2, 4, 6]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlim(-11, 7)
ax.set_ylim(-11, 7)
rsq = np.corrcoef(acrbc['f_score'], acrvar['f_score'])[0,1]**2
ax.annotate(f'$R^2$ = {rsq:.2f}', xy=(0.85, 0.35), xycoords='axes fraction',
            ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='black', lw=1))

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS01_barcoding/acriflavine_corr.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

