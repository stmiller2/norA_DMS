#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 21:01:53 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import scoring_utils as ngs

### Define input library associations
# Selections performed on different days have separate input libraries used 
# for pre-selection counts
sample_dict ={}
sample_dict['F'] = ['acr10pH7', 'nor0048pH7']
sample_dict['N'] = ['h2opH6', 'h2opH7', 'acr10pH6', 'nor0048pH6']

def get_inp(sel):
    for key in sample_dict.keys():
        for item in sample_dict[key]:
            if item in sel:
                return key

### Define barcode constant regions and lookup table
prebc = 'GTAACCGCCACG'
postbc = 'TGAGATCCGGCT'
lt = pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/bc_lookup.csv')

### Set sample_names list
sample_names = ['inpF-rep1', 'inpF-rep2',
                'inpN-rep1', 'inpN-rep2',
                'nor0048pH6-sel-rep1', 'nor0048pH6-sel-rep2',
                'nor0048pH7-sel-rep1', 'nor0048pH7-sel-rep2',
                'acr10pH6-sel-rep1', 'acr10pH6-sel-rep2',
                'acr10pH7-sel-rep1', 'acr10pH7-sel-rep2',
                'h2opH6-sel-rep1', 'h2opH6-sel-rep2',
                'h2opH7-sel2-rep1', 'h2opH7-sel2-rep2']

### Get raw barcode read counts
counts = {}
for name in sample_names:
    print(f'Processing sample {name}')
    #Load reads
    ngs.start_time('Importing csv...')
    raw = pd.read_csv(f'{base_dir}/data/processed/efficiency/good_reads/good_reads_{name}.csv', names = ['read'])
    ngs.end_time()
    
    #Extract barcodes
    raw['bc'] = ngs.extract_subsequence(raw['read'], prebc, postbc)
    print(f'\t{sum(raw["bc"].isna()):,} of {len(raw):,} did not contain both constants')
    
    #Count barcodes
    bc_counts = ngs.count_barcodes(raw['bc'], lt)
    counts[name] = ngs.sum_barcode_counts(bc_counts) #Sum barcode counts by variant and save to dict
    
    #Save csvs
    bc_counts.to_csv(f'{base_dir}/data/processed/efficiency/barcode_counts/bc_counts_{name}.csv', index = False)
    counts[name].to_csv(f'{base_dir}/data/processed/efficiency/variant_counts/var_counts_{name}.csv', index = False)
    
    print(f'Sample {name} complete', end = '\n\n')

#%%OPTIONAL: load variant count csvs from step 1
counts = {}
for name in sample_names:
    counts[name] = pd.read_csv(f'{base_dir}/data/processed/efficiency/variant_counts/var_counts_{name}.csv')
    
### Check replicate frequency correlations
conditions = sorted(list(set([item[:-5] for item in sample_names])))

for condition in conditions:
    fig = ngs.correlation_scatter(counts[f'{condition}-rep1']['counts']/counts[f'{condition}-rep1']['counts'].sum(),
                                  counts[f'{condition}-rep2']['counts']/counts[f'{condition}-rep2']['counts'].sum(),
                                  title = f'{condition} replicate correlation (frequency)')
    
    #Save replicate frequency correlation plots
    fig.savefig(f'{base_dir}/data/processed/efficiency/scoring_figs/frequency_correlations/frequency_corr_{condition}.png',
                dpi = 200, bbox_inches = 'tight')

### Get input distributions and plot
for input_lib in [con for con in conditions if 'inp' in con]:
    fig, ax = plt.subplots()
    r1 = counts[f'{input_lib}-rep1']['counts']/counts[f'{input_lib}-rep1']['counts'].sum()
    r2 = counts[f'{input_lib}-rep2']['counts']/counts[f'{input_lib}-rep1']['counts'].sum()
    bins = np.arange(0,.008, .00004)
    ax.hist(r1, bins = bins, alpha = 0.5, label = 'Rep. 1')
    ax.hist(r2, bins = bins, alpha = 0.5, label = 'Rep. 2')
    s1 = 3 * (np.mean(r1) - np.median(r1)) / np.std(r1)
    s2 = 3 * (np.mean(r2) - np.median(r2)) / np.std(r2)
    cv1 = np.std(r1) / np.mean(r1)
    cv2 = np.std(r2) / np.mean(r2)
    ax.annotate(f'Skewness\n  Rep1: {s1:.2f}\n  Rep2: {s2:.2f}\nC.V.\n  Rep1: {cv1:.2f}\n  Rep2: {cv2:.2f}', (0.8, 0.2), 
                  xycoords = 'axes fraction', 
                  bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8, ec='black'))
    ax.legend()
    ax.set_title(f'Input library {input_lib[-1]} distribution')
    ax.set_ylabel('Variant count')
    ax.set_xlabel('Proportion of total reads')
    plt.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.9,hspace=0.3,wspace=0.2)
    plt.show()
    
    #Save input distribution plots
    fig.savefig(f'{base_dir}/data/processed/efficiency/scoring_figs/input_distributions/input_dist_{input_lib}.png',
                dpi = 200, bbox_inches = 'tight')
    
    ax.set_xscale('log')
    fig.savefig(f'{base_dir}/data/processed/efficiency/scoring_figs/input_distributions/input_dist_log_{input_lib}.png',
                dpi = 200, bbox_inches = 'tight')

### Calculate functional scores
pc = 0.01 #Set pseudocount for dealing with complete dropouts (zero observations)

scores = {}
selections = [con for con in conditions if 'inp' not in con]

for selection in selections:
    input_lib = f'inp{get_inp(selection)}'
    scores[f'{selection}-rep1'] = ngs.variant_fn_scores(counts[f'{input_lib}-rep1'],
                                                        counts[f'{selection}-rep1'],
                                                        ref = 'WT',
                                                        pc = pc,
                                                        log = 'log2')
    scores[f'{selection}-rep2'] = ngs.variant_fn_scores(counts[f'{input_lib}-rep2'],
                                                        counts[f'{selection}-rep2'],
                                                        ref = 'WT',
                                                        pc = pc,
                                                        log = 'log2')
    
    #Save functional score csvs
    scores[f'{selection}-rep1'].to_csv(f'{base_dir}/data/processed/efficiency/replicate_functional_scores/{selection}-rep1.csv', index = False)
    scores[f'{selection}-rep2'].to_csv(f'{base_dir}/data/processed/efficiency/replicate_functional_scores/{selection}-rep2.csv', index = False)

#OPTIONAL: load functional score csvs from step 3
#conditions = sorted(list(set([item[:-5] for item in sample_names])))
#selections = [con for con in conditions if 'inp' not in con]
#scores = {}
#for selection in selections:
#    scores[f'{selection}-rep1'] = pd.read_csv(f'{base_dir}/data/processed/efficiency/replicate_functional_scores/{selection}-rep1.csv')
#    scores[f'{selection}-rep2'] = pd.read_csv(f'{base_dir}/data/processed/efficiency/replicate_functional_scores/{selection}-rep2.csv')

### Check replicate functional score correlations
for selection in selections:
    fig = ngs.correlation_scatter(scores[f'{selection}-rep1']['f_score'],
                                  scores[f'{selection}-rep2']['f_score'],
                                  title = f'{selection} replicate correlation (functional score)')
    #Save replicate score correlation plots
    fig.savefig(f'{base_dir}/data/processed/efficiency/scoring_figs/score_correlations/score_corr_{selection}.png',
                dpi = 200, bbox_inches = 'tight')

### Combine replicates with Enrich2 restricted maximum likelihood algorithm
final = {}
for selection in selections:
    #Collapse replicates
    final[selection] = ngs.collapse_replicates(scores[f'{selection}-rep1'],
                                               scores[f'{selection}-rep2'])
    
    #Save final scores csv
    final[selection].to_csv(f'{base_dir}/data/processed/efficiency/final_functional_scores/{selection}.csv', index = False)
    
    