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
import barcode_mapping_utils as bcm #For translate function

### Define barcode constant regions and lookup table
prebc = 'GTAACCGCCACG'
postbc = 'TGAGATCCGGCT'
lt = pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/T2_bc_lookup.csv')[['bc','mut']]

### Define reading frame constant regions and reflib
prerf = 'TTAATTATATGT'
postrf = 'CATCGTATGCCA'
ref_lib = pd.read_csv(f'{base_dir}/data/processed/barcoding/ref_libs/norA_T2_ref_lib.csv', header = None)

### Set sample_names list
sample_names = ['bc_inp_rep1', 'bc_inp_rep2',
                'var_inp_rep1', 'var_inp_rep2',
                'bc_acr20-sel_rep1', 'bc_acr20-sel_rep2',
                'var_acr20-sel_rep1', 'var_acr20-sel_rep2',
                'bc_nor0064-sel_rep1', 'bc_nor0064-sel_rep2',
                'var_nor0064-sel_rep1', 'var_nor0064-sel_rep2']

### Get raw read counts
counts = {}
for name in sample_names:
    print(f'Processing sample {name}')
    #Load reads
    ngs.start_time('Importing csv...')
    raw = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/good_reads/good_reads_T2-{name}.csv', names = ['read'])
    ngs.end_time()
    
    #Process barcode sequencing
    if 'bc' in name:
        #Extract barcodes
        raw['bc'] = ngs.extract_subsequence(raw['read'], prebc, postbc)
        print(f'\t{sum(raw["bc"].isna()):,} of {len(raw):,} did not contain both constants')
    
        #Count barcodes
        bc_counts = ngs.count_barcodes(raw['bc'], lt)
        counts[name] = ngs.sum_barcode_counts(bc_counts) #Sum barcode counts by variant and save to dict
    
        #Save csvs
        bc_counts.to_csv(f'{base_dir}/data/processed/empirical_barcode_validation/barcode_counts/bc_counts_{name}.csv', index = False)
        counts[name].to_csv(f'{base_dir}/data/processed/empirical_barcode_validation/variant_counts/var_counts_{name}.csv', index = False)
        
    #Process direct variant sequencing
    if 'var' in name:
        #Translate variants
        raw = bcm.translate_variants(raw, prerf, postrf)
        
        #Count variants
        counts[name] = ref_lib.merge(raw['translation'].value_counts().reset_index(), left_on = 1, right_on = 'translation', how = 'left').fillna({'count': 0})[[0, 'count']]
        counts[name].columns = ['mut', 'counts']
        
        #Save csv
        counts[name].to_csv(f'{base_dir}/data/processed/empirical_barcode_validation/variant_counts/var_counts_{name}.csv', index = False)
        
    print(f'Sample {name} complete', end = '\n\n')

#%%OPTIONAL: load variant count csvs from step 1
#counts = {}
#for name in sample_names:
#    counts[name] = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/variant_counts/var_counts_{name}.csv')
    
### Check replicate frequency correlations
conditions = sorted(list(set([item[:-5] for item in sample_names])))

for condition in conditions:
    fig = ngs.correlation_scatter(counts[f'{condition}_rep1']['counts']/counts[f'{condition}_rep1']['counts'].sum(),
                                  counts[f'{condition}_rep2']['counts']/counts[f'{condition}_rep2']['counts'].sum(),
                                  title = f'{condition} replicate correlation (frequency)')
    
    #Save replicate frequency correlation plots
    fig.savefig(f'{base_dir}/data/processed/empirical_barcode_validation/scoring_figs/frequency_correlations/frequency_corr_{condition}.png',
                dpi = 200, bbox_inches = 'tight')

### Get input distributions and plot
for input_lib in [con for con in conditions if 'inp' in con]:
    fig, ax = plt.subplots()
    r1 = counts[f'{input_lib}_rep1']['counts']/counts[f'{input_lib}_rep1']['counts'].sum()
    r2 = counts[f'{input_lib}_rep2']['counts']/counts[f'{input_lib}_rep2']['counts'].sum()
    bins = np.arange(0,.008, .00004)
    ax.hist(r1, bins = bins, alpha = 0.5, label = 'Rep. 1')
    ax.hist(r2, bins = bins, alpha = 0.5, label = 'Rep. 2')
    s1 = 3 * (np.mean(r1) - np.median(r1)) / np.std(r1)
    s2 = 3 * (np.mean(r2) - np.median(r2)) / np.std(r2)
    cv1 = np.std(r1) / np.mean(r1)
    cv2 = np.std(r2) / np.mean(r2)
    ax.annotate(f'Skewness\n  Rep1: {s1:.2f}\n  Rep2: {s2:.2f}\nC.V.\n  Rep1: {cv1:.2f}\n  Rep2: {cv2:.2f}', (0.7, 0.2), 
                  xycoords = 'axes fraction', 
                  bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8, ec='black'))
    ax.legend()
    ax.set_title('Input library distribution')
    ax.set_ylabel('Variant count')
    ax.set_xlabel('Proportion of total reads')
    plt.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.9,hspace=0.3,wspace=0.2)
    plt.show()
    
    #Save input distribution plots
    fig.savefig(f'{base_dir}/data/processed/empirical_barcode_validation/scoring_figs/input_distributions/input_dist_{input_lib}.png',
                dpi = 200, bbox_inches = 'tight')
    
    ax.set_xscale('log')
    fig.savefig(f'{base_dir}/data/processed/empirical_barcode_validation/scoring_figs/input_distributions/input_dist_log_{input_lib}.png',
                dpi = 200, bbox_inches = 'tight')

#%% Calculate functional scores
pc = 0.01 #Set pseudocount for dealing with complete dropouts (zero observations)

scores = {}
selections = [con for con in conditions if 'inp' not in con]

for selection in selections:
    input_lib = 'bc_inp' if 'bc' in selection else 'var_inp'
    scores[f'{selection}_rep1'] = ngs.variant_fn_scores(counts[f'{input_lib}_rep1'],
                                                        counts[f'{selection}_rep1'],
                                                        ref = 'WT',
                                                        pc = pc,
                                                        log = 'log2')
    scores[f'{selection}_rep2'] = ngs.variant_fn_scores(counts[f'{input_lib}_rep2'],
                                                        counts[f'{selection}_rep2'],
                                                        ref = 'WT',
                                                        pc = pc,
                                                        log = 'log2')
    
    #Save functional score csvs
    scores[f'{selection}_rep1'].to_csv(f'{base_dir}/data/processed/empirical_barcode_validation/replicate_functional_scores/{selection}_rep1.csv', index = False)
    scores[f'{selection}_rep2'].to_csv(f'{base_dir}/data/processed/empirical_barcode_validation/replicate_functional_scores/{selection}_rep2.csv', index = False)

#OPTIONAL: load functional score csvs from step 3
#conditions = sorted(list(set([item[:-5] for item in sample_names])))
#selections = [con for con in conditions if 'inp' not in con]
#scores = {}
#for selection in selections:
#    scores[f'{selection}-rep1'] = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/replicate_functional_scores/{selection}-rep1.csv')
#    scores[f'{selection}-rep2'] = pd.read_csv(f'{base_dir}/data/processed/empirical_barcode_validation/replicate_functional_scores/{selection}-rep2.csv')

### Check replicate functional score correlations
for selection in selections:
    fig = ngs.correlation_scatter(scores[f'{selection}_rep1']['f_score'],
                                  scores[f'{selection}_rep2']['f_score'],
                                  title = f'{selection} replicate correlation (functional score)')
    #Save replicate score correlation plots
    fig.savefig(f'{base_dir}/data/processed/empirical_barcode_validation/scoring_figs/score_correlations/score_corr_{selection}.png',
                dpi = 200, bbox_inches = 'tight')

### Combine replicates with Enrich2 restricted maximum likelihood algorithm
final = {}
for selection in selections:
    #Collapse replicates
    final[selection] = ngs.collapse_replicates(scores[f'{selection}_rep1'],
                                               scores[f'{selection}_rep2'])
    
    #Save final scores csv
    final[selection].to_csv(f'{base_dir}/data/processed/empirical_barcode_validation/final_functional_scores/{selection}.csv', index = False)
    
    