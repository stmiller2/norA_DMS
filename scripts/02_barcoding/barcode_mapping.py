# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 10:18:24 2023

@author: Silas
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import barcode_mapping_utils as bcm

subpools = ['T1', 'T2', 'T3', 'T4', 'T5']
prebc = 'GTAACCGCCACG' #Sequence immediately preceding the 20N barcode
postbc = 'TGAGATCCGGCT' #Sequence immediately following the 20N barcode
prerfs = ['AGGAAAACTAGT', 'TTAATTATATGT', 'GCAGAAGTTTCA', 'TATACATCGTAT', 'ATGATAAGACCT'] #Sequence immediately preceding the reading frame of each subpool
postrfs = ['ATAGGGTTAATT', 'CATCGTATGCCA', 'AAAGTAAATTAT', 'GCTATTACAAAT', 'TAGCTCGCTTGT'] #Sequence immediately following the reading frame of each subpool

for subpool, prerf, postrf in zip(subpools, prerfs, postrfs):
    ref_lib = pd.read_csv(f'{base_dir}/data/processed/barcoding/ref_libs/norA_{subpool}_ref_lib.csv', header = None)
    seq_df = pd.read_csv(f'{base_dir}/data/processed/barcoding/good_reads/norA_{subpool}_good_reads.csv', names = ['read'])

    ### Extract barcodes
    seq_df = bcm.extract_barcodes(seq_df, prebc, postbc)

    ### Translate variants
    seq_df = bcm.translate_variants(seq_df, prerf, postrf)

    ### Group barcodes by levenshtein distance = 1
    bc_groups = bcm.group_barcodes(seq_df)

    ### Collapse barcode groups
    seq_df = bcm.collapse_groups(seq_df, bc_groups)

    ### Sample read cutoffs and filter
    pass_pc = bcm.sample_cutoffs_and_filter(seq_df, 
                                            cutoff_sample = 25, 
                                            filepath = f'{base_dir}/data/processed/barcoding/mapping_figs/{subpool}')

    ### Map barcodes
    mapped_barcodes = bcm.map_barcodes(pass_pc)

    ### View chastity values and filter
    mapped_barcodes = bcm.sample_chastity_and_filter(mapped_barcodes, 
                                                     filepath = f'{base_dir}/data/processed/barcoding/mapping_figs/{subpool}')

    ### Name variants
    mapped_barcodes = bcm.id_mutants(mapped_barcodes, ref_lib).dropna().reset_index(drop = True)

    ### Save barcode mapping dataframes
    mapped_barcodes.to_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/{subpool}_bc_lookup.csv', index = False)

### Combine subpool barcode lookup tables into master lookup table with all bcs
lookup_tables = {}
for subpool in subpools:
    lookup_tables[subpool] = pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/{subpool}_bc_lookup.csv')
    
combined_df = pd.concat(lookup_tables.values(), ignore_index=True)[['bc', 'mut']]
combined_df.to_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/bc_lookup.csv')