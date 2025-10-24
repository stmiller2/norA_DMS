#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 11:04:21 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp

#Load data
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=-3,
                                return_format='dataframe', 
                                wt_policy='exclude')
ligcols = data.columns

#Get cluster assignments
cluster_assignments = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv', index_col = 'mut')['cluster']
data['cluster'] = cluster_assignments

#Set up minimum activity mask
active_mask = (data[ligcols] > -1).sum(axis = 1)>0 #Only variants with F>-1
data = data[active_mask] #Remove variants that may be misfolded (below activity cutoff) for functional analysis

#Define clusters of interest
clusters_of_interest = {37: 'Cation impairment',
                        31: 'Puromycin tolerance',
                        32: 'High-TPSA tolerance',
                        12: 'TPP tolerance',
                        25: 'Pentamidine impairment',
                        2: 'Ethidium tolerance',
                        14: 'Universal enrichment'}

#Print cluster data for copy into Table S3
for cluster, phenotype in clusters_of_interest.items():
    cluster_data = data[data['cluster'] == cluster]
    print(f'Cluster #{cluster} ({phenotype}):')
    print(f'{len(cluster_data)} total members')
    print(', '.join(cluster_data.index))
