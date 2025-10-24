#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 21:01:53 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import pandas as pd; import numpy as np
import os
import pickle

### Pickle ΔpH sensitivity data
os.chdir(f'{base_dir}/data/processed/specificity')

files = ['acr10pH6-sel.csv', 'acr10pH7-sel.csv',
         'nor0048pH6-sel.csv', 'nor0048pH7-sel.csv',
         'h2opH6-sel.csv', 'h2opH7-sel2.csv']

names = ['acr_pH60', 'acr_pH70',
         'nor_pH60', 'nor_pH70', 
         'h2o_pH60', 'h2o_pH70']

missing = pd.read_csv(f'{base_dir}/data/processed/efficiency/final_functional_scores/missing_variants.csv')[['mut', 'f_score', 'SE']] #Variants not present in the original library

scores = {}
for i in range(len(files)):
    scores[names[i]] = pd.read_csv(f'{base_dir}/data/processed/efficiency/final_functional_scores/{files[i]}')[['mut', 'f_score', 'SE']]
    scores[names[i]] = pd.concat([scores[names[i]], missing], ignore_index=True)
    
    #Remove scores with error > half a standard deviation
    hlfstd = scores[names[i]]['f_score'].std()/2
    scores[names[i]]['f_hiq'] = [scores[names[i]]['f_score'][j] 
                                 if scores[names[i]]['SE'][j] < hlfstd 
                                 else np.nan 
                                 for j in range(len(scores[names[i]]))]
    
    #Set index to use mutations
    scores[names[i]].index = scores[names[i]]['mut']
    scores[names[i]] = scores[names[i]].drop(['mut'], axis = 1)
    
with open(f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1', 'wb') as file:
    pickle.dump(scores, file)

