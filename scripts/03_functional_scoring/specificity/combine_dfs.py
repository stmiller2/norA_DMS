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

### Pickle specificity data
os.chdir(f'{base_dir}/data/processed/specificity')

files = ['acr20-sel.csv', 'eth6-sel.csv', 'nor0048-sel.csv', 'ofl0008-sel.csv',
         'pen25-sel.csv', 'ppa156-sel2.csv', 'pur48-sel.csv', 'tpp9-sel.csv',
         'h2o-sel.csv']

names = ['acr', 'eth', 'nor', 'ofl', 'pen', 'ppa', 'pur', 'tpp', 'h2o']

missing = pd.read_csv(f'{base_dir}/data/processed/specificity/final_functional_scores/missing_variants.csv')[['mut', 'f_score', 'SE']] #Variants not present in the original library

scores = {}
for i in range(len(files)):
    scores[names[i]] = pd.read_csv(f'{base_dir}/data/processed/specificity/final_functional_scores/{files[i]}')[['mut', 'counts_input', 'counts_selected', 'f_score', 'SE']]
    scores[names[i]] = pd.concat([scores[names[i]], missing], ignore_index=True)
    
    #Remove scores with error > half a standard deviation
    hlfstd = scores[names[i]]['f_score'].std()/2
    scores[names[i]]['f_hiq'] = [scores[names[i]]['f_score'][j] 
                                 if scores[names[i]]['SE'][j] < hlfstd 
                                 else np.nan 
                                 for j in range(len(scores[names[i]]))]
    
    #Normalize using dataset's own standard deviation
    scores[names[i]]['f_hiq_norm'] = scores[names[i]]['f_hiq']/scores[names[i]]['f_hiq'].std()
    scores[names[i]]['SE_norm'] = scores[names[i]]['SE']/scores[names[i]]['f_hiq'].std()
    
    #Set index to use mutations
    scores[names[i]].index = scores[names[i]]['mut']
    scores[names[i]] = scores[names[i]].drop(['mut'], axis = 1)
    
with open(f'{base_dir}/data/processed/specificity/specificity_clean.pk1', 'wb') as file:
    pickle.dump(scores, file)
