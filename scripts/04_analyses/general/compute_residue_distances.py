#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 15:37:42 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp

data = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'exclude') #Data only needed for index (all mutations)

""" Compute distance to the nearest biologically relevant residue """
cryoEM_inward_bindingsite = [12, 15, 16, 18, 19, 20, 23, 44, 45, 48, 51, 98, 105,
                             106,109, 132, 133, 135, 136, 137, 138, 140, 141, 211, 
                             214, 215, 218, 219, 222, 223, 252, 303, 306, 307, 
                             310, 313, 329, 332, 333, 336, 337, 340, 344] #Exposed to cavity in PyMol - see /results/general/9b3m_bindingsite.pse
cryoEM_outward_bindingsite = [12, 15, 16, 19, 20, 22, 47, 48, 51, 98, 105, 109, 
                              113, 117, 129, 132, 133, 136, 137, 140, 211, 214, 
                              218, 219, 222, 223, 255, 306, 307, 310, 311, 313, 
                              314, 317, 318, 325, 328, 329, 332, 333, 336, 337, 
                              340, 341, 344] #Exposed to cavity in PyMol - see /results/general/7lo8_bindingsite.pse
coupling_residues = [222, 307, 98] #From literature


dists = pd.DataFrame()
dists.index = data.index
dists['position'] = dists.index.str[1:-1].astype(str)

#Binding site residues
for pos in data['position'].unique():
    min_dist = float('inf')
    for res in np.unique(cryoEM_inward_bindingsite + cryoEM_outward_bindingsite): #For residues not resolved in cryoEM structure, find the nearest BS residue from either structure in the alphafold model
        res_dist = sp.min_residue_distance(f'{base_dir}/data/external/norA_AF.pdb', 'A', int(pos), int(res))
        min_dist = min(res_dist, min_dist)
    print(f'{pos} complete')
    dists.loc[dists['position'] == pos, 'alphafold_bindingsite'] = min_dist
for pos in data['position'].unique():
    min_dist = float('inf')
    for res in cryoEM_inward_bindingsite:
        res_dist = sp.min_residue_distance(f'{base_dir}/data/external/9b3m.pdb', 'A', int(pos), int(res))
        min_dist = min(res_dist, min_dist)
    print(f'{pos} complete')
    dists.loc[dists['position'] == pos, '9b3m_bindingsite'] = min_dist
for pos in data['position'].unique():
    min_dist = float('inf')
    for res in cryoEM_outward_bindingsite:
        res_dist = sp.min_residue_distance(f'{base_dir}/data/external/7lo8.pdb', 'Z', int(pos), int(res))
        min_dist = min(res_dist, min_dist)
    print(f'{pos} complete')
    dists.loc[dists['position'] == pos, '7lo8_bindingsite'] = min_dist

#Coupling residues
for pos in data['position'].unique():
    min_dist = float('inf')
    for res in coupling_residues:
        res_dist = sp.min_residue_distance(f'{base_dir}/data/external/norA_AF.pdb', 'A', int(pos), int(res))
        min_dist = min(res_dist, min_dist)
    print(f'{pos} complete')
    dists.loc[dists['position'] == pos, 'alphafold_coupling'] = min_dist
for pos in data['position'].unique():
    min_dist = float('inf')
    for res in coupling_residues:
        res_dist = sp.min_residue_distance(f'{base_dir}/data/external/9b3m.pdb', 'A', int(pos), int(res))
        min_dist = min(res_dist, min_dist)
    print(f'{pos} complete')
    dists.loc[dists['position'] == pos, '9b3m_coupling'] = min_dist
for pos in data['position'].unique():
    min_dist = float('inf')
    for res in coupling_residues:
        res_dist = sp.min_residue_distance(f'{base_dir}/data/external/7lo8.pdb', 'Z', int(pos), int(res))
        min_dist = min(res_dist, min_dist)
    print(f'{pos} complete')
    dists.loc[dists['position'] == pos, '7lo8_coupling'] = min_dist

dists.to_csv(f'{base_dir}/results/general/residue_distances.csv')