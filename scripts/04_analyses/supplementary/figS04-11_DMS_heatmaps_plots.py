#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 12:26:56 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np 
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import plot_utils as pl
from plot_utils import spec_cmap, spec_vmin, spec_vmax

#%%
""" FIG S04 - S11: DMS HEATMAPS FOR EACH SUBSTRATE """

#Load data
scores = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                  remove_low_quality=True, 
                                  subtract_h2o=True, 
                                  normalize_by_std=True,
                                  fill_dropouts=-100,  #arbitrary fill for dark heatmap color
                                  return_format='dictionary', 
                                  wt_policy='mut_names')

#Make all heatmaps
drugs = ['Acriflavine', 'Ethidium', 'Norfloxacin', 'Ofloxacin',
         'Pentamidine', 'Pipemidic acid', 'Puromycin', 'Tetraphenylphosphonium']
concs = [20, 6, 0.048, 0.008, 25, 1.56, 4.8, 9]

for drug, conc, key in zip(drugs, concs, scores.keys()):
    data = scores[key]
    data['mutation'] = data.index.str[-1]
    data['position'] = data.index.str[1:-1].astype(int)
    data.loc[data['f_hiq_norm'].isna(), 'SE_norm'] = np.nan
    data.loc[data['f_hiq_norm'] == -100, 'SE_norm'] = np.nan
    title =f'{drug} ({conc}ug/mL)'
    fig = pl.dms_heatmap(data = data, valcol = 'f_hiq_norm', secol = 'SE_norm', 
                         cmap = spec_cmap, vmin = spec_vmin, vmax = spec_vmax, title = title)
    fig.savefig(f'{base_dir}/results/supplementary/figS04-11_DMS_heatmaps/{drug.lower()}.png',
                bbox_inches = 'tight', dpi = 300, transparent = True)