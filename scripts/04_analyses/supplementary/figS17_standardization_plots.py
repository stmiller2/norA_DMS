#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 09:17:37 2023

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import matplotlib.pyplot as plt
import seaborn as sns
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp

#%%
""" FIG S17A: UNNORMALIZED VIOLINS """

#Load data
data = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = False, #Normalize set to False
                                fill_dropouts = False,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')
data.columns = ['Acriflavine', 'Ethidium', 'Norfloxacin', 'Ofloxacin', 'Pentamidine', 'Pipemidic acid', 'Puromycin', 'TPP']
ligcols = data.columns

#Plot
fig, ax = plt.subplots()
violin = sns.violinplot(data[ligcols].values, linecolor = 'black', linewidth = 1)
ax.set_xticklabels(ligcols, rotation = 30, ha = 'right', rotation_mode = 'anchor')
ax.set_ylabel('Raw functional score')
ax.set_title('Before standardization', fontsize = 14)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS17_standardization/unnormed.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S17B: NORMALIZED VIOLINS """

#Load data
data = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True, #Normalize set to True
                                fill_dropouts = False,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')
data.columns = ['Acriflavine', 'Ethidium', 'Norfloxacin', 'Ofloxacin', 'Pentamidine', 'Pipemidic acid', 'Puromycin', 'TPP']
ligcols = data.columns

#Plot
fig, ax = plt.subplots()
violin = sns.violinplot(data[ligcols].values, linecolor = 'black', linewidth = 1)
ax.set_xticklabels(ligcols, rotation = 30, ha = 'right', rotation_mode = 'anchor')
ax.set_ylabel('Standardized functional score')
ax.set_title('After standardization', fontsize = 14)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS17_standardization/normed.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)
