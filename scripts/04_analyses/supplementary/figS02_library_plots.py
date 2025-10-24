#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 16:44:43 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%
""" FIG S2A: LIBRARY COVERAGE """

sizes = [7731, 30]
labels = ['Represented', 'Missing']
colors = ['#4F8DA7', '#EE842B']

fig, ax = plt.subplots(figsize=(4, 4))
wedges, _ = ax.pie(sizes, colors=[c for c in colors], startangle=310, wedgeprops={'alpha': .8})
ax.legend(wedges, labels, loc='upper right')
angle0 = (wedges[0].theta2 + wedges[0].theta1) / 2
x0 = 0.7 * np.cos(np.radians(angle0))
y0 = 0.7 * np.sin(np.radians(angle0))
ax.text(x0, y0, str(sizes[0]), ha='center', va='center', fontsize = 13)
angle1 = (wedges[1].theta2 + wedges[1].theta1) / 2
x1 = 1.15 * np.cos(np.radians(angle1))
y1 = 1.15 * np.sin(np.radians(angle1))
ax.text(x1, y1, str(sizes[1]), ha='center', va='center', fontsize = 13)
ax.set_title('Library is well-covered', fontsize = 13)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS02_library/pie.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S2B: INPUT LIBRARY DISTRIBUTION """

#Load variant counts for representative input library (inpA rep1)
data = pd.read_csv(f'{base_dir}/data/processed/specificity/variant_counts/var_counts_inpA-rep1.csv')['counts']
data = data / data.sum()

#Make plot
fig, ax = plt.subplots( figsize = (4,4))
ax.hist(data, bins=150, color= '#4F8DA7')
ax.set_xlabel('Frequency within population')
ax.set_ylabel('Number of variants')
ax.set_title('Input distribution is unbiased', pad = 15)
ax.set_xlim(-0.0001, 0.001)
ax.set_xticks([0, .4e-3, .8e-3])

#Annotate coefficient of variation and 90/10 skew
sk = data.quantile(.9)/data.quantile(.1)
cv = data.std() / data.mean()
ax.annotate(f'C.V. = {cv:.2f}\n90/10 skew = {sk:.2f}', 
            xy=(0.3, 0.7), xycoords='axes fraction', 
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8, ec='black'))

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS02_library/input.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S2C: BARCODES PER VARIANT """

#Get barcode associations for all subpools and compute num. barcodes per variant
bpv = pd.concat([pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/T1_bc_lookup.csv')['mut'].value_counts(),
                pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/T2_bc_lookup.csv')['mut'].value_counts(),
                pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/T3_bc_lookup.csv')['mut'].value_counts(),
                pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/T4_bc_lookup.csv')['mut'].value_counts(),
                pd.read_csv(f'{base_dir}/data/processed/barcoding/lookup_tables/T5_bc_lookup.csv')['mut'].value_counts()])

#Plot
fig, ax = plt.subplots(figsize=(4,4))
ax.hist(bpv, bins=50, color = '#4F8DA7')

# Annotate mean and median
mn = bpv.mean()
md = bpv.median()
ax.annotate(f'Mean = {mn:.1f}\nMedian = {round(md)}', 
            xy=(0.3, 0.7), xycoords='axes fraction', 
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8, ec='black'))

ax.set_title('Each variant is represented\nby many barcodes')
ax.set_xlabel('Barcode count')
ax.set_ylabel('Number of variants')

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS02_library/bpv.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)
