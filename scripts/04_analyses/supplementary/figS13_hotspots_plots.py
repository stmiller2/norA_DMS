#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:48:59 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import plot_utils as pl

#%%
""" FIG S13A: CONSERVATION """

#Create colorbar
cmap_struct = mpl.colors.LinearSegmentedColormap.from_list("", ['#a0275f', '#ffffff', '#1b7d82'])
fig, ax = plt.subplots(figsize = (3, .2))
norm = plt.Normalize(-1, 1)
sm = plt.cm.ScalarMappable(cmap = cmap_struct, norm = norm)
sm.set_array([])
cbar = plt.colorbar(sm, cax = ax, orientation = 'horizontal')
ax.text(-1, -2.5, 'Conserved', fontsize = 8, ha = 'center', style = 'italic')
ax.text(0, -2.5, 'Average', fontsize = 8, ha = 'center', style = 'italic')
ax.text(1, -2.5, 'Variable', fontsize = 8, ha = 'center', style = 'italic')
cbar.set_label('Consurf Score', labelpad = 15)
cbar.set_ticks([-1, 0, 1])
cbar.outline.set_linewidth(0.5)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS13_hotspots/conservation_colorbar.png', 
            bbox_inches = 'tight', dpi = 300, transparent = True)

#Get consurf data
consurf_scores = pd.read_excel(f'{base_dir}/data/external/consurf_scores.xlsx', header = None).iloc[29:417, :2].reset_index(drop = True) #Filter out consurf annotations, just keep data
consurf_scores.columns = ['residue', 'conservation']

#Define structural categories
support_helices = [(67, 84), (157, 176), (268, 288), (357, 382)]
core_helices = [(3, 33), (37, 65), (92, 119), (125, 149), (203, 232), (240, 264), (292, 320), (325, 353)]
sh_mask = consurf_scores['residue'].apply(lambda x: any(lo <= x <= hi for lo, hi in support_helices)) #Support helices
ch_mask = consurf_scores['residue'].apply(lambda x: any(lo <= x <= hi for lo, hi in core_helices)) #Core helices
ah_mask = sh_mask | ch_mask #All helices
al_mask = ~ah_mask #All loops

#Process consurf data
categories = ['All\nhelices', 'Core\nhelices', 'Support\nhelices', 'Loops']
scores = []
for mask in [ah_mask, ch_mask, sh_mask, al_mask]:
    scores.append(consurf_scores.loc[mask, 'conservation'].mean())

#Plot
fig, ax = plt.subplots(figsize = (3.5,6))
ax.bar(categories, scores, color=sm.to_rgba(scores), edgecolor='black')
ax.set_ylabel('ConSurf score')
ax.set_title('Sequence conservation')
ax.axhline(0, color='black', linewidth=0.8)
ax.text(-.65, -.55, 'More\nconserved', style = 'italic', ha = 'right', fontsize = 8)
ax.text(-.65, .65, 'More\nvariable', style = 'italic', ha = 'right', fontsize = 8)

#Adjust style
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(axis='x', length=0)
ax.axhline(0, color='black', linewidth=0.8)
ax.axvline(-0.55, color='black', linewidth=0.5)
ax.set_ylim(-0.55,0.75)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS13_hotspots/conservation_bar.png', 
            bbox_inches = 'tight', dpi = 300, transparent = True)


#%%
""" FIG S13B: ROSETTA & THERMOMPNN DDG PREDICTIONS BY CLUSTER & CONFORMATION """

# Load data
data = sp.load_specificity_data(data_path=f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality=True, 
                                subtract_h2o=True, 
                                normalize_by_std=True,
                                fill_dropouts=False,
                                return_format='dataframe', 
                                wt_policy='exclude')

#Get cluster assignments
cluster_assignments = pd.read_csv(f'{base_dir}/results/clustering/cluster_assignments.csv', index_col = 'mut')['cluster']
data['cluster'] = cluster_assignments

#Clean and reformat data
data.columns = ['acriflavine', 'ethidium', 'norfloxacin', 'ofloxacin', 'pentamidine', 'pipemidic acid', 'puromycin', 'TPP', 'cluster']
data = data.dropna()
data.index = data.index.str[1:]

#Get predicted stability data
stability_df = pd.read_csv(f'{base_dir}/data/external/computational_stability.csv',
                           index_col=0)
stability_df['cluster'] = data['cluster']

#Plot
fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

#Rosetta subplot
melted_df = stability_df.melt(id_vars=['cluster'], value_vars=['rosetta_7lo8', 'rosetta_9b3l', 'rosetta_9b3m'],
                              var_name='state', value_name='ddG')
melted_df['cluster'] = melted_df['cluster'].map({26: 'Univ. permitted', 7: 'Univ. disabling'})
melted_df['conformation'] = melted_df['state'].map({'rosetta_7lo8': 'Outward open',
                                                    'rosetta_9b3l': 'Occluded',
                                                    'rosetta_9b3m': 'Inward open'})
melted_df['ddG'] = melted_df['ddG'].clip(upper=25) #Clip Rosetta values at a max of 25 for better plotting (Rosetta produces some very high ddG outliers)

ax = axes[0]
sns.violinplot(data=melted_df, x='conformation', y='ddG', hue='cluster',
               split=True, inner='quartile', palette=['#1E86B3', '#dadcdf'], ax=ax)

#Annotate misfolding cutoffs for each conformation
ax.axhline(stability_df['rosetta_7lo8'].std(), xmin=.025, xmax=0.325, color='darkred', linestyle='--', linewidth=2)
ax.axhline(stability_df['rosetta_9b3l'].std(), xmin=.35, xmax=0.65, color='darkred', linestyle='--', linewidth=2)
ax.axhline(stability_df['rosetta_9b3m'].std(), xmin=.675, xmax=0.975, color='darkred', linestyle='--', linewidth=2)

#Label & style subplot
ax.set_ylabel('Predicted ΔΔG', labelpad = 10)
ax.text(-0.0725, 0.5, "clipped at 25 REU",
        fontsize=7, rotation=90, va='center', ha='center', transform=ax.transAxes)
ax.set_xlabel(None)
ax.tick_params(axis='x', length=0)
ax.axvline(-.5, color='black', linewidth=0.8)
ax.text(-.01, .1, 'More\nstable', style='italic', ha='right', va='top', fontsize=8, transform=ax.transAxes)
ax.text(-.01, .9, 'Less\nstable', style='italic', ha='right', va='bottom', fontsize=8, transform=ax.transAxes)
ax.set_ylim(-31, 32)
tix = [-20, -10, 0, 10, 20]
ax.set_yticks(tix, tix)
ax.set_title('Rosetta', ha='left', x=0.05, y=1)
for spine in ax.spines:
    if spine != 'bottom':
        ax.spines[spine].set_visible(False)

#ThermoMPNN subplot
melted_df = stability_df.melt(id_vars=['cluster'], value_vars=['thermo_7lo8', 'thermo_9b3l', 'thermo_9b3m'],
                              var_name='state', value_name='ddG')
melted_df['cluster'] = melted_df['cluster'].map({26: 'Univ. permitted', 7: 'Univ. disabling'})
melted_df['conformation'] = melted_df['state'].map({'thermo_7lo8': 'Outward open',
                                                    'thermo_9b3l': 'Occluded',
                                                    'thermo_9b3m': 'Inward open'})

ax = axes[1]
sns.violinplot(data=melted_df, x='conformation', y='ddG', hue='cluster',
               split=True, inner='quartile', palette=['#1E86B3', '#dadcdf'], ax=ax)

#Annotate misfolding cutoffs for each conformation
ax.axhline(stability_df['thermo_7lo8'].std(), xmin=.025, xmax=0.325, color='darkred', linestyle='--', linewidth=2)
ax.axhline(stability_df['thermo_9b3l'].std(), xmin=.35, xmax=0.65, color='darkred', linestyle='--', linewidth=2)
ax.axhline(stability_df['thermo_9b3m'].std(), xmin=.675, xmax=0.975, color='darkred', linestyle='--', linewidth=2)

#Label & style subplot
ax.set_ylabel('Predicted ΔΔG', labelpad = 18)
ax.set_xlabel(None)
ax.tick_params(axis='x', length=0)
ax.axvline(-.5, color='black', linewidth=0.8)
ax.text(-.01, .1, 'More\nstable', style='italic', ha='right', va='top', fontsize=8, transform=ax.transAxes)
ax.text(-.01, .9, 'Less\nstable', style='italic', ha='right', va='bottom', fontsize=8, transform=ax.transAxes)
tix = [-1, 0, 1, 2, 3]
ax.set_yticks(tix, tix)
ax.set_title('ThermoMPNN', ha='left', x=0.05, y=1)
for spine in ax.spines:
    if spine != 'bottom':
        ax.spines[spine].set_visible(False)

#Overall plot style
plt.suptitle('ΔΔG predictions by cluster & conformation', ha = 'left', x = 0.125, fontsize = 14)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right')
for ax in axes:
    ax.get_legend().remove()

plt.show()
plt.savefig(f'{base_dir}/results/supplementary/figS13_hotspots/rosetta_thermo_violins.png',
            bbox_inches='tight', dpi=300, transparent=True)



#%%
""" FIG S13C: FINDING FUNCTIONAL HOTSPOTS """

#Classify variants as folded or misfolded in each conformation
stability_df['oo_folded'] = (stability_df['rosetta_7lo8'] < stability_df['rosetta_7lo8'].std()) & (stability_df['thermo_7lo8'] < stability_df['thermo_7lo8'].std()) #Must meet cutoff in both Rosetta & ThermoMPNN predictions
stability_df['oc_folded'] = (stability_df['rosetta_9b3l'] < stability_df['rosetta_9b3l'].std()) & (stability_df['thermo_9b3l'] < stability_df['thermo_9b3l'].std()) #Must meet cutoff in both Rosetta & ThermoMPNN predictions
stability_df['io_folded'] = (stability_df['rosetta_9b3m'] < stability_df['rosetta_9b3m'].std()) & (stability_df['thermo_9b3m'] < stability_df['thermo_9b3m'].std()) #Must meet cutoff in both Rosetta & ThermoMPNN predictions

#Identify mutations that are folded & universally disabling
udis = stability_df[stability_df['cluster'] == 7] #Universally disabling mutations
udis_stable = udis[udis['oo_folded'] | udis['oc_folded'] | udis['io_folded']] #Universally disabling mutations that are folded in at least one conformation
udis_stable_positions = udis_stable.index.str[:-1].astype(int).value_counts() #Count number of mutations at each position fitting the above criteria

#Plot
fig, ax = plt.subplots(figsize = (4,5))
n, bins, patches = ax.hist(udis_stable_positions, bins=np.arange(0.5, 20.5, 1), edgecolor='black')
for patch, bin_edge in zip(patches, bins[:-1]):
    if bin_edge >= 9:
        patch.set_facecolor("#E44B3A")
    elif 4 <= bin_edge < 10:
        patch.set_facecolor("#EE842B")
    elif 2 <= bin_edge < 4:
        patch.set_facecolor("#FFD9A7")
    elif 0 <= bin_edge < 2:
        patch.set_facecolor("#747474")
ax.set_ylim(0, 60)
ax.set_xticks([0, 3, 6, 9, 12, 15, 18])
ax.set_xlabel('Num. substitutions causing loss of function')
ax.set_ylabel('Number of positions')
ax.set_title('Universally disabling mutants\npredicted to remain stable', ha = 'left', x = 0.075, y = .9)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.axvline(.1, color='black', linewidth=0.8)
ax.axhline(0.1, color='black', linewidth=0.8)

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS13_hotspots/hotspots_bar.png', 
            bbox_inches = 'tight', dpi = 300, transparent = True)

#Print PyMol commands for coloring structure figure
print('color white')
print(f'color 0xE44B3A, resi {"+".join(udis_stable_positions[udis_stable_positions>9].index.astype(str).tolist())}')
print(f'color 0xEE842B, resi {"+".join(udis_stable_positions[(udis_stable_positions<=9) & (udis_stable_positions>4)].index.astype(str).tolist())}')
print(f'color 0xFFD9A7, resi {"+".join(udis_stable_positions[(udis_stable_positions<=4) & (udis_stable_positions>2)].index.astype(str).tolist())}')
print(f'color 0x747474, resi {"+".join(udis_stable_positions[udis_stable_positions<=2].index.astype(str).tolist())}')

# Save csv of hotspots for manual annotation
wt = 'MNKQIFVLYFNIFLIFLGIGLVIPVLPVYLKDLGLTGSDLGLLVAAFALSQMIISPFGGTLADKLGKKLIICIGLILFSVSEFMFAVGHNFSVLMLSRVIGGMSAGMVMPGVTGLIADVSPSHQKAKNFGYMSAIINSGFILGPGIGGFMAEVSHRMPFYFAGALGILAFIMSVVLIHDPKKSTTSGFQKLEPQLLTKINWKVFITPAILTLVLAFGLSAFETLYSLYTSYKVNYSPKDISIAITGGGIFGALFQIYFFDKFMKYFSELTFIAWSLIYSVIVLVLLVIADGYWTIMVISFVVFIGFDMIRPAITNYFSNIAGDRQGFAGGLNSTFTSMGNFIGPLIAGALFDVHIEAPIYMAIGVSLAGVVIVLIEKQHRAKLKEQNM'
df = pd.DataFrame()
df['wt_res'] = [wt[i-1] for i in udis_stable_positions.index]
df['position'] = udis_stable_positions.index
df['num_disabling'] = udis_stable_positions.reset_index(drop = True)

df.to_csv(f'{base_dir}/results/general/hotspots.csv', index = None)

#%%
""" FIG S13D: SINGLE-RESIDUE HEATMAP FOR SELECTED HOTSPOT Y278 """

fig = pl.residue_heatmap(resi = 278, cbar = False, return_fig = True)
fig.savefig(f'{base_dir}/results/supplementary/figS13_hotspots/Tyr278.png', 
            bbox_inches = 'tight', dpi = 300, transparent = True)

#%%
""" FIG S13E: HIBIT MEMBRANE-LOCALIZED ABUNDANCE OF SELECT HOTSPOTS """

#Load data
df = pd.read_csv(f'{base_dir}/data/external/hibit_data.csv')
df = df[['Y278F','Y278I','L30N','T211W','P110S','P311M','Q255W', #Hotspot residues
         'WT_batch1', 'WT_batch3', 'WT_batch4', #WT
         'no_hibit', 'P144_G163del']] #Controls

#HiBiT experiments were performed in batches across five days. To normalize 
#day-to-day variability, wild type was included in every batch, and values are
#normalized to the wild type value recorded on that particular day. 
batches = {1: ["P144_G163del", "L30N", "P110S", "T211W", "Y278F"],
           3: ["Q255W", "Y278I", "I298D", "F303Y", "F306G"],
           4: ["M150E", "R310N", "P311M", "N332Q", "T336E", "H354Q", "E356P"]}

#Average technical replicates per biological replicate
bio_reps = ["rep1", "rep2", "rep3"]
bio_values = []
for i, rep in enumerate(bio_reps):
    tech_rows = df.iloc[i*3:(i+1)*3]
    bio_values.append(tech_rows.mean(axis=0))
bio_df = pd.DataFrame(bio_values, index=bio_reps)

#Subtract background (no hibit) per replicate
bio_df = bio_df.subtract(bio_df["no_hibit"], axis=0)
bio_df = bio_df.drop(columns="no_hibit")

#Normalize each replicate by its own batch's WT
for batch_num, variant_cols in batches.items():
    variant_cols = [col for col in variant_cols if col in bio_df.columns]
    wt_col = f"WT_batch{batch_num}"
    bio_df[[wt_col] + variant_cols] = (
        bio_df[[wt_col] + variant_cols].div(bio_df[wt_col], axis=0)
    )

#Keep only WT_batch1 as "WT" for plotting purposes
bio_df = bio_df.drop(columns=[c for c in bio_df.columns if c.startswith("WT_batch") and c != "WT_batch1"])
bio_df = bio_df.rename(columns={"WT_batch1": "WT"})

#Compute means and errors across normalized replicates
means = bio_df.mean(axis=0)
errs = bio_df.std(axis=0, ddof=1)

#Sort for plotting
means = means.sort_values()
errs = errs.loc[means.index]

labels = [r"$\bf{WT}$" if l == "WT" else
          r"$\bf{P144\_G163del}$" if l == "P144_G163del" 
          else l for l in 
          [label.replace("WT_batch1", "WT") for label in means.index]]

colors = ['#69B572' if l == r"$\bf{WT}$" else
          '#E1B23A' if l == r"$\bf{P144\_G163del}$" else 
          '#4F8DA7' for l in labels]

#Plot
fig =  plt.figure(figsize=(3, 4))
plt.bar(labels, means, yerr=errs, capsize=3, color=colors,
        error_kw=dict(elinewidth=1))
x_positions = np.arange(len(labels))
plt.axhline(1, linestyle='--', color='gray', linewidth=1)
plt.fill_between([-1,len(means)], 1/3, 3, color='gray', alpha=0.2, zorder=0)
plt.text(-.5, 2.25, '3-fold WT expression', fontsize = 9)
plt.fill_between([-1,len(means)], 1/2, 2, color='gray', alpha=0.2, zorder=0)
plt.text(-.5, 1.5, '2-fold WT expression', fontsize = 9)

#plot replicate scatter points
for i, label in enumerate(means.index):
    if label not in bio_df.columns:
        continue
    for rep in bio_df.index:
        y = bio_df.loc[rep, label]
        plt.scatter(i + np.random.uniform(-0.15, 0.15), y,
                    color='black', alpha=.3, s=15)

#Adjust style
plt.xlim(-1, len(means))
plt.yscale('log')
plt.ylim(np.log(1.028), np.log(1000))
plt.xticks(rotation=45, ha='right', rotation_mode = 'anchor', fontsize = 9)
plt.yticks(fontsize = 10)
plt.ylabel("Relative protein abundance (log)")
plt.tight_layout()

plt.show()
fig.savefig(f'{base_dir}/results/supplementary/figS13_hotspots//hibit_hotspots.png',
            bbox_inches = 'tight', dpi = 300, transparent = True)
