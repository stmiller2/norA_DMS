#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 15:14:27 2024

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#%%
""" FIG 6A: CARTOON SURFACE RELATING ENERGY AVAILABLE FOR TRANSPORT TO ΔPH AND η """

show_mut = True
show_wt = True

#Set up toy data
pmf_range = (0,10) 
eta_range = (0, 1)
delta_pH, efficiency = np.meshgrid(np.linspace(pmf_range[0], pmf_range[1], 30), 
                                   np.linspace(eta_range[0], eta_range[1], 30))
available_energy = 3 * np.sqrt(delta_pH) * np.sqrt(efficiency) #sqrt for curvature

#Plot
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(delta_pH, efficiency, available_energy, 
                cmap='viridis', edgecolor = 'white', linewidth = .3)

#WT and mut lines
delta_pH_vals = np.linspace(pmf_range[0], pmf_range[1], 30)
if show_mut and show_wt:
    eta_vals = [0.65, 0.2]
    labels = ['WT', 'Mut']
    colors = [mpl.colors.to_hex(c) for c in mpl.colormaps['viridis'](eta_vals)]
elif show_wt and not show_mut:
    eta_vals = [0.65]
    labels = ['WT']
    colors = [mpl.colors.to_hex(mpl.colormaps['viridis'](eta_vals))]

if show_wt or show_mut:
    for eta, color, label in zip(eta_vals, colors, labels):
        energy_vals = 3 * np.sqrt(delta_pH_vals) * np.sqrt(eta)
        ax.plot(delta_pH_vals, [eta] * len(delta_pH_vals), energy_vals, 
                color='black', linewidth=1.5, label=label, zorder = 3)
        ax.plot([10, 10], [eta, 1], [3*np.sqrt(10)*np.sqrt(eta), 3*np.sqrt(10)*np.sqrt(eta)],
                color='black', linestyle='dotted', linewidth=1.5)
        ax.text(8, eta, 3*np.sqrt(10)*np.sqrt(eta), label, fontsize=12, weight='bold')

#Label and style
ax.set_xlabel('PMF (ΔG$_{proton}$)', fontsize = 12, weight = 'bold')
ax.set_ylabel('Efficiency (η)', fontsize = 12, weight = 'bold')
ax.set_zlabel('Energy available\nfor transport\n', fontsize = 12, weight = 'bold')
ax.set_xticks([0,2,4,6,8,10],['Lower', '', '', '', '', 'Higher'])
ax.set_yticks([0,.2,.4,.6,.8,1.0],['0', '', '', '', '', '1'])
ax.set_zticks([0,2,4,6,8,10],[])
ax.set_box_aspect((4, 4, 3), zoom=0.8)

plt.show()
fig.savefig(f'{base_dir}/results/fig6/pmf_eta_budget.png',
            dpi = 300, bbox_inches = 'tight', transparent = True)

#%%
""" FIG 6B: BAR PLOT ILLUSTRATING INEFFICIENT VARIANTS BECOME LESS PROMISCUOUS """

#Set up toy data
labs = ['ΔG$_{proton}$', 'WT', 'Mut', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
bars = np.array([11, 8, 5, 7.5, 3, 2, 7, 4, 6.5, 4.5, 6])
bars = np.array([11, 8, 5, 2, 3, 4, 4.5, 6, 6.5, 7, 7.5])
xpos = np.array([0, 2.5, 3.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5])
sub_cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', ['#ffffff', '#ec684b'])
colors = [mpl.colors.to_hex(c) for c in mpl.colormaps['viridis'](bars/12)[:3]] + [mpl.colors.to_hex(c) for c in sub_cmap(bars/8)[3:]]
w = .75

#Plot
fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(xpos, bars, color = colors, width = w)
ax.hlines(y = bars[1], xmin = xpos[1]-(w/2), xmax = xpos[-1]+(w/2), 
          color = colors[1], linestyles='--', linewidth=1)
ax.hlines(y = bars[2], xmin = xpos[2]-(w/2), xmax = xpos[-1]+(w/2), 
          color = colors[2], linestyles='--', linewidth=1)
ax.fill_between((xpos[1]+(w/2), xpos[-1]+(w/2)), 0, bars[1], color = colors[1], alpha = 0.1, zorder = 0)
ax.fill_between((xpos[2]+(w/2), xpos[-1]+(w/2)), 0, bars[2], color = colors[2], alpha = 0.2, zorder = 0)

#Label and style
ax.text(xpos[[1,2]].mean(), 9.3, 'Energy available\nfor transport', 
        fontsize=12, weight='bold', ha='center')
ax.text(xpos[[1,2]].mean(), 8.6, '$ΔG_{available} = f(η, ΔG_{proton})$', 
        fontsize=10, style='italic', ha='center')
ax.text(xpos[3:].mean(), 9.3, 'Energetic cost\nof transport', 
        fontsize=12, weight='bold', ha='center')
ax.text(xpos[3:].mean(), 8.6, 'F(substrate identity)', 
        fontsize=10, style='italic', ha='center')
ax.text(xpos[[2,3]].mean(), 6.25, 'Substrates of\n$WT$ only',
        fontsize = 9, ha = 'center')
ax.text(xpos[[2,3]].mean(), 3, 'Substrates of\n$WT$ and $Mut$',
        fontsize = 9, ha = 'center')
ax.text(xpos[1],-1.5, 'η$_{WT}$', fontsize=10, style='italic', ha='center')
ax.text(xpos[[1,2]].mean(),-1.5, '>', fontsize=10, style='italic', ha='center')
ax.text(xpos[2],-1.5, 'η$_{mut}$', fontsize=10, style='italic', ha='center')
ax.set_xticks(xpos, labs, fontsize = 12)
ax.set_yticks([],[])
ax.set_ylim(0, bars.max()*1.15)
ax.set_ylabel('Energy', fontsize = 12, weight = 'bold')
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
ax.spines.left.set_bounds((0, bars.max()))
    
plt.show()
fig.savefig(f'{base_dir}/results/fig6/energy_bars.png',
            dpi = 300, bbox_inches = 'tight', transparent = True)

