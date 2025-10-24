#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:26:33 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from time import time
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import specificity_utils as sp
import efficiency_utils as ef


aa_color_dict = {'R': '#E44B3A', 'H': '#E44B3A', 'K': '#E44B3A', 'D': '#4F8DA7', 'E': '#4F8DA7',
                 'S': '#EE842B', 'T': '#EE842B', 'N': '#EE842B', 'Q': '#EE842B', 'C': '#EE842B',
                 'G': '#A87BB7', 'P': '#A87BB7', 'A': '#747474', 'V': '#747474', 'I': '#747474',
                 'L': '#747474', 'M': '#747474', 'F': '#69B572', 'Y': '#69B572', 'W': '#69B572', 
                 '*': 'lightgray'}
#https://coolors.co/4f8da7-2bb1a1-86b26e-e1b23a-ee842b-e45c3a-c66c79-a87bb7

sub_colors = ['#4f8da7', '#2bb1a1', '#86B26E', '#e1b23a', '#ee842b', '#e45c3a', '#C66C79', '#A87BB7']

spec = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                remove_low_quality = True,
                                subtract_h2o = True,
                                normalize_by_std = True,
                                fill_dropouts = -3,
                                return_format = 'dataframe',
                                wt_policy = 'exclude')
effc = ef.load_efficiency_data(data_path = f'{base_dir}/data/processed/efficiency/efficiency_clean.pk1',
                               remove_low_quality=True,
                               subtract_h2o=True,
                               normalize_by_std=True,
                               return_format='scores_dataframe',
                               wt_policy = 'exclude')

""" Specificity custom colormap """
spec_vmin = spec.min().min(); spec_vmax = spec.max().max()
limit = max(abs(spec_vmin), abs(spec_vmax))
spec_norm = mpl.colors.Normalize(-limit,limit)
colors = [[spec_norm(-limit), "#1E86B3"],
          [spec_norm(0), '#dadcdf'],
          [spec_norm(limit), "#E82121"]]
spec_cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)
spec_cmap.set_under('#155E7D')

""" Efficiency custom colormaps """
#Norfloxacin
nor_effc_vmin = min(effc['nor_dFpH']); nor_effc_vmax = max(effc['nor_dFpH'])
limit = max(abs(nor_effc_vmin*1.1), abs(nor_effc_vmax*1.1))
effc_norm = mpl.colors.Normalize(-limit, limit)
sorted_indices = effc['nor_dFpH'].sort_values(key = abs).index

colors = [[effc_norm(-limit), "#573280"],
          [effc_norm(-1*effc['nor_dFpH'].std()), '#dadcdf'],
          [effc_norm(1*effc['nor_dFpH'].std()), '#dadcdf'],
          [effc_norm(limit), "#F58300"]]
nor_effc_cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)

#Acriflavine
acr_effc_vmin = min(effc['acr_dFpH']); acr_effc_vmax = max(effc['acr_dFpH'])*1.1
limit = max(abs(acr_effc_vmin*1.1), abs(acr_effc_vmax*1.1))
norm = mpl.colors.Normalize(-limit, limit)
sorted_indices = effc['acr_dFpH'].sort_values(key = abs).index

colors = [[norm(-limit), "#573280"],
          [norm(-1*effc['acr_dFpH'].std()), '#dadcdf'],
          [norm(1*effc['acr_dFpH'].std()), '#dadcdf'],
          [norm(limit), "#F58300"]]
acr_effc_cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)

#%%
def start_time(message = 'Starting task...'):
    print(message, end = '')
    global STIME
    STIME = time()
    
def end_time():
    ETIME = time()
    runtime = ETIME-STIME
    if runtime > 60:
        print(f' Done, {round(runtime/60, 2)} minutes', end = '\n')
    else:
        print(f' Done, {round(runtime, 2)} seconds', end = '\n')

def stack_plot(cluster_data, 
               wt_label_count, 
               mut_label_count, 
               title = None,
               bar_width = 0.8, label_size = 12, text_x = -0.25):
    wt = pd.Series([mut[0] for mut in cluster_data.index]).value_counts()
    mut = pd.Series([mut[-1] if '_' not in mut else '*' for mut in cluster_data.index]).value_counts()
    text_v_adjust = len(cluster_data) * 0.055
    
    fig, ax = plt.subplots(figsize = (2,4))
    h = len(cluster_data)
    for i, (aa, f) in enumerate(wt.items()):
        ax.bar(0, h, color = aa_color_dict[aa], edgecolor = 'white', width = bar_width)
        if i < wt_label_count:
            ax.text(text_x, h-text_v_adjust, aa, ha = 'center', color = 'white', weight = 'bold', fontsize = label_size)
        h = h - f
    h = len(cluster_data)
    for i, (aa, f) in enumerate(mut.items()):
        ax.bar(1, h, color = aa_color_dict[aa], edgecolor = 'white', width = bar_width)
        if i < mut_label_count:
            ax.text(1+text_x, h-text_v_adjust, aa, ha = 'center', color = 'white', weight = 'bold', fontsize = label_size)
        h = h - f
    ax.set_xticks([0,1])
    ax.set_xticklabels(['WT', 'Mut'], fontsize = 12)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_aspect(2.5/ax.get_data_ratio())
    if title:
        ax.set_title(title)
    return fig

def cluster_heatmap_small(data, cmap, vmin, vmax, labels = False, figsize = (4, .5)):
    fig, ax = plt.subplots(figsize = figsize)
    sns.heatmap(data,
                cmap = cmap,
                center = 0,
                vmin = vmin, vmax = vmax,
                cbar = False)
    ax.set_xticks([])
    ax.set_yticks([])
    for x in np.arange(1, 8, 1):
        ax.axvline(x, color = 'black', linewidth = 0.1)
    for x, lig in enumerate(['acr','eth','nor','ofl','pen','ppa','pur','tpp']):
        ax.text((x/8)+.06, -.4, lig, ha = 'center', transform = ax.transAxes)
    if labels:
        ax.set_ylabel('Mutations')
        ax.set_xlabel('Substrates')
    else:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_visible(True)
    return fig

def dms_heatmap(data, vmin, vmax, title, mutcol = 'mutation', poscol = 'position',
                valcol = 'f_norm_hiq', secol = 'SE_hiq', error_limit = 1.5, 
                num_sections = 3, aa_order = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W','_'],
                wt = 'MNKQIFVLYFNIFLIFLGIGLVIPVLPVYLKDLGLTGSDLGLLVAAFALSQMIISPFGGTLADKLGKKLIICIGLILFSVSEFMFAVGHNFSVLMLSRVIGGMSAGMVMPGVTGLIADVSPSHQKAKNFGYMSAIINSGFILGPGIGGFMAEVSHRMPFYFAGALGILAFIMSVVLIHDPKKSTTSGFQKLEPQLLTKINWKVFITPAILTLVLAFGLSAFETLYSLYTSYKVNYSPKDISIAITGGGIFGALFQIYFFDKFMKYFSELTFIAWSLIYSVIVLVLLVIADGYWTIMVISFVVFIGFDMIRPAITNYFSNIAGDRQGFAGGLNSTFTSMGNFIGPLIAGALFDVHIEAPIYMAIGVSLAGVVIVLIEKQHRAKLKEQNM',
                cmap = 'drug', cbar_label = None):
    """
    Generate a wrapped heatmap for displaying deep mutational scanning data

    Parameters
    ----------
    data : pd.DataFrame()
        Pandas dataframe with your data. At a minimum, it must contain columns
        with the mutant ID (one-letter AA code), position (integer), and 
        functional score or other value you'd like to plot.
    vmin : float
        The minimum value for scaling the colorbar. Set this to the overall
        minimum for the plots you want to directly compare.
    vmax : float
        The maximum value for scaling the colorbar. Set this to the overall
        maximum for the plots you want to directly compare.
    title : string
        The title of the plot; e.g., "Acriflavine (20ug/mL)".
    mutcol : string, optional
        The name of the column in the 'data' dataframe in which mutant IDs are
        stored. The default is 'mutation'.
    poscol : string, optional
        The name of the column in the 'data' dataframe in which mutation
        positions are stored. The default is 'position'.
    valcol : string, optional
        The name of the column in the 'data' dataframe in which the value you'd
        like to plot is stored. The default is 'f_norm_hiq'.
    secol : string, optional
        The name of the column in the 'data' dataframe in which the error
        values are stored. Set to "None" if you wish to not display any erorr
        bars / don't have that data. The default is 'SE_hiq'. 
    error_limit : float, optional
        The maximum error value for which datapoints should still be displayed.
        Any datapoint with an error higher than this will be masked. The 
        default is 1.5.
    num_sections : integer, optional
        The number of sections to split the heatmap into, for easier fitting of
        large datasets in one figure. The default is 3.
    aa_order : list, optional
        The order in which to display the amino acids. Changing the default is
        not recommended; The default is ordered by chemical properties.
    wt : string, optional
        The wild type sequence of your protein of interest. The default is the
        sequence of the multidrug transporter NorA.
    cmap : matplotlib colormap, optional
        Colormap to use for the heatmap. The default is drug, which uses custom
        colors developed for NorA DMS paper drug screens. pH uses purple/orange
        colormap used in NorA DMS paper pH screens. Setting to None will use
        default coolwarm with missing values set to white. Can also set this to
        be a cmap of your choosing.

    Returns
    -------
    fig : matplotlib.pyplot figure
        A heatmap.

    """
    start_time(f'Plotting heatmap for {title.lower()}...')
    if cmap == None:
        cmap = mpl.colormaps.get_cmap('coolwarm').copy()
        cmap.set_extremes(bad = 'white')
    if cmap == 'drug':
        limit = max(abs(vmin), abs(vmax))
        norm = mpl.colors.Normalize(-limit,limit)
        colors = [[norm(-limit), "#1E86B3"],
                  [norm(0), '#dadcdf'],
                  [norm(limit), "#E82121"]]
        cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)
        cmap.set_under('#155E7D')
    if cmap == 'pH':
        limit = max(abs(vmin), abs(vmax))
        norm = mpl.colors.Normalize(-limit,limit)
        sigma = data[valcol].std()
        colors = [[norm(-limit), "#573280"],
                  [norm(-1*sigma), '#dadcdf'],
                  [norm(1*sigma), '#dadcdf'],
                  [norm(limit), "#F58300"]]
        cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)
    #Make pivot table
    pivot = data.pivot(index = mutcol, columns = poscol, values = valcol).reindex(aa_order, axis = 0)
    if secol != None:
        pivot_se = data.pivot(index = mutcol, columns = poscol, values = secol).reindex(aa_order, axis = 0)
    #Calculate number of rows per section
    rows_per_section = np.shape(pivot)[1] // num_sections
    if secol != None:
        #Mask data with high error
        mask = pivot_se > error_limit
        pivot_masked = pivot.mask(mask)
        pivot_se_masked = pivot_se.mask(mask)
    else:
        pivot_masked = pivot
    #Create subplots for each section
    fig, axes = plt.subplots(num_sections, 1, figsize = (15, 2.5*num_sections))
    cbar_ax = fig.add_axes([.81, .18, .015, .65])
    #Plot each section
    for i, ax in enumerate(axes):
        start_row = i * rows_per_section
        end_row = start_row + rows_per_section
        if i == num_sections-1:
            end_row = np.shape(pivot_masked)[1]
        section_data = pivot_masked.iloc[:, start_row:end_row]
        section_data = section_data.replace(-np.inf, -999)
        if secol != None:
            section_se = pivot_se_masked.iloc[:, start_row:end_row]
        #Plot heatmap
        sns.heatmap(section_data, square = True, cmap = cmap, vmin = vmin, vmax = vmax,
                    cbar = (i == 0), cbar_ax = None if i else cbar_ax, 
                    cbar_kws = None if i else {'label': cbar_label},
                    center = 0, ax = ax)
        #Get WT positions
        section_wt = wt[start_row:end_row]
        wt_coords = []
        for j, aa in enumerate(section_wt):
            ypos = aa_order.index(aa)
            xpos = j
            wt_coords.append((xpos, ypos))
        #Highlight WT cells
        for coords in wt_coords:
            dot = plt.Circle((coords[0]+0.5, coords[1]+0.5), 0.2, color = 'gray')
            ax.add_artist(dot)
        if secol != None:
            #Add error bars
            for k in range(section_data.shape[0]):
                for l in range(section_data.shape[1]):
                    x = l + 0.5 #x-coordinate of cell center
                    y = k + 0.5 #y-coordinate of cell center
                    length = section_se.iloc[k, l] / error_limit #length of line as fraction of cell width
                    ax.plot([x - length / 2, x + length / 2], [y - length / 2, y + length / 2],
                            color = 'black', linewidth = 0.5, transform = ax.transData,
                            clip_on = False, alpha = 0.5, zorder = 3, linestyle = 'solid')
        #Set up labels
        x_tick_labels = [str(m+1+(i*rows_per_section)) if (m+1+(i*rows_per_section)) % 5 == 0 else '' for m in range(len(section_wt))]
        begin_ticks = [(m+1+(i*rows_per_section)) for m in range(len(section_wt)) if (m+1+(i*rows_per_section)) % 5 == 0][0] - (i*rows_per_section) - 1
        ax.set_xticks(range(begin_ticks, len(section_wt), 5)) #Add an x tick every 5 cells
        ax.set_xticklabels([]) #Remove x tick labels
        ax.set_xlabel('') #Remove x label
        ax.tick_params(axis='x', which='both', length=8, direction='out', top=True, bottom=False)
        for n, label in enumerate(x_tick_labels):
            ax.text(n+0.25, -0.2, label, ha='left', fontsize=8, alpha=0.6) #Add custom tick labels as text
        ax.set_yticks(np.arange(0.5, 21.5, 1))
        ax.set_yticklabels(aa_order[0:-1]+['*'], rotation = 0, fontsize = 6)
        ax.set_ylabel('')
        if i == 0:
            ax.set_title(title, fontsize = 15, y = 1.075)
        if cbar_label:
            cbar_ax.yaxis.set_label_coords(3, 0.5) #Sets position manually so figsize is the same regardless of cbar tick label length
    #Adjust spacing
    plt.subplots_adjust(top=0.925,
                        bottom=0.125,
                        left=0.05,
                        right=0.8,
                        hspace=0.2,
                        wspace=0.2)
    end_time()
    plt.show()
    return fig

def residue_heatmap(resi, cbar = True, return_fig = False):
    plt.rcParams['figure.dpi'] = 300
    data_dict = sp.load_specificity_data(data_path = f'{base_dir}/data/processed/specificity/specificity_clean.pk1',
                                         remove_low_quality = True,
                                         subtract_h2o = True,
                                         normalize_by_std = True,
                                         fill_dropouts = -999, #Arbitrary fill for dark heatmap color
                                         return_format = 'dictionary',
                                         wt_policy = 'mut_names')
    
    data = pd.DataFrame()
    for key, df in data_dict.items():
        data[key] = df['f_hiq_norm']
        data[f'{key}_se'] = df['SE_norm']
    ligcols = [col for col in data.columns if '_se' not in col]
    errcols = [col for col in data.columns if '_se' in col]
    vmin = data[ligcols].replace(-999, np.nan).min().min(); vmax = data[ligcols].max().max()
    limit = max(abs(vmin), abs(vmax))
    norm = mpl.colors.Normalize(-limit,limit)
    colors = [[norm(-limit), "#1E86B3"],
              [norm(0), '#dadcdf'],
              [norm(limit), "#E82121"]]
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)
    cmap.set_under('#155E7D')
    
    error_limit = 0.5
    aa_order = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W','_']
    aa_names = {'R': 'Arg','H': 'His','K': 'Lys','D': 'Asp','E': 'Glu','S': 'Ser','T': 'Thr','N': 'Asn','Q': 'Gln','C': 'Cys','G': 'Gly','P': 'Pro','A': 'Ala','V': 'Val','I': 'Ile','L': 'Leu','M': 'Met','F': 'Phe','Y': 'Tyr','W': 'Trp'}
    
    data['position'] = data.index.str[1:-1].astype(int)
    data['mut'] = data.index
    
    resi_data = data[data['position'] == resi]
    resi_data.index = resi_data['mut'].str[-1]
    resi_data = resi_data.reindex(aa_order, axis = 0)
    
    fig, ax = plt.subplots()
    if cbar:
        cbar_ax = fig.add_axes([.65, .18, .03, .65])
    
    if cbar:
        sns.heatmap(resi_data[ligcols], cmap=cmap, center = 0, 
                    vmax = vmax, vmin = vmin,
                    ax = ax, cbar_ax = cbar_ax, square = True,
                    cbar_kws={'label': 'Standardized functional score'})
        cbar_ax.figure.axes[-1].yaxis.label.set_rotation(270)
        cbar_ax.figure.axes[-1].yaxis.label.set_verticalalignment('bottom')
    else:
        sns.heatmap(resi_data[ligcols], cmap=cmap, center = 0, 
                    vmax = vmax, vmin = vmin,
                    ax = ax, cbar = None, square = True)
        
    #Add error bars
    for k in range(len(aa_order)):
        for l in range(len(ligcols)):
            x = l + 0.5 #x-coordinate of cell center
            y = k + 0.5 #y-coordinate of cell center
            length = resi_data[errcols].iloc[k, l] / error_limit #length of line as fraction of cell width
            if length < 1:
                ax.plot([x - length / 2, x + length / 2], [y - length / 2, y + length / 2],
                        color = 'black', linewidth = 0.5, transform = ax.transData,
                        clip_on = False, alpha = 0.5, zorder = 3, linestyle = 'solid')
    
    #Labels & styling
    ax.set_xticklabels(ligcols, rotation=45, rotation_mode = 'anchor',
                       ha='right', fontsize = 8)
    ax.set_yticks(np.arange(0.5, 21.5, 1))
    ax.set_yticklabels(aa_order[0:-1]+['*'], rotation = 0, fontsize = 8)
    ax.set_ylabel(None)
    title = aa_names[resi_data["mut"].iloc[-1][0]] + resi_data["mut"].iloc[-1][1:-1]
    ax.set_title(title)
    plt.subplots_adjust(top=0.92,
                        bottom=0.15,
                        left=0.11,
                        right=0.86,
                        hspace=0.2,
                        wspace=0.2)
    plt.show()
    if return_fig:
        return fig
