#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 09:45:08 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

import numpy as np
import pandas as pd
import pickle
from Bio.Align import substitution_matrices

def load_specificity_data(data_path,
                          remove_low_quality = True,
                          subtract_h2o = True,
                          normalize_by_std = True,
                          fill_dropouts = -3,
                          return_format = 'dictionary',
                          wt_policy = 'WT',
                          water_policy = 'exclude'):
    
    """
    Load and preprocess specificity data.

    Parameters:
        remove_low_quality (bool): Mask scores with high error. Default True.
        subtract_h2o (bool): Subtract water reference scores from f_score. Default True.
        normalize_by_std (bool): Normalize scores and errors by their standard deviation. Default True.
        fill_dropouts (float or int): Replace zero-count dropouts with this value. Default -3.
        return_format (str): 'dictionary' to return per-ligand DataFrames, 'dataframe' to return a combined DataFrame. Default 'dictionary'.
        wt_policy (str): How to handle wild-type row: 'WT' (keep), 'exclude' (drop), 'mut_names' (expand into mutation names). Default 'WT'.
        data_path (str): Path to the pickled specificity data. Default provided.
        water_path (str): Path to the pickled water reference data. Default provided.

    Returns:
        dict or pd.DataFrame: Processed data in the specified format.
    """
    
    specificity_dataframes = pickle.load(open(data_path, 'rb'))
       
    #Copy the dataframe; remove normed columns (will readd them later if req'd)
    processed_dict = {key: value[value.columns.drop(['f_hiq','f_hiq_norm','SE_norm'])].copy() for key, value in specificity_dataframes.items()}
    
    for ligand in processed_dict.keys():
        df = processed_dict[ligand]
        final_col = 'f_score'
        if subtract_h2o:
            water_scores = processed_dict['h2o']
            #Subtract water scores from f_score column
            df['f_score'] -= water_scores['f_score']
        if remove_low_quality:
            #Copy f_score column to f_hiq column...
            df['f_hiq'] = df['f_score']
            #... then, block out scores with error > 1/2 stdev
            mask = df['SE'] > 0.5 * df['f_score'].std()
            df.loc[mask, 'f_hiq'] = np.nan
            final_col = 'f_hiq'
        if normalize_by_std:
            divisor = df[final_col].std()
            #Divide f_hiq and SE by f_hiq's stdev. Store in f_hiq_norm and SE_norm
            df[f'{final_col}_norm'] = df[final_col] / divisor
            df['SE_norm'] = df['SE'] / divisor
            final_col = f'{final_col}_norm'
        if fill_dropouts:
            mask = (df['counts_selected'] == 0) & (df['counts_input'] != 0)
            df.loc[mask, final_col] = fill_dropouts
        if wt_policy == 'exclude':
            df.drop(index='WT', errors='ignore', inplace = True)
        elif wt_policy == 'mut_names':
            wt_row = df.loc['WT']
            df.drop(index='WT', errors='ignore', inplace = True)
            wt_sequence = 'MNKQIFVLYFNIFLIFLGIGLVIPVLPVYLKDLGLTGSDLGLLVAAFALSQMIISPFGGTLADKLGKKLIICIGLILFSVSEFMFAVGHNFSVLMLSRVIGGMSAGMVMPGVTGLIADVSPSHQKAKNFGYMSAIINSGFILGPGIGGFMAEVSHRMPFYFAGALGILAFIMSVVLIHDPKKSTTSGFQKLEPQLLTKINWKVFITPAILTLVLAFGLSAFETLYSLYTSYKVNYSPKDISIAITGGGIFGALFQIYFFDKFMKYFSELTFIAWSLIYSVIVLVLLVIADGYWTIMVISFVVFIGFDMIRPAITNYFSNIAGDRQGFAGGLNSTFTSMGNFIGPLIAGALFDVHIEAPIYMAIGVSLAGVVIVLIEKQHRAKLKEQNM'
            for i, aa in enumerate(wt_sequence, start=1):
                mut_name = f"{aa}{i}{aa}"
                df.loc[mut_name] = wt_row
    
    if water_policy == 'exclude' and 'h2o' in processed_dict:
        del processed_dict['h2o']
    
    if return_format == 'dictionary':
        return processed_dict
    
    if return_format == 'dataframe':
        processed_df = pd.DataFrame()
        for ligand in processed_dict.keys():
            processed_df[ligand] = processed_dict[ligand][final_col]
        return processed_df
    
def physicochemical_distance(mutation, matrix = 'BLOSUM90'):
    """
    Returns the substitution matrix score for a given mutation string (e.g., 'E222A').
    Default is blosum90
    """  
    original_residue = mutation[0]
    mutated_residue = mutation[-1]
    
    mtx = substitution_matrices.load(matrix)
    
    # Look up the score in the matrix (order of keys does not matter)
    pair = (original_residue, mutated_residue)
    if pair not in mtx:
        pair = (mutated_residue, original_residue)
    
    return mtx.get((original_residue, mutated_residue), np.nan)
    
    