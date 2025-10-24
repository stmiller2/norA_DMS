#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 09:45:26 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""

import numpy as np
import pandas as pd
import pickle

def load_efficiency_data(data_path,
                         remove_low_quality = True,
                         subtract_h2o = True,
                         normalize_by_std = True,
                         return_format = 'dFpH_dataframe',
                         wt_policy = 'WT'):

    efficiency_dataframes = pickle.load(open(data_path, 'rb'))
       
    #Copy the dataframe; remove normed columns (will readd them later if req'd)
    processed_dict = {key: efficiency_dataframes[key][['f_score','SE']].copy() for key in ['h2o_pH60', 'h2o_pH70', 'acr_pH60', 'acr_pH70', 'nor_pH60', 'nor_pH70',]}
    for key, df in processed_dict.items():
        score_col = 'f_score'
        err_col = 'SE'
        if subtract_h2o:
            if 'h2o' not in key:
                water_scores = processed_dict[f'h2o_{key[-4:]}']
                df['f_score'] -= water_scores['f_score']
        if remove_low_quality:
            df['f_hiq'] = df['f_score']
            mask = df['SE'] > 0.5 * df['f_score'].std()
            df.loc[mask, 'f_hiq'] = np.nan
            score_col = 'f_hiq'
        if normalize_by_std:
            divisor = df[score_col].std()
            df[f'{score_col}_norm'] = df[score_col] / divisor
            df[f'{err_col}_norm'] = df[err_col] / divisor
            score_col = f'{score_col}_norm'
            err_col = f'{err_col}_norm'
        if wt_policy == 'exclude':
            df.drop(index='WT', errors='ignore', inplace = True)
        elif wt_policy == 'mut_names':
            wt_row = df.loc['WT']
            df.drop(index='WT', errors='ignore', inplace = True)
            wt_sequence = 'MNKQIFVLYFNIFLIFLGIGLVIPVLPVYLKDLGLTGSDLGLLVAAFALSQMIISPFGGTLADKLGKKLIICIGLILFSVSEFMFAVGHNFSVLMLSRVIGGMSAGMVMPGVTGLIADVSPSHQKAKNFGYMSAIINSGFILGPGIGGFMAEVSHRMPFYFAGALGILAFIMSVVLIHDPKKSTTSGFQKLEPQLLTKINWKVFITPAILTLVLAFGLSAFETLYSLYTSYKVNYSPKDISIAITGGGIFGALFQIYFFDKFMKYFSELTFIAWSLIYSVIVLVLLVIADGYWTIMVISFVVFIGFDMIRPAITNYFSNIAGDRQGFAGGLNSTFTSMGNFIGPLIAGALFDVHIEAPIYMAIGVSLAGVVIVLIEKQHRAKLKEQNM'
            for i, aa in enumerate(wt_sequence, start=1):
                mut_name = f"{aa}{i}{aa}"
                df.loc[mut_name] = wt_row
    dFpH_dict = {}
    for lig in ['acr', 'nor']:
        dFpH_dict[lig] = pd.DataFrame()
        dFpH_dict[lig]['f_pH60'] = processed_dict[f'{lig}_pH60'][score_col]
        dFpH_dict[lig]['SE_pH60'] = processed_dict[f'{lig}_pH60'][err_col]
        dFpH_dict[lig]['f_pH70'] = processed_dict[f'{lig}_pH70'][score_col]
        dFpH_dict[lig]['SE_pH70'] = processed_dict[f'{lig}_pH70'][err_col]
        dFpH_dict[lig]['dFpH'] = dFpH_dict[lig]['f_pH70'] - dFpH_dict[lig]['f_pH60']
        dFpH_dict[lig]['dFpH_SE'] = np.sqrt(dFpH_dict[lig]['SE_pH70']**2 + dFpH_dict[lig]['SE_pH60']**2)
        if remove_low_quality:
            mask = dFpH_dict[lig]['dFpH_SE'] > 0.5 * dFpH_dict[lig]['dFpH'].std()
            dFpH_dict[lig].loc[mask, 'dFpH'] = np.nan
    if return_format == 'scores_dictionary':
        return processed_dict
    if return_format == 'scores_dataframe':
        scores_df = pd.DataFrame()
        for lig in ['acr', 'nor']:
            scores_df[f'{lig}_pH60'] = dFpH_dict[lig]['f_pH60']
            scores_df[f'{lig}_pH60_SE'] = dFpH_dict[lig]['SE_pH60']
            scores_df[f'{lig}_pH70'] = dFpH_dict[lig]['f_pH70']
            scores_df[f'{lig}_pH70_SE'] = dFpH_dict[lig]['SE_pH70']
            scores_df[f'{lig}_dFpH'] = dFpH_dict[lig]['dFpH']
            scores_df[f'{lig}_SE'] = dFpH_dict[lig]['dFpH_SE']
        return scores_df
    if return_format == 'dFpH_dictionary':
        return dFpH_dict
    if return_format == 'dFpH_dataframe':
        dFpH_df = pd.DataFrame()
        for lig in ['acr', 'nor']:
            dFpH_df[f'{lig}_dFpH'] = dFpH_dict[lig]['dFpH']
            dFpH_df[f'{lig}_SE'] = dFpH_dict[lig]['dFpH_SE']
        return dFpH_df
