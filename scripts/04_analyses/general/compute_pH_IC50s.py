#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 10:11:04 2025

@author: Silas Miller - Raman Lab 
         stmiller2@wisc.edu
"""


base_dir = '/Volumes/sraman4/General/publication_repository/norA_DMS'

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os
os.chdir(f'{base_dir}/scripts/00_utils')
import scoring_utils as sc #for RML estimator to collapse replicates

data = pd.read_csv(f'{base_dir}/data/external/pH_ic50_data.csv')

color_dict = {6: 'tab:green', 7: 'tab:purple'}

def unique_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

def logistic(x, mn, mx, IC50, hill):
    return mn + (mx - mn) / (1 + (x / IC50) ** hill)

def ic50_regression(strain_data):
    fig, ax = plt.subplots()
    for replicate in [1,2,3]:
        rep_mask = strain_data['replicate'] == replicate
        for pH in [6, 7]:
            pH_mask = strain_data['pH'] == pH
            x = strain_data.loc[rep_mask & pH_mask, 'norfloxacin_conc_ng/mL']
            y = strain_data.loc[rep_mask & pH_mask, 'OD600']

            ax.scatter(x, y, color=color_dict[pH])
            ax.set_xscale('log')
            initial_guess = [min(y), max(y), np.median(x), 3.5]
            bounds = (0, [np.inf, np.inf, np.inf, 10])
            popt, pcov = curve_fit(logistic, x, y, p0=initial_guess, bounds=bounds, maxfev=1000)
            mn_fit, mx_fit, IC50_fit, hill_fit = popt
            x_fit = np.logspace(0, np.log10(max(x)), 100)
            y_fit = logistic(x_fit, mn_fit, mx_fit, IC50_fit, hill_fit)
            ax.plot(x_fit, y_fit, label=f'pH {pH}', color=color_dict[pH])

        strain = data['strain'].unique()[0]
        ax.set_title(f'{strain}, norfloxacin IC50')
        ax.set_xlabel('[norfloxacin] (ng/mL)')
        ax.set_ylabel('OD600')
        ax.set_ylim(0, 1.21)
        unique_legend(ax)
    return fig

""" View all IC50 curves """
for strain in ['WT', 'E222A', 'F188C', 'T336E', 'M150E', 
               'P121R', 'I76M', 'H354Q', 'A105E', 'F47D', 
               'S133E', 'F159K', 'H89L', 'N332Q', 'I298D']:
    strain_mask = data['strain'] == strain
    fig = ic50_regression(data.loc[strain_mask])
    plt.show()

#%%
""" Compute IC50s for each replicate separately """
def get_IC50(x, y):
    initial_guess = [min(y), max(y), np.median(x), 3.5]
    bounds = (0, [np.inf, np.inf, np.inf, 10])
    popt, pcov = curve_fit(logistic, x, y, p0=initial_guess, bounds=bounds, maxfev=1000)
    mn_fit, mx_fit, IC50_fit, hill_fit = popt
    err = np.sqrt(np.diag(pcov))[2]
    return IC50_fit, err

ic50s = pd.DataFrame(columns = ['strain', 'rep', 'ic50_6', 'std_6', 'ic50_7', 'std_7'])

i=0
for strain in ['WT', 'E222A', 'F188C', 'T336E', 'M150E', 
               'P121R', 'I76M', 'H354Q', 'A105E', 'F47D', 
               'S133E', 'F159K', 'H89L', 'N332Q', 'I298D']:
    strain_mask = data['strain'] == strain
    strain_data = data.loc[strain_mask]

    for rep in [1, 2, 3]:
        rep_mask = strain_data['replicate'] == rep
        ic50_6, std_6 = np.nan, np.nan
        ic50_7, std_7 = np.nan, np.nan
        for pH in [6, 7]:
            pH_mask = strain_data['pH'] == pH
            x = strain_data.loc[rep_mask & pH_mask, 'norfloxacin_conc_ng/mL']
            y = strain_data.loc[rep_mask & pH_mask, 'OD600']
            ic50, std = get_IC50(x, y)
            if pH == 6:
                ic50_6, std_6 = ic50, std
            elif pH == 7:
                ic50_7, std_7 = ic50, std
        ic50s.loc[i] = {'strain': strain,
                        'rep': rep,
                        'ic50_6': ic50_6,
                        'std_6': std_6,
                        'ic50_7': ic50_7,
                        'std_7': std_7}
        i+=1


#%% 
""" Use RML estimator to find max likelihood IC50 and error """
ic50_rml = pd.DataFrame(columns = ['strain', 'ic50_6', 'ic50_7', 'err_6', 'err_7', 'eps_6', 'eps_7'])

i = 0
for strain in ['WT', 'E222A', 'F188C', 'T336E', 'M150E', 
               'P121R', 'I76M', 'H354Q', 'A105E', 'F47D', 
               'S133E', 'F159K', 'H89L', 'N332Q', 'I298D']:
    strain_mask = ic50s['strain'] == strain
    strain_data = ic50s.loc[strain_mask]
    ic50_6, ic50_7, err_6, err_7, eps_6, eps_7 = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    for pH in [6, 7]:
        betaML, var_betaML, eps = sc.rml_estimator(strain_data[f'ic50_{pH}'], strain_data[f'std_{pH}']**2)
        if pH == 6:
            ic50_6 = betaML
            err_6 = np.sqrt(var_betaML)
            eps_6 = eps
        elif pH == 7:
            ic50_7 = betaML
            err_7 = np.sqrt(var_betaML)
            eps_7 = eps
    
    ic50_rml.loc[i] = {'strain': strain,
                       'ic50_6': ic50_6,
                       'ic50_7': ic50_7,
                       'err_6': err_6,
                       'err_7': err_7,
                       'eps_6': eps_6,
                       'eps_7': eps_7}
    i+=1

ic50_rml.to_csv(f'{base_dir}/results/general/pH_ic50s.csv',
                index = None)

#%%     
""" Get clonal dFpH scores """
clonal_fscores = ic50_rml.copy()[['strain']]
ref = 'WT'
ref_mask = ic50_rml['strain'] == ref

for pH in [6, 7]:
    ref_ic50 = ic50_rml.loc[ref_mask, f'ic50_{pH}'].iloc[0]
    ref_err = ic50_rml.loc[ref_mask, f'err_{pH}'].iloc[0]
    
    ic50 = ic50_rml[f'ic50_{pH}']
    err = ic50_rml[f'err_{pH}']
    
    norm_ic50 = ic50 / ref_ic50
    norm_err = norm_ic50 * np.sqrt((err / ic50)**2 + (ref_err / ref_ic50)**2)
    
    clonal_fscores[f'f_{pH}'] = np.log2(norm_ic50)
    clonal_fscores[f'err_{pH}'] = norm_err / (norm_ic50 * np.log(2))
    
clonal_fscores['dFpH'] = clonal_fscores['f_7'] - clonal_fscores['f_6']
clonal_fscores['err_dFpH'] = np.sqrt(clonal_fscores['err_6']**2 + clonal_fscores['err_7']**2)
clonal_fscores = clonal_fscores[['strain', 'f_6', 'f_7', 'err_6', 'err_7', 'dFpH', 'err_dFpH']] 

clonal_fscores.to_csv(f'{base_dir}/results/general/clonal_dFpH.csv',
                      index = None)


