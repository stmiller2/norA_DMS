# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:39:18 2023
@author: Silas
Version 2.0
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from time import time

codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate(seq):
    """
    Translate a DNA sequence to protein. Stop codons are denoted with '_'.
    Sequences containing uncalled bases ('N') return NaN.
    """
    if 'N' in seq:
        return np.nan
    return ''.join(codontable[seq[i:i+3]] for i in range(0, len(seq)-len(seq)%3, 3))

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

def extract_subsequence(seq_col, preseq, postseq):
    """
    Extract a subsequence from a dataframe column (series) that occurs between
    two other known subsequences. Returns a pandas series containing the 
    subsequence, or NaN if either constant region is missing.
    """    
    start_time('Extracting subsequences...')
    
    subseqs = seq_col.apply(
        lambda x: x[x.find(preseq)+len(preseq):x.find(postseq)]
        if (preseq in x and postseq in x)
        else float('nan')
        )
    
    end_time()
    return subseqs

def count_barcodes(bc_col, lookup_table):
    """
    Counts barcode occurences and merges counts with a barcode-to-variant
    lookup table. Returns a copy of the lookup table with 'counts' column
    added.
    """
    start_time('Counting barcodes and merging with lookup table...')
    
    lookup_table['bc'] = lookup_table['bc'].fillna('')
    counts = pd.merge(bc_col.value_counts().reset_index(),
                      lookup_table,
                      how = 'right',
                      left_on = 'bc',
                      right_on = 'bc')

    counts = counts.rename(columns={'count': 'counts',})
    cols = ['bc','mut','counts']+[col for col in counts.columns if col not in ['bc','mut','counts']]
    counts = counts[cols]
    counts['counts'] = counts['counts'].replace(np.nan, 0)
    
    end_time()
    return counts

def sum_barcode_counts(df):
    """
    Sums barcode counts on a per-variant basis. Returns a new df with unique
    muts and the total number of counts in the df input.
    """
    start_time('Summing barcode counts...')
    
    sums = df.groupby('mut')['counts'].sum().reset_index()
    cols = [col for col in df.columns if col not in ['bc', 'counts']]
    sums = sums.merge(df[cols].drop_duplicates(), on='mut', how = 'left')
    
    end_time()
    return sums

def correlation_scatter(values_1, values_2, title = '', xlab = 'Rep. 1', ylab = 'Rep. 2', 
                        corr = 'Spearman', size = (5,5), xlim = None, ylim = None,
                        xscale = None, yscale = None, show = True):
    """
    Compares two datasets on a scatterplot and calculates R-squared and p value
    Returns a scatterplot pyplot fig
    """
    # Clean invalid values
    data = pd.DataFrame({'x': values_1, 'y': values_2}).replace([np.inf, -np.inf], np.nan).dropna()
    if data.empty:
        raise ValueError("No valid data points after removing NaN/Inf.")
    
    fig, ax = plt.subplots(figsize=size if size else None)
    ax.scatter(data['x'], data['y'], s=20, alpha=0.5)
    
    if corr in ('Pearson', 'Spearman'):
        func = pearsonr if corr == 'Pearson' else spearmanr
        corr_coef, _ = func(data['x'], data['y'])
        ax.annotate(f'{corr} R = {corr_coef:.2f}', (0.1, 0.8),
                    xycoords='axes fraction',
                    bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8, ec='black'))
    elif corr is not None:
        raise ValueError("corr must be 'Pearson', 'Spearman', or None")
    
    ax.set(xlabel=xlab, ylabel=ylab, title=title)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if xscale: ax.set_xscale(xscale)
    if yscale: ax.set_yscale(yscale)
    
    if show: plt.show()
    return fig

def variant_fn_scores(input_counts, selected_counts, ref = 'WT', pc = 0.5, log = 'log2'):
    """
    Merges input and selected variant counts into a dataframe with computed
    functional scores and standard errors. Output is a new df with columns 
    'bc', 'mut', 'counts_input', 'counts_selected', 'f_score', and 'SE'. 
    'f_score' is the reference-normalized log ratio of post-selection counts to 
    pre-selection counts. 'SE' is calculated for each f-score by Poisson 
    assumptions (see Enrich2 paper)
    """
    start_time('Getting functional scores...')
    
    merged_counts = input_counts.merge(selected_counts, on='mut', 
                                       suffixes=('_input', '_selected'))
    avg_ref_counts_inp = merged_counts[merged_counts['mut'] == ref]['counts_input'].mean()
    avg_ref_counts_sel = merged_counts[merged_counts['mut'] == ref]['counts_selected'].mean()
    if log == 'log2':
        merged_counts['f_score'] = np.log2((merged_counts['counts_selected']+pc)/(avg_ref_counts_sel+pc)) - np.log2((merged_counts['counts_input']+pc)/(avg_ref_counts_inp+pc))
    if log == 'ln':
        merged_counts['f_score'] = np.log((merged_counts['counts_selected']+pc)/(avg_ref_counts_sel+pc)) - np.log((merged_counts['counts_input']+pc)/(avg_ref_counts_inp+pc))
    merged_counts['SE'] = np.sqrt((1/(merged_counts['counts_input']+pc))+(1/(avg_ref_counts_inp+pc))+(1/(merged_counts['counts_selected']+pc))+(1/(avg_ref_counts_sel+pc)))
    
    end_time()
    return merged_counts


def rml_estimator(y, sigma2i, iterations=100):
    """
    This function from Erich2 source code
    
    Implementation of the robust maximum likelihood estimator.
        ::
            @book{demidenko2013mixed,
              title={Mixed models: theory and applications with R},
              author={Demidenko, Eugene},
              year={2013},
              publisher={John Wiley \& Sons}
            }
    """
    w = 1 / sigma2i
    sw = np.sum(w, axis=0)
    beta0 = np.sum(y * w, axis=0) / sw

    sigma2ML = np.sum((y - np.mean(y, axis=0)) ** 2 / (beta0 - 1), axis=0)
    eps = np.zeros(beta0.shape)
    betaML = None
    for _ in range(iterations):
        w = 1 / (sigma2i + sigma2ML)
        sw = np.sum(w, axis=0)
        sw2 = np.sum(w ** 2, axis=0)
        betaML = np.sum(y * w, axis=0) / sw
        sigma2ML_new = (
            sigma2ML
            * np.sum(((y - betaML) ** 2) * (w ** 2), axis=0)
            / (sw - (sw2 / sw))
        )
        eps = np.abs(sigma2ML - sigma2ML_new)
        sigma2ML = sigma2ML_new
    var_betaML = 1 / np.sum(1 / (sigma2i + sigma2ML), axis=0)
    return betaML, var_betaML, eps

def collapse_replicates(*reps):
    """Collapse multiple replicate DataFrames using the RML estimator."""
    start_time('Collapsing replicates...')
    merged = pd.concat(reps, ignore_index=True)

    def process_group(group):
        y, sigma2i = group['f_score'].values, group['SE'].values**2
        betaML, var_betaML, eps = rml_estimator(y, sigma2i)
        return pd.Series({
            'f_score': betaML,
            'counts_input': group['counts_input'].sum(),
            'counts_selected': group['counts_selected'].sum(),
            'SE': np.sqrt(var_betaML),
            'eps': eps,
            'mut': group['mut'].iloc[0]
        })

    master = merged.groupby('mut').apply(process_group).reset_index(drop=True)
    
    end_time()
    return master

