# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 09:15:04 2023
@author: Silas
Version 1.0

For mapping barcodes to variants for NGS reads

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import ceil
from time import time
from IPython import get_ipython, display

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

def extract_barcodes(seq_df, prebc, postbc):
    """
    Extract barcode sequences from NGS reads in a DataFrame.
    Adds 'bc' column with the barcode between prebc and postbc,
    and 'bc_len' column with the barcode length. Missing barcodes are NaN.
    """
    seq_df_out = seq_df.copy()
    
    ### Extract barcodes
    stime = time()
    print('Extracting barcodes... ', end='')
    seq_df_out['bc'] = seq_df_out['read'].apply(
        lambda x: x[x.find(prebc)+len(prebc):x.find(postbc)]
        if (prebc in x and postbc in x)
        else float('nan')
        )
    print(f'Done, {round(time()-stime, 2)} seconds\n')

    ### Get barcode lengths
    stime = time()
    print('Calculating barcode lengths... ', end='')
    seq_df_out['bc_len'] = seq_df_out['bc'].str.len()
    print(f'Done, {round(time()-stime, 2)} seconds\n')

    return seq_df_out

def translate_variants(seq_df, prerf, postrf):
    """
    Extract translated variants from NGS reads in a DataFrame.
    Returns a copy with a 'translation' column containing the AA sequence
    between prerf and postrf. Untranslated/unbarcoded reads are removed,
    and the 'read' column is dropped to save space.
    """
    seq_df_out = seq_df.copy()
    
    stime = time()
    print('Translating variants... ', end='')
    seq_df_out['translation'] = seq_df_out['read'].apply(
        lambda x: translate(x[x.find(prerf)+len(prerf):x.find(postrf)])
        if (prerf in x and postrf in x)
        else float('nan')
        )
    seq_df_out = seq_df_out.dropna().drop('read', axis = 1)
    print(f'Done, {round(time()-stime, 2)} seconds')

    return seq_df_out


def group_barcodes(seq_df):
    """
    Groups barcodes by Levenshtein distance = 1. Only barcodes seen >= 2 times 
    are considered. Most abundant barcodes are parents; Levenshtein-1 barcodes 
    become children.
    """
    stime = time()
    print('Grouping barcodes according to Levenshtein distance = 1... ', end='')

    #Filter barcodes seen at least twice and check for NaNs
    counts = seq_df['bc'].value_counts()
    counts = counts[counts >= 2]
    assert counts.index.isna().sum() == 0, 'NaN values in bc column - replace with ""!'

    #Function to generate all possible Levenshtein distance = 1 sequences
    def lev1_bcs(bc):
        output = set()
        for i in range(len(bc)):
            output.add(bc[:i] + bc[i+1:]) #One-base deletion
            for base in 'ATGC': #One-base substitution
                output.add(bc[:i] + base + bc[i+1:])
            for base in 'ATGC': #One-base insertion
                output.add(bc[:i] + base + bc[i:])
        output.discard(bc) #Remove original sequence
        return output

    groups = {}
    searchable = {bc for bc in counts.index}
    found = set()
    
    for bc in counts.index:
        strings = lev1_bcs(bc)
        for string in strings:
            if string in searchable and string not in found:
                if bc in groups:
                    groups[bc].add(string)
                    found.add(string)
                else:
                    groups[bc] = {string}
                    found.add(bc)
                    found.add(string)

    print(f'Done, {round(time()-stime,2)} seconds', end = '\n\n')
    return groups


def collapse_groups(seq_df, groups):
    """
    Collapse child barcodes onto their parent barcode within Levenshtein-1 groups.
    Parent barcodes remain unchanged; children are overwritten with the parent sequence.
    """
    stime = time()
    print('Collapsing groups onto parent barcode... ', end='')

    seq_df_out = seq_df.copy()

    #Create mapping from child -> parent
    reverse_dict = {child: parent for parent, children in groups.items() for child in children}

    #Replace barcodes with parent if in mapping
    seq_df_out['bc'] = seq_df_out['bc'].map(lambda x: reverse_dict.get(x, x))

    print(f'Done, {round(time()-stime,2)} seconds\n')

    #Update barcode lengths
    stime = time()
    print('Resetting barcode lengths... ', end='')
    seq_df_out['bc_len'] = seq_df_out['bc'].str.len()
    print(f'Done, {round(time()-stime,2)} seconds\n')

    return seq_df_out

def sample_cutoffs_and_filter(seq_df, cutoff_sample=25, filepath='../figures'):
    """
    Sample potential pair read count cutoffs and plot surviving unique pairs.
    User selects a cutoff; dataframe is filtered to include only pairs at or 
    above that cutoff. Returns filtered dataframe with added 'bc_freq' column.
    """
    stime = time()
    print(f'Sampling pair cutoffs from 0-{cutoff_sample}... ', end = '')
    
    seq_df_out = seq_df.copy()
    seq_df_out['pair'] = seq_df_out['bc'] + '-' + seq_df_out['translation']
    counts = seq_df_out['pair'].value_counts()
    
    #Surviving pairs for each cutoff
    surviving = [counts.gt(cutoff).sum() for cutoff in range(cutoff_sample)]

    #Plot
    get_ipython().run_line_magic('matplotlib', 'inline')
    fig, ax = plt.subplots()
    ax.plot(range(cutoff_sample), surviving)
    ax.set_ylim(0, surviving[-1]*10)
    ax.set_xlabel('Pair cutoff')
    ax.set_ylabel('Unique barcode-variant pairs')
    ax.set_title('Surviving barcode-variant pairs vs. pair cutoff')
    display.display(fig)
    
    print(f'Done, {round(time()-stime,2)} seconds\n\n')
    print('Please see the cutoff sampling graph in the Plots pane to determine an appropriate pair_cutoff.')
    
    pair_cutoff = int(input('Enter desired integer pair_cutoff: '))
    
    #Annotate cutoff and save plot
    ax.plot(pair_cutoff-1, surviving[pair_cutoff-1], 'o', color = 'tab:blue')
    ax.text(pair_cutoff-1, surviving[pair_cutoff-1]+(.4*surviving[len(surviving)-1]), 
             f'Chosen cutoff: {pair_cutoff} \n'
             f'Remaining unique pairs: {surviving[pair_cutoff-1]:,}')
    fig.savefig(f'{filepath}/pair_cutoffs.png', dpi = 200, bbox_inches = 'tight')
    print(f'Saved pair cutoffs plot to {filepath}/pair_cutoffs.png. ')
    
    #Filter dataframe
    print(f'Filtering for pairs observed at or above pair_cutoff ({pair_cutoff})... ', end = ''); stime = time()
    pairs_above_cutoff = counts[counts >= pair_cutoff].index
    seq_df_out = seq_df_out[seq_df_out['pair'].isin(pairs_above_cutoff)].drop(columns='pair')

    # Add barcode frequency column
    bc_counts = seq_df_out['bc'].value_counts().rename('bc_freq')
    seq_df_out = seq_df_out.merge(bc_counts, left_on='bc', right_index=True)
    
    print(f'Done, {round(time()-stime,2)} seconds\n', end = '\n\n')
    return seq_df_out

def map_barcodes(seq_df):
    """
    Map barcodes to their variants.
    Returns a dataframe with one row per unique barcode-variant pair,
    including 'map_freq' (pair count), 'bc_len', and 'bc_freq'.
    """
    stime = time()
    print('Mapping barcodes... ', end = '')
    
    mapped_barcodes = seq_df.groupby(['bc','translation']).size().reset_index(name = 'map_freq')
    mapped_barcodes = mapped_barcodes.merge(seq_df[['bc', 'bc_len', 'bc_freq']].drop_duplicates(), 
                                            on = 'bc', how = 'left')
    mapped_barcodes = mapped_barcodes[['bc', 'bc_len', 'translation', 'bc_freq', 'map_freq']]
    
    print(f'Done, {round(time()-stime,2)} seconds', end = '\n\n')
    return mapped_barcodes

def sample_chastity_and_filter(seq_df, filepath='../figures'):
    """
    For each barcode, calculate chastity value (counts for most common mapping
    divided by counts for second-most common mapping), plot the distribution, 
    get a user-defined chastity cutoff, and filter the dataframe.
    """ 
    stime = time()
    print('Calculating chastity values and plotting distribution... ', end = '')
    
    #Compute chastity values per barcode
    chastity_values = pd.DataFrame()
    chastity_values['bc'] = seq_df['bc'].unique()
    chastity_values['chastity_value'] = seq_df.groupby('bc').apply(
        lambda mmdata: mmdata['map_freq'].nlargest(2).iloc[0] /
        mmdata['map_freq'].nlargest(2).sum()
        ).tolist()
    chastity_values['fraction_of_total'] = seq_df.groupby('bc').apply(
        lambda mmdata: mmdata['map_freq'].max() / 
        mmdata['map_freq'].sum()
        ).tolist()
    
    #Barcodes that do not map perfectly
    imperfect_chastity_values = chastity_values[
        chastity_values['chastity_value'] != 1.0
        ]['chastity_value']
    #Merge with main dataframe
    seq_df_out = pd.merge(seq_df, chastity_values, on='bc', how='left')
    
    #Plot histogram of chastity values
    fig, ax = plt.subplots()
    y, x, _ = ax.hist(imperfect_chastity_values, bins = ceil(len(imperfect_chastity_values)/50))
    perfectmaps_uniqct = len(chastity_values)-len(imperfect_chastity_values)
    ax.text(x[0]+(x[len(x)-1]-x[0])/15, y.max()*.9, 
             f'{perfectmaps_uniqct:,} barcodes map perfectly',
             bbox = dict(facecolor='none', edgecolor='black', boxstyle='round,pad=.25'))
    ax.set_xlabel('Chastity value')
    ax.set_ylabel('Count (number of barcodes)')
    ax.set_title('Distribution of chastity values among multimapping barcodes')
    display.display(fig)
    
    print(f'Done, {round(time()-stime,2)} seconds\n\n')
    print('Please see the chastity value distribution in the plots pane to determine an appropriate chastity_cutoff.')
    
    #Get user-defined cutoff
    chastity_cutoff = float(input('Enter desired float chastity_cutoff: '))
    
    #Annotate chosen cutoff and save figure
    ax.vlines(chastity_cutoff, 0, y.max(), linestyle='dashed', color='black')
    ax.text(chastity_cutoff - 0.2, y.max()/2, f'Chosen cutoff: {chastity_cutoff}')
    fig.savefig(f'{filepath}/chastity_value_dist.png', dpi=200, bbox_inches='tight')
    print(f'Saved chastity value distribution to {filepath}/chastity_value_dist.png')
    
    #Filter by chastitity cutoff
    print(f'Filtering chastity value >= {chastity_cutoff}... ', end = ''); stime = time()
    seq_df_out = seq_df_out[seq_df_out['chastity_value'] >= chastity_cutoff].reset_index(drop = True)
    max_map_freqs = seq_df_out.groupby('bc')['map_freq'].max()
    seq_df_out = seq_df_out.loc[seq_df_out.set_index(['bc', 'map_freq']).index.isin(
            max_map_freqs.reset_index().set_index(['bc', 'map_freq']).index
            )].reset_index(drop = True)
    print(f'Done, {round(time()-stime,2)} seconds', end = '\n\n')
    
    return seq_df_out

def id_mutants(seq_df, ref_lib):
    """
    Identifies mutants and names them with standard mutant notation in a new
    column called 'mut'. Uses ref_lib csv input.

    """
    stime = time()
    print('Identifying expected mutants... ', end = '')
    seq_df_out = seq_df.merge(ref_lib, left_on = 'translation', 
                              right_on = 1, how = 'left')
    seq_df_out = seq_df_out.rename(columns = {0: 'mut'}).drop(1, axis = 1)
    print(f'Done, {round(time()-stime,2)} seconds')
    
    return seq_df_out
