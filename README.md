# Energetic and Structural Control of Polyspecificity in a Multidrug Transporter

This repository contains scripts and data for reproducing the analyses published in our PNAS publication *Energetic and structural control of polyspecificity in a multidrug transporter*.  

## Directory Overview

- **[`data/`](data/)** : Raw and processed sequencing data; external data collected for validation experiments and computational predictions.  
- **[`scripts/`](scripts/)** : Code used to process sequencing data and reproduce all analyses and figures in the paper.  
- **[`results/`](results/)** : Output of all analyses used to create the figures in the paper.  

Please see the `README.md` files in each directory for more detailed information.  

## Generalized Workflow

Our high-throughput analyses are based on functional scores, which are calculated from NGS read counts. Sample data are processed according to the following workflow:

1. **Raw sequencing data**  
   Raw paired-end sequencing `.fastq.gz` files come off the sequencer (see [`data/raw/`](data/raw/)).  

2. **Merge and filter reads**  
   Forward and reverse reads are merged using PEAR and filtered for quality (see [`scripts/01_merge_filter_NGS/`](scripts/01_merge_filter_NGS/)).  

3. **Barcode counting and lookup**  
   DNA barcodes are counted and summed to get the total number of observations for each library variant. Barcode-variant associations are determined in [`scripts/02_barcoding/`](scripts/02_barcoding/).  

4. **Functional score calculation**  
   Variant counts from pre- and post-selection samples are compared to compute functional scores. Data pertaining to specificity and efficiency analyses are processed and bundled into Python dictionaries, then exported as `.pk1` files for easy import into later analyses (see [`scripts/03_functional_scoring/`](scripts/03_functional_scoring/)).  

5. **Analysis of functional scores**  
   Data are imported into figure-specific analysis scripts (see [`scripts/04_analyses/`](scripts/04_analyses/)). Depending on the analysis, additional data processing may occur here, such as normalization or filtering out datapoints with high error.  
