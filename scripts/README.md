## Directory Overview

- **[`00_utils/`](00_utils/)** : Helper functions called by other scripts for barcode mapping, functional scoring, clustering, importing data, and creating plots for the manuscript.
- **[`01_merge_filter_NGS/`](01_merge_filter_NGS/)** : NGS sequencing data processing pipeline used for merging paired-end reads and filtering for quality in all sequencing-based analyses.
- **[`02_barcoding/`](02_barcoding/)** : Scripts used for mapping barcodes to variants and generating functional scores for validating variant-barcode associations.
- **[`03_functional_scoring/`](03_functional_scoring/)** : Scripts for counting barcodes and computing functional scores, and performing data post-processing steps that are common to all analyses.
- **[`04_analyses/`](04_analyses/)** : Figure-by-figure analysis scripts for recreating each plot in the manuscript main text and supplementary figures. Our pipeline for computing predicted ΔΔG values using the Rosetta membrane framework is also provided here.
