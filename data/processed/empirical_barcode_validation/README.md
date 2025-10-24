## Directory Overview

In this analysis, we ran test selections using a subset of the NorA library that spanned only ~240bp, and thus could be analyzed by directly sequencing the variable region of the NorA gene. By computing functional scores from counts derived from both direct variant sequencing and barcoding, we determined that our methods of barcoding, barcode mapping, quality filtering, and functional score calculation report the same functional information as the gold-standard of direclty sequencing library variants.

- **[`good_reads/`](good_reads/)** : Sequencing reads after merging and filtering for quality. These files are too large to host on GitHub, but are available on request.
- **[`barcode_counts/`](barcode_counts/)** : CSVs containing raw read counts per mapped barcode (only for barcode-sequenced samples).
- **[`variant_counts/`](variant_counts/)** : CSVs containing raw read counts per variant (Either directly counted [for variant-sequenced samples] or summed from `barcode_counts` [for barcode-sequenced samples]).
- **[`replicate_functional_scores/`](replicate_functional_scores/)** : Functional scores computed for individual replicates.
- **[`scoring_figs/`](scoring_figs/)** : PNG images illustrating scoring QC like distribution of pre-selection libraries, correlation of variant frequencies between replicates, and correlation of functional scores between replicates.
- **[`final_functional_scores/`](final_functional_scores/)** : Functional scores for each selection after merging replciates.

Analysis scripts comparing barcode- and variant-derived functional scores can be found in Figure S1 of the paper (see relevant script [here](../../../scripts/04_analyses/supplementary/figS01_barcoding_plots.py) and results [here](../../../results/supplementary/figS01_barcoding))
