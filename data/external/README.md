## Directory Overview

This directory contains non-sequencing data used for clonal validations and other analyses in the publication. The identity and use of each file is provided below.

#### NorA structures:
Used for [calculating residue distances](../../scripts/04_analyses/general/compute_residue_distances.py) (7lo8, 9b3m, alphafold) and for computing [Rosetta](../../scripts/04_analyses/rosetta_ddG_pipeline/)- and ThermoMPNN-predicted ΔΔG values (7lo8, 9b3m, and 9b3l).
- **`7lo8.pdb`** : NorA structure in the outward-open conformation (see [Brawley 2022](https://doi.org/10.1038/s41589-022-00994-9)).
- **`9b3m.pdb`** : NorA structure in the inward-open conformation (see [Xie 2025](https://doi.org/10.1038/s41467-024-54986-5)).
- **`9b3l.pdb`** : NorA structure in the occluded conformation (see [Xie 2025](https://doi.org/10.1038/s41467-024-54986-5)).
- **`norA_AF.pdb`** : NorA structural model predicted by alphafold (inward-open conformation).
#### ΔΔG scanning data:
Rosetta modeling output for full DMS [ΔΔG scans](../../scripts/04_analyses/rosetta_ddG_pipeline/), used for predicting variant stability to identify potential functional hotspots. Scanning was performed in all conformations to differentiate stability mutants (unstable in all conformations) from conformational dynamics mutants (stable, but unable to access some conformations).
- **`7lo8_ddG_results.csv`** : Outward-open ΔΔG scan summary 
- **`9b3m_ddG_results.csv`** : Inward-open ΔΔG scan summary
- **`9b3l_ddG_results.csv`** : Occluded ΔΔG scan summary
> Full scanning data (e.g., Rosetta files for individual mutations) are too large to host here and are available on request.
- **`computational_stability.csv`** : Combined csv containing ΔΔG values for each conformation, predicted by both Rosetta and ThermoMPNN. See relevant analysis script [here](../../scripts/04_analyses/supplementary/figS13_hotspots_plots.py); Results are shown in [Supplementary Fig. S13](../../results/supplementary/figS13_hotspots/) and [Supplementary Table S2](../../results/supplementary/tableS02_hotspots/).
#### Other data:
- **`consurf_scores.xlsx`** : Output from the ConSurf web server showing evolutionary conservation for each position of NorA. Conservation is analyzed [here](../../scripts/04_analyses/supplementary/figS13_hotspots_plots.py) and shown in [Supplementary Fig. 13A](../../results/supplementary/figS13_hotspots/)
- **`ethidium_efflux_data.xlsx`** : Platereader data showing ethidium fluorescence over time, analyzed [here](../../scripts/04_analyses/fig5_plots.py) to biophysically validate growth-based functional scores in [Fig. 5A and 5B](../../results/fig5/)
- **`hibit_data.csv`** : Platereader data showing HiBiT luminescence of selected HiBiT-tagged variants, analyzed [here](../../scripts/04_analyses/fig5_plots.py) to measure membrane-localized protein abundance of library variants, shown in [Fig. 5E](../../results/fig5/) and Supplementary Figures [13E](../../results/supplementary/figS13_hotspots/), [16B](../../results/supplementary/figS16_pH_ic50_validations/), and [18](../../results/supplementary/figS18_sensitivity_analysis/).
- **`substrates_ic50_data.csv`** : Platereader data showing OD600 values after 16 hours of growth in varying concentrations of each tested drug, for both WT NorA and the catalytically inactive mutant E222A, analyzed [here](../../scripts/04_analyses/supplementary/figS03_substrate_ic50s_plots.py) and shown in [Supplementary Fig. 3](../../results/supplementary/figS03_substrate_ic50s/).
- **`pH_ic50_data.csv`** : Platereader data showing OD600 values after 16 hours of growth in varying concentrations of norfloxacin at pH 6.0 or pH 7.0 for several library variants, analyzed [here](../../scripts/04_analyses/fig5_plots.py) and [here](../../scripts/04_analyses/supplementary/figS16_pH_IC50_validations_plots.py), and shown in [Fig. 5C and 5D](../../results/fig5/) and [Supplementary Fig. 16](../../results/supplementary/figS16_pH_ic50_validations/).
