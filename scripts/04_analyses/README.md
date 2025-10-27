## Directory Overview

These scripts can be used to reproduce all analyses included in the paper. Each script accesses data from [`data/processed/`](../data/processed) and produces output in the [`results/`](results/) directory.

#### General analysis scripts
- **[`general/`](general/)** : Analyses used in multiple figures (regressions of pH-perturbed IC50 assays, computing distances of each residue from the binding site or coupling residues). 
- **[`rosetta_ddG_pipeline/`](rosetta_ddG_pipeline/)** : Scripts for running computational ΔΔG predictions using the Rosetta membrane framework.
- **[`agglomerative_clustering.py`](agglomerative_clustering.py)** : Script for selecting clustering parameters (cophenetic correlation and gap statistic analysis) and performing hierarchical clustering to produce [`cluster_assignments.csv`](../results/clustering/cluster_assignments.csv), which is used in further analyses. See the results directory [here](../results/clustering/).
#### Figure-by-figure analysis scripts
- **[`fig1_plots.py`](fig1_plots.py)** : Script for reproducing all analyses shown in Figure 1. See results directory [here](../results/fig1).
- **[`fig2_plots.py`](fig2_plots.py)** : Script for reproducing all analyses shown in Figure 2. See results directory [here](../results/fig2).
- **[`fig3_plots.py`](fig3_plots.py)** : Script for reproducing all analyses shown in Figure 3. See results directory [here](../results/fig3).
- **[`fig4_plots.py`](fig4_plots.py)** : Script for reproducing all analyses shown in Figure 4. See results directory [here](../results/fig4).
- **[`fig5_plots.py`](fig5_plots.py)** : Script for reproducing all analyses shown in Figure 5. See results directory [here](../results/fig5).
- **[`fig6_plots.py`](fig6_plots.py)** : Script for reproducing all analyses shown in Figure 6. See results directory [here](../results/fig6).
- **[`supplementary/`](supplementary/)** : Scripts for reproducing all analyses shown in supplementary figures and tables. See results directories [here](../results/supplementary).
