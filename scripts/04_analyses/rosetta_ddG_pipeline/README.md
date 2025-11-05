## RosettaMP ΔΔG scanning pipeline

This pipeline is based upon Raman lab scripts for soluble ΔΔG scans and [Tiemann 2023](https://doi.org/10.1016/j.bpj.2022.12.031) and is written for use on the HTCondor-based high-throughput computing environment available through the UW-Madison biochemistry department.
Please see [`usage.txt`](usage.txt) for more detailed instructions on operating the pipeline, or visit the dedicated [GitHub repository](https://github.com/stmiller2/membrane_ddG_scan) for future updates. In short:

1. The protein of interest is oriented in the membrane with the [PPM web server](https://opm.phar.umich.edu/ppm_server3_cgopm).
2. [`./prep_inputs.sh`](prep_inputs.sh) cleans the PDB ([`scripts/clean_pdb3.py`](scripts/clean_pdb3.py)), generates a membrane spanfile (`scripts/spanfile_from_pdb.linuxgccrelease`), and energy-minimizes the structure with a cartesian fastrelax protocol (rosettascripts [`scripts/mp_cart_relax.xml`](scripts/mp_cart_relax.xml) using [`scripts/f19_cart_1.5.wts`](scripts/f19_cart_1.5.wts) weights)
   - PDB cleaning and spanfile generation are run locally. Relaxation is run on the compute cluster using the `rosettacommons/rosetta:latest` Docker image maintained by Rosetta Commons.
3. [`./mp_cartddG_pipeline.sh`](mp_cartddG_pipeline.sh) creates mutfiles, bash scripts to run the `cartesian_ddG` Rosetta application, and HTCondor submit files for each mutation.
4. When all jobs are complete, [`./parse_results.sh`](parse_results.sh) takes the average WT and MUT energies for each iteration and subtracts them to get the predicted ΔΔG.
5. Parsed ΔΔG predictions are saved to an `output.csv` file. After running this pipeline using NorA structures in each major conformation (7LO8, 9B3L, 9B3M), we combined these outputs into a single file along with ThermoMPNN predictions found at [`data/external/computational_stability.csv`](../../data/external/computational_stability.csv)
