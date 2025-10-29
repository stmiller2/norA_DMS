## NGS Data merge and filter pipeline

This NGS data processing pipeline makes use of a HTCondor-based high-throughput computing environment available through the UW-Madison biochemistry department. 
Please see [`usage.md`](usage.md) for more detailed instructions on operating the pipeline and the [dedicated GitHub repo](https://github.com/stmiller2/merge_and_filter_starterpack/) for future updates. In short:

1. Raw gzipped fastq files are copied to the [`Fastq/`](Fastq/) directory
2. Processing parameters are set in the [`params.env`](params.env) file
    - E.g., filepaths, minimum overlap for merging with PEAR, CPUs & RAM to request from HTCondor, quality filtering cutoffs, sample names
    - The parameters set in [`params.env`](params.env) were used for all sequencing analyses in this publication.
3. `./process_ngs.sh params.env` reads the parameters, creates & submits a condor submit file, and creates a static-compiled C++ script for quality filtering.
4. For each sample, the pipeline outputs `good_reads.csv` containing all successfully merged sequencing reads passing the quality filters set in [`params.env`](params.env). `good_reads.csv` is the input for all barcoding and functional scoring analyses.
