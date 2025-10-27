## NGS Data merge and filter pipeline

This NGS data processing pipeline makes use of a HTCondor-based high-throughput computing environment available through the UW-Madison biochemistry department. 
Please see [`usage.txt`](usage.txt) for more detailed instructions on operating the pipeline. In short:

1. Raw gzipped fastq files are copied to the [`Fastq/`](Fastq/) directory
2. Processing parameters are set in the [`params.csv`](params.csv) file
    - E.g., filepaths, minimum overlap for merging with PEAR, CPUs & RAM to request from HTCondor, quality filtering cutoffs, sample names
    - The parameters set in [`params.csv`](params.csv) were used for all sequencing analyses in this publication.
3. `./process_ngs.sh params.csv` reads the parameters, creates & submits a condor submit file, and creates a static-compiled C++ script for quality filtering.
    - The quality filtering script is created and compiled locally on-the-fly by this bash script due to quirks of our specific g++ setup. The filtering code can be viewed within [`process_ngs.sh`](process_ngs.sh) or [`example_quality_filtering.cpp`](example_quality_filtering.cpp) for adaptation into your own pipeline, if desired.
4. For each sample, the pipeline outputs `good_reads.csv` containing all successfully merged sequencing reads passing the quality filters set in [`params.csv`](params.csv). `good_reads.csv` is the input for all barcoding and functional scoring analyses.
