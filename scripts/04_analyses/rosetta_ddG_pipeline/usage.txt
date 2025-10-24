--------------------------------------------------------------------------------
    Silas' Membrane Protein ddG Scan Pipeline - Version 2, Updated 02/14/2025
                           Contact: stmiller2@wisc.edu
--------------------------------------------------------------------------------

Overview:
---------
This pipeline facilitates computing Rosetta ΔΔG scores for all possible single
missense mutations of a membrane protein. It takes inspiration and code from 
Tony Meger's soluble protein ddG pipeline, as well as the Tiemann 2023 
Biophysical Journal paper. More detailed notes about what was lifted or changed
can be found in comments throughout the code.
This was written for the Biochemistry Compute Cluster, but should be easy to 
adapt to any system with HTCondor.

Usage Instructions:
-------------------

**PREPARE FILES**
    **STEP 1:** Copy the entire ddG pipeline directory and rename it to reflect 
      your experiment.

    **STEP 2:** Orient your protein in the membrane with the OPM database
      https://opm.phar.umich.edu/ #OPM
      https://opm.phar.umich.edu/ppm_server3_cgopm #PPM 3.0
          - You can do this either by searching your protein in the database
            to see if it's already oriented, or you can upload a PDB file to
            the PPM 3.0 webserver for orientation
          - It's OK if there's multiple chains (e.g., an antibody fragment or
            cryo-EM fiducial). You will extract only the chain of interest later
          - Once oriented, download the transformed file and open it in PyMol
            to confirm that the membrane looks good
          - Copy the transformed pdb prot_tr.pdb to the "inputs" directory of
            the ddG pipeline

    **STEP 3:** Copy the directory over to BCC
          - cp -r ddG_v2 /scratch/{username}/
          
    **STEP 4:** Enter the "prep_inputs.sh" script with your preferred 
      commandline file editor and set the user variables
      
    **STEP 5:** Run prep_inputs.sh to clean the pdb, generate a spanfile, and
      relax the protein with Rosetta.
          - ./prep_inputs.sh
          - Double-check the generated spanfile prot_tr_A.span to be sure it
            makes sense (i.e., correct number of transmembrane domains). Note
            that the numbering in this file is Rosetta numbering, so they may
            not line up exactly with your original PDB (Rosetta numbering skips
            missing density)
          - The relaxation might take like 30 minutes, depending on how large
            your protein is.
          - When the relaxation finishes, I like to copy the relaxed pdb
            prot_tr_A_0001.pdb back to my computer and open it in PyMol just to
            make sure Rosetta didn't do anything stupid. Open the original PDB 
            and relaxed PDB in the same PyMol session, align if necessary, and
            confirm that the relaxation didn't mess with the conformation too 
            much. Generally the backbone should not have moved much at all
            
**RUN DDG SCAN**
    **STEP 6:** Enter the "mp_cartddG_pipeline.sh" script with your preferred
      commandline file editor and set the user variables
          - It's recommended to start with FULL_RUN=false, which will only do
            M1A, to confirm everything is set up correctly. Then, you can set
            FULL_RUN=true and do the entire scan (which takes a long time and
            submits many HTCondor jobs)
            
    **STEP 7:** Run mp_cartddG_pipeline.sh to create a directory and mutfile 
      for each mutation, create a bash script to run the cartesian_ddg Rosetta
      application, and submit HTCondor jobs for each mutation's ddG calculation
          - ./mp_cartddG_pipeline.sh
          - Depending on the size of your protein, this may take several hours.
            You can check the status of your jobs with condor_q, or 
            condor_tail -f {jobid} to see how a specific run is going
            
**PARSE AND SAVE RESULTS**
    **STEP 8:** When all your HTCondor jobs are complete, enter the 
      "parse_results.sh" script with your preferred commandline file editor and
      set the user variables
      
    **STEP 9:** Run parse_results.sh to loop through all the ddG output files,
      take the average WT and MUT total energies from each of the iterations
      cartesian_ddg does, and compute the difference (ddG). All these values are
      saved in the output csv
          - ./parse_results.sh
          - This might take a few minutes depending on the size of your protein
          
    **STEP 10:** Save data and clean up. 
          - Recommended: Copy the csv over to the fileserver for further 
          analysis. Then, tar the entire ddG directory for safekeeping. This has
          all your raw data, scripts, log files, etc.          
            tar -czf ddG_scan.tar.gz ddG_v2/
          Then move the tar file to the fileserver and remove everything from
          the scratch directory
    
    **STEP 11:** Profit????

--------------------------------------------------------------------------------