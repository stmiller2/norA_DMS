#!/bin/bash

#User-defined variables here:
#-----------------------------------------------------------------------------
  path="/scratch/stmiller2/ddG_v2"
  prot="9b3m"
  relaxed_pdb="${path}/inputs/9b3m_tr_A_0001.pdb" #Relaxed with prep_inputs.sh
  fasta="${path}/inputs/9b3m_tr_A.fasta" #Comes from cleanpdb3.py
  spanfile="${path}/inputs/9b3m_tr_A.span" #Did you check for errors?
  
  #Safety: Set this to true when you're ready to pull the trigger
  FULL_RUN=false #If false, will only do pos_1 -> Ala (for testing purposes)
  
#-----------------------------------------------------------------------------

#Dev
#-----------------------------------------------------------------------------
  DDG_APP="/usr/local/bin/cartesian_ddg.cxx11threadserialization.linuxgccrelease"
  DOCKER_ROSETTADB="/usr/local/database" #RosettaCommons docker database
  WT_SEQ=$(grep -v "^>" "$fasta" | tr -d '\n') #Get WT sequence from fasta
  PROTEIN_LENGTH=${#WT_SEQ} #Get sequence length
  AAs=("A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y")
#-----------------------------------------------------------------------------

# Function for creating mutfiles
write_ddg_files() {
	local pos=$1
	local mut_aa=$2
	local wt_aa=${WT_SEQ:pos-1:1} #Get wild type residue at position
	local mutation_path="${path}/muts/pos_${pos}/${mut_aa}"
	
	# Write mutfile
	cat > "${mutation_path}/mut_${pos}${mut_aa}.mutfile" <<EOL
total 1
1
$wt_aa $pos $mut_aa
EOL
	echo ${pos},${mut_aa} >> ${path}/muts/list
}

# Write bash script to run cartesian_ddg application
cat > "${path}/scripts/${prot}_ddG.sh" <<EOL
#!/bin/bash
${DDG_APP} -database ${DOCKER_ROSETTADB} -s $(basename "$relaxed_pdb") -score:weights f19_cart_1.5.wts -in:membrane -mp:setup:spanfiles $(basename "$spanfile") -mp:lipids:composition DLPC -has_pore false -ddg:mut_file mut_\${1}\${2}.mutfile -ddg:legacy true -ddg::dump_pdbs false -ddg:frag_nbrs 4 -ddg:optimize_proline true -ddg:cartesian -ddg:bbnbrs 1 -ddg:iterations 3 -force_iterations false -fa_max_dis 9.0 -missing_density_to_jump 
EOL
	
# Write submit file to send ddg calculations to cluster
cat > "${path}/scripts/${prot}_ddG.sub" <<EOL
Universe = docker
docker_image = rosettacommons/rosetta:latest
log = condor.log
error = condor.err
output = ${prot}_ddG.log

arguments = \$(position) \$(identity)
executable = ${prot}_ddG.sh

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
on_exit_remove = ExitCode =?= 0
on_exit_hold = JobRunCount > 5
periodic_release = (JobStatus == 5) && ((CurrentTime - EnteredCurrentStatus) > 300) && (JobRunCount < 10) && (HoldReasonCode =!= 1) && (HoldReasonCode =!= 6) && (HoldReasonCode =!= 12) && (HoldReasonCode =!= 13) && (HoldReasonCode =!= 14) && (HoldReasonCode =!= 21) && (HoldReasonCode =!= 22)

initialdir = ${path}/muts/pos_\$(position)/\$(identity)/
transfer_input_files = ${path}/scripts/${prot}_ddG.sh, ${relaxed_pdb}, ${path}/scripts/f19_cart_1.5.wts, ${spanfile}, ${path}/muts/pos_\$(position)/\$(identity)/mut_\$(position)\$(identity).mutfile

request_cpus = 1
request_memory = 500M
request_disk = 250M

+WantFlocking = true

queue position identity from ${path}/muts/list
EOL

#1) Set up directories and files
if [[ "$FULL_RUN" == true ]]; then
	for pos in $(seq 1 "$PROTEIN_LENGTH"); do
		pos_dir="${path}/muts/pos_${pos}"
		mkdir -p "$pos_dir"
		for aa in "${AAs[@]}"; do
			mkdir -p "${pos_dir}/${aa}"
			write_ddg_files "$pos" "$aa"
		done
	done
else
	#Only run M1A (for testing)
	mkdir -p "${path}/muts/pos_1/A"
	write_ddg_files 1 "A"
fi

#2) Submit jobs
cd ${path}/scripts
condor_submit ${prot}_ddG.sub

