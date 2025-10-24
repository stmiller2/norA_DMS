#!/bin/bash

#User-defined variables here:
#-----------------------------------------------------------------------------
  path="/scratch/stmiller2/ddG_v2"
  prot="9b3m"
  pdb="${path}/inputs/9b3m_tr.pdb" #PDB file (transformed into membrane!) Use https://opm.phar.umich.edu/
  chain="A"
  auto_spanfile=true #generate spanfile from pdb? (set SPAN_PATH in Dev if false)
  auto_cleanpdb=true #run clean_pdb3.py (set CLEANEDPDB_PATH in Dev if false)
#-----------------------------------------------------------------------------

#Dev:
#-----------------------------------------------------------------------------
  LOCAL_ROSETTADB="/scratch/ameger/software/Rosetta/main/database/"
  DOCKER_ROSETTADB="/usr/local/database"
  ROSETTASCRIPTS="/usr/local/bin/rosetta_scripts.cxx11threadserialization.linuxgccrelease"
  SPAN_PATH="${path}/inputs/9b3m_tr_A.span" #only necessary if auto_spanfile=false
  CLEANEDPDB_PATH="${path}/inputs/9b3m_tr_A.pdb" #only necessary if auto_cleanpdb=false
#-----------------------------------------------------------------------------

#Prepare pdb and spanfile
#-----------------------------------------------------------------------------
cd "${path}/inputs"

#Clean pdb
if $auto_cleanpdb; then
	${path}/scripts/clean_pdb3.py ${pdb} ${chain}
	cleaned_pdb="$(basename "$pdb" .pdb)_${chain}.pdb"
	echo "PDB cleaned"
else
	cleaned_pdb="$(basename "$CLEANEDPDB_PATH")"
	echo "Using ${CLEANEDPDB_PATH} for cleaned pdb"
fi

#Make spanfile
if $auto_spanfile; then
	${path}/scripts/spanfile_from_pdb.linuxgccrelease -database ${LOCAL_ROSETTADB} -in:file:s ${cleaned_pdb}
	spanfile="$(basename "$cleaned_pdb" .pdb).span"
	echo "Spanfile created, remember to check for errors!"
else
	spanfile="$(basename "$SPAN_PATH")"
	echo "Using ${SPAN_PATH} for spanfile"
fi
#-----------------------------------------------------------------------------

#Cartesian MP relax 
#-----------------------------------------------------------------------------

#Build executable and submit files to send relax to cluster
cat > "${path}/inputs/relax_${prot}.sh" <<EOL
#!/bin/bash
${ROSETTASCRIPTS} -database ${DOCKER_ROSETTADB} -parser:protocol mp_cart_relax.xml -parser:script_vars repeats=5 energy_func=f19_cart_1.5.wts energy_fawtb=0 -in:file:s ${cleaned_pdb} -optimization::default_max_cycles 200 -mp:setup:spanfiles ${spanfile} -mp:scoring:hbond -mp:lipids:composition DLPC -mp::thickness 15 -relax:jump_move true -relax:coord_constrain_sidechains -relax:constrain_relax_to_start_coords -nstruct 1 -fa_max_dis 9.0 -ignore_unrecognized_res true -packing:pack_missing_sidechains false -ex1 -ex2 -flip_HNQ -missing_density_to_jump -score:weights f19_cart_1.5.wts -out:pdb -out:file:scorefile ${prot}_relax_scores.sc
EOL

cat > "${path}/inputs/relax_${prot}.sub" <<EOL
Universe = docker
docker_image = rosettacommons/rosetta:latest
log = condor.log
error = condor.err
output = relax_${prot}.log

executable = relax_${prot}.sh

should_transfer_files = YES
when_to_transfer_output = ON_EXIT

transfer_input_files = ${path}/inputs/relax_${prot}.sh, ${path}/inputs/${cleaned_pdb}, ${path}/scripts/f19_cart_1.5.wts, ${path}/inputs/${spanfile}, ${path}/scripts/mp_cart_relax.xml

request_cpus = 8
request_memory = 4G
request_disk = 4G

queue 1
EOL

#Submit file
cd ${path}/inputs
condor_submit relax_${prot}.sub

#Command for running locally:
#$/scratch/ameger/software/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -database ${LOCAL_ROSETTADB} -parser:protocol ${path}/scripts/mp_cart_relax.xml -parser:script_vars repeats=5 energy_func=${path}/scripts/f19_cart_1.5.wts energy_fawtb=0 -in:file:s ${path}/inputs/${cleaned_pdb} -optimization::default_max_cycles 200 -mp:setup:spanfiles ${path}/inputs/${spanfile} -mp:scoring:hbond -mp:lipids:composition DLPC -mp::thickness 15 -relax:jump_move true -relax:coord_constrain_sidechains -relax:constrain_relax_to_start_coords -nstruct 1 -fa_max_dis 9.0 -ignore_unrecognized_res true -packing:pack_missing_sidechains false -ex1 -ex2 -flip_HNQ -missing_density_to_jump -score:weights ${path}/scripts/f19_cart_1.5.wts -out:pdb -out:file:scorefile ${path}/inputs/relax_scores.sc

#Some notes on the above parameters:
#-----------------------------------------------------------------------------
#Cartesian-ddG application docs say you need to relax in cartesian space, so most things
#here are from the GitHub linked in the Tiemann 2023 Biophysical Journal paper. The
#fa19_cart_1.5.wts file is also from this paper (cartesian implementation of the 
#franklin2019 scorefunction, I guess). The one change I made is adjusting the max cycles 
#to 200, down from default 2000, because someone on RosettaCommons forums suggested it 
#saying it was good enough and would reduce the frequency of the "Inaccurate G!" error.
#
#Another note: rmoretti says the "FUNC NOT CONVERGED" warning is not something to be 
#concerned about, because the minimizer is applied repeatedly, so even if it's not 
#converged on that particular iteration, it probably got pretty close, and it'll have an
#opportunity to get closer on the next round when the minimizer is applied again.
#-----------------------------------------------------------------------------
