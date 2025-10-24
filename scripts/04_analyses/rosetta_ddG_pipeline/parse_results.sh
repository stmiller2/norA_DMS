#!/bin/bash

#User-defined variables here:
#-----------------------------------------------------------------------------
  path="/scratch/stmiller2/ddG_v2"
  output_csv="${path}/ddG_results.csv"
#-----------------------------------------------------------------------------

#Combine all ddG values into a single csv file:
#-----------------------------------------------------------------------------
echo "pos,aa,wtavg,mutavg,ddG" > "$output_csv" #Initialize CSV with header

#Loop through all ddg files in the muts directory
for ddg_file in $(find "${path}/muts" -name "*.ddg"); do
    # Extract position and mutation from the filename
    filename=$(basename "$ddg_file")
    pos=$(echo "$filename" | sed -E 's/mut_([0-9]+)([A-Z]+).ddg/\1/')
    aa=$(echo "$filename" | sed -E 's/mut_([0-9]+)([A-Z]+).ddg/\2/')

    # Extract average WT and MUT total energies
    wt_avg=$(grep "COMPLEX:   Round" "$ddg_file" | grep " WT" | awk '{sum+=$4; count++} END {print sum/count}')
    mut_avg=$(grep "COMPLEX:   Round" "$ddg_file" | grep " MUT" | awk '{sum+=$4; count++} END {print sum/count}')
    # Compute ddG
    ddG=$(echo "scale=6; $mut_avg - $wt_avg" | bc)
    # Save energies to output csv
    echo "$pos,$aa,$wt_avg,$mut_avg,$ddG" >> "$output_csv"
done
#-----------------------------------------------------------------------------

echo "ddG calculations complete. Results saved in $output_csv"