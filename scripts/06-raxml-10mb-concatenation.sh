#!/bin/bash

DATADIR="../data/Wheat_Relative_History_Data_Glemin_et_al/Concatenation10Mb_OneCopyGenes"

mkdir ../results/
mkdir ../results/RAxML/
mkdir ../results/RAxML/10Mb-concatenation

# Store the ones we already did here
completed_filenames="../results/RAxML/10Mb-concatenation/concat_windows_log_completed_files.txt"
touch "$completed_filenames"

for file in "$DATADIR"/*; do
	# Check if RAxML was already run on the file
	if grep "$file" "$completed_filenames"; then
		echo "Already did ${file}"
	else
		raxml-ng --msa $file --model GTR+G4 --threads 5
    		mv ${file}**raxml** ../results/RAxML/10Mb-concatenation
		echo "$file" >> "$completed_filenames"
	fi
done
