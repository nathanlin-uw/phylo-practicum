#!/bin/bash

# Change these each time
DATADIR="../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes"
run_description="all_genes"

mkdir ../results/
mkdir ../results/RAxML/
mkdir ../results/RAxML/gene_trees

# We'll store the ones we already did here
completed_filenames="../results/RAxML/gene_trees/${run_description}_log_completed_files.txt"
touch "$completed_filenames"

for file in "$DATADIR"/*; do
	# Check to see if we have run RAxML on the file already 
		# Exit status of 0 means true, and if there is a match the exit status will be 0
	if grep "$file" "$completed_filenames"; then
		# Do nothing
		echo "Already finished ${file}"
	else
		# Run RAxML on the file and putting the results in our results folder
		raxml-ng --msa $file --model GTR+G4
    		mv ${file}**raxml** ../results/RAxML/gene_trees
		# Add the file to the completed filenames record
		echo "$file" >> "$completed_filenames"
	fi
done
