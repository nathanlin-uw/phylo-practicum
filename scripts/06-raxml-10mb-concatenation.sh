#!/bin/bash

DATADIR="../data/Wheat_Relative_History_Data_Glemin_et_al/Concatenation10Mb_OneCopyGenes"

mkdir ../results/
mkdir ../results/RAxML/
mkdir ../results/RAxML/10Mb-concatenation

for file in "$DATADIR"/*; do
    raxml-ng --msa $file --model GTR+G4 --threads 5
    mv ${file}**raxml** ../results/RAxML/10Mb-concatenation
done
