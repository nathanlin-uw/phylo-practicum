#!/bin/bash

MAPFILE="../results/07-species_mapping.txt"
DATADIR="../results/10concatenation10Mb_OneCopy-phylip-ch3"
OUTDIR="../results/HyDe/10Mb-concatenation-ch3"

mkdir ../results/HyDe/
mkdir ../results/HyDe/10Mb-concatenation-ch3

for file in "$DATADIR"/*; do
	echo "Processing file: $file

    # 1. Extract the filename from the path
    filename=$(basename "$file")
    base="${filename%.*}"


    header=$(head -n 1 "$file")
    n=$(echo $header | awk '{print $1}')
    s=$(echo $header | awk '{print $2}')

    run_hyde.py \
        -i "$file" \
        -m "$MAPFILE" \
        -o H_vulgare \
        -n "$n" \
        -t 17 \
        -s "$s" \
        --prefix "$OUTDIR/$base"

done
