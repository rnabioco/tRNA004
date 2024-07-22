#!/usr/bin/env bash

#BSUB -J calculate_metrics 
#BSUB -eo logs/calculate_metrics.%J.errout
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# used to get alignment metrics at per read level for subsampling for figure s2

# Define input directory, score threshold, and output file
#input_dir="$HOME/tRNAworkshop/rebasecalled/alignedbams/full_length/HF"
#output_csv="$HOME/tRNAworkshop/rebasecalled/metadata/full_length"

input_dir="$HOME/tRNAworkshop/rebasecalled/alignedbams/allreads/HF"
output_csv="$HOME/tRNAworkshop/rebasecalled/metadata/allreads/human.csv"


# Run the Python script
python3 calculate_alignment_metrics.py "$input_dir" "$output_csv"

