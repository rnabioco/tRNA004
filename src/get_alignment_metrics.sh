#!/usr/bin/env bash

#BSUB -J calculate_metrics 
#BSUB -eo logs/calculate_metrics.%J.errout
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# Define input directory, score threshold, and output file
input_dir="$HOME/tRNAworkshop/rebasecalled/alignedbams/allreads"
output_csv="$HOME/tRNAworkshop/rebasecalled/alignedbams/allreads/alignmentmetrics.csv"

# Run the Python script
python3 calculate_alignment_metrics.py "$input_dir" "$output_csv"

