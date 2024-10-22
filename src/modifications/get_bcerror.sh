#! /usr/bin/env bash

#BSUB -J bcerr_loop 
#BSUB -eo logs/bcerr.%J.erreut
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# Define input and output directories
input_dir="/beevol/home/whitel/tRNAworkshop/rebasecalled/alignedbams/full_length/subset"
output_dir="/beevol/home/whitel/tRNAworkshop/bcerror/RNA004v5/hac"
ref="/beevol/home/whitel/tRNAworkshop/ref/tthermophila-numbered-mature-tRNAs.fa"

# Iterate through each BAM file in the input directory
for bam_file in $input_dir/*.bam; do
    base_name=$(basename "$bam_file" .bwa_full_length.bam)
    output_file="$output_dir/${base_name}.tsv"
    /beevol/home/whitel/tRNAworkshop/src/get_bcerror_freqs.py "$bam_file" "$ref" "$output_file"

 done


