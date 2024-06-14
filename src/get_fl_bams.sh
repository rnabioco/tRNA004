#! /usr/bin/env bash

#BSUB -J bam_filter 
#BSUB -eo logs/bam_filter.%J.errout
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# Define input and output directories
input_dir="/beevol/home/whitel/tRNAworkshop/rebasecalled/alignedbams"
output_dir="/beevol/home/whitel/tRNAworkshop/rebasecalled/alignedbams/full_length/phage_subset"

# Iterate through each BAM file in the input directory
for bam_file in $input_dir/*.bam; do
    base_name=$(basename "$bam_file" .bam)
    output_file="$output_dir/${base_name}_full_length.bam"
    /beevol/home/whitel/tRNAworkshop/src/filter_reads.py "$bam_file" "$output_file"

    samtools index $output_file
 done


