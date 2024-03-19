#!/usr/bin/env bash

in="$HOME/tRNAworkshop/rebasecalled/alignedbams/full_length"
dest="$HOME/tRNAworkshop/bedgraphs/full_length"

# Ensure the output directory exists
mkdir -p "$dest"

set -o nounset -o pipefail -o errexit -x

# Iterate through each BAM file in the input directory
for bam_file in $in/*.bam; do
    # Extract the base name without the path and extension
    base_name=$(basename "$bam_file" .bwa.bam)

    # Calculate temporary scale factor for normalization by CPM
    tmpscale=$(samtools view -c "$bam_file" | awk '{print 1000000/$1}')

    # Generate scaled bedgraph 
    bedtools genomecov -ibam "$bam_file" -scale $tmpscale -bg > "${dest}/${base_name}.fl.bg"
done

