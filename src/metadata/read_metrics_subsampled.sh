#!/bin/bash

# used to generate subsampled data for Fig S2

bam_dir="$HOME/tRNAworkshop/rebasecalled/alignedbams/full_length"
output_dir="$HOME/tRNAworkshop/alignedbams/full_length"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Number of reads to sample (e.g., 10000)
sample_size=10000

# Loop through each BAM file in the directory
for bam_file in "$bam_dir"/*.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "Skipping, not a file: $bam_file"
        continue
    fi
    
    # Extract filename without path
    filename=$(basename "$bam_file" .bam)
    
    # Output file for extracted data
    output_file="$output_dir/${filename}_alignment_data.csv"
    
    # Subsample reads using samtools and awk, then extract alignment scores, MAPQ, read lengths, and average quality scores
    samtools view -s 0.1 "$bam_file" | awk -v sample_size="$sample_size" '
    BEGIN {
        srand()
    }
    {
        if (rand() <= sample_size/NR) {
            as = ""
            for (i = 12; i <= NF; i++) {
                if (substr($i, 1, 2) == "AS") {
                    as = substr($i, 6)
                }
            }
            mapq = $5
            read_length = length($10)
            qual_string = $11
            qual_sum = 0
            for (j = 1; j <= length(qual_string); j++) {
                qual_sum += ord(substr(qual_string, j, 1)) - 33
            }
            avg_qual = qual_sum / read_length
            print as "," mapq "," read_length "," avg_qual
        }
    }
    function ord(c) {
        return c ? c : 0
    }' > "$output_file"
    
    echo "Extracted data for $filename to $output_file"
done

