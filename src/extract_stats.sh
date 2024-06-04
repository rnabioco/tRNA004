#!/bin/bash

# Output file
output_file="crossalignments_samtools_stats_summary.txt"

#output locatoin
outloc="/beevol/home/whitel/tRNAworkshop/metadata"

# Header for the output file
echo -e "Filename\tReadsMapped\tErrorRate\tAverageLength\tAverageQuality" > "$outloc/$output_file"

# Directory containing your BAM files
bam_dir="/beevol/home/whitel/tRNAworkshop/rebasecalled/alignedbams/crossalignments"

# Loop through each BAM file in the directory
for bam_file in "$bam_dir"/*.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "Skipping, not a file: $bam_file"
        continue
    fi
    
    # Extract filename without path
    filename=$(basename "$bam_file")
    
    # Extract statistics using samtools
    stats=$(samtools stats "$bam_file" | head -n 40)
    
    # Parse specific statistics, refining to ignore anything after a '#' and trimming whitespace
    reads_mapped=$(echo "$stats" | awk -F":\t" '/reads mapped:/ {print $2}' | cut -d' ' -f1)
    error_rate=$(echo "$stats" | awk -F":\t" '/error rate:/ {print $2}' | awk '{print $1}')
    average_length=$(echo "$stats" | awk -F":\t" '/average length:/ {print $2}' | awk '{print $1}')
    average_quality=$(echo "$stats" | awk -F":\t" '/average quality:/ {print $2}' | awk '{print $1}')
    
    # Append the extracted statistics to the output file
    echo -e "$filename\t$reads_mapped\t$error_rate\t$average_length\t$average_quality" >> "$outloc/$output_file"
done
