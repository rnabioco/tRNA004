#!/bin/bash
#BSUB -J count_bam_reads
#BSUB -o count_bam_reads_out.%J
#BSUB -e count_bam_reads_err.%J
#BSUB -n 1
#BSUB -q rna

# Output file
output_file="bam_read_counts.txt"

# Output location
outloc="/beevol/home/whitel/tRNAworkshop/metadata"

# Ensure the output directory exists
mkdir -p "$outloc"

# Header for the output file
echo -e "Filename\tTotal Reads" > "$outloc/$output_file"

# Directory containing your BAM files
bam_dir="/beevol/home/whitel/tRNAworkshop/rebasecalled/alignedbams/allreads"

# Loop through each BAM file in the directory
for bam_file in "$bam_dir"/*.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "Skipping, not a file: $bam_file"
        continue
    fi
    
    # Extract filename without path
    filename=$(basename "$bam_file")
    
    # Count the number of reads in the BAM file
    total_reads=$(samtools view -c "$bam_file")
    
    # Append the filename and total reads to the output file
    echo -e "$filename\t$total_reads" >> "$outloc/$output_file"
done

