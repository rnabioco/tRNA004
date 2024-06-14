#!/bin/bash

# Output file
output_file="fastq_read_counts.txt"

# Output location
outloc="/beevol/home/whitel/tRNAworkshop/metadata"

# Ensure the output directory exists
mkdir -p "$outloc"

# Header for the output file
echo -e "Filename\tTotal Reads\tFL Reads" > "$outloc/$output_file"

# Directory containing your FASTQ files
fastq_dir="/beevol/home/whitel/tRNAworkshop/rebasecalled/fastqs"

# Loop through each FASTQ file in the directory
for fastq_file in "$fastq_dir"/*.fastq; do
    if [ ! -f "$fastq_file" ]; then
        echo "Skipping, not a file: $fastq_file"
        continue
    fi

    # Extract filename without path
    filename=$(basename "$fastq_file")

    # Count the number of lines and divide by 4 to get the number of total reads
    total_reads=$(wc -l < "$fastq_file" | awk '{print $1/4}')

    # Count the number of reads that are 105 nt or longer
    fl_reads=$(awk '{if(NR%4==2) print length($0)}' "$fastq_file" | awk '$1 >= 105' | wc -l)

    # Append the filename, total reads, and fl reads to the output file
    echo -e "$filename\t$total_reads\t$fl_reads" >> "$outloc/$output_file"
done

