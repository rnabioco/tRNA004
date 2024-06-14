#! /usr/bin/env bash

#BSUB -J bam_to_fq
#BSUB -oe logs.%J.errout
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# Loop through all BAM files in the current directory and convert to fastqs
# use case: Dorado rebasecalling produces unaligned bams
for file in ./*.bam; do
  # Get the filename without the extension
  base_name=$(basename "$file" .bam)
  
  # Convert BAM to FASTQ using samtools
  samtools fastq "$file" > "./$base_name.fastq"
done
