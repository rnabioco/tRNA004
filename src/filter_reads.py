#!/usr/bin/env python

import pysam

# Open the BAM file
bamfile = pysam.AlignmentFile("your_file.bam", "rb")

# Adapter lengths
adapter_5p_len = 24  # Total length of the 5' tRNA adapter
adapter_3p_len = 30  # Total length of the 3' tRNA adapter
truncation_max = 15  # Maximum truncation allowed at the 5' end (as helicase loses contact with the RNA)

# Filtering logic to ID "full length" reads that align to adapters but permit 5' end truncation
with pysam.AlignmentFile("filtered_alignments.bam", "wb", template=bamfile) as outfile:
    for read in bamfile:
        # Assuming the reference sequence includes the 5' adapter, RNA body, and 3' adapter,
        # and that the adapters are at the very beginning and end of the reference:
        
        # Calculate the expected minimum alignment start (considering truncation)
        expected_min_start = 0  # Start of the reference (5' adapter) allowing for some truncation

        # Calculate the expected minimum end position considering the 3' adapter
        # (read should align through or very close to the end of the RNA body, into the 3' adapter)
        ref_length = bamfile.get_reference_length(read.reference_name)
        expected_min_end = ref_length - adapter_3p_len  # Read should extend into the 3' adapter
        
        # Check if the read meets the criteria
        if read.reference_start <= expected_min_start + truncation_max and read.reference_end >= expected_min_end:
            outfile.write(read)

bamfile.close()
