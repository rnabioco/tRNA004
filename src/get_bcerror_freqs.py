import argparse
import pysam
import pandas as pd

"""
This script processes a BAM file to calculate per-nucleotide error frequencies, 
considering only alignments on the positive strand. The script computes metrics 
such as mismatch, insertion, and deletion frequencies, along with the mean quality score
for each position in the reference sequence.

The script requires a BAM file and a corresponding FASTA file as input.
The results are outputted as a TSV file.

Usage:
    python get_bcerror_freqs.py [bam_file] [fasta_file] [output_tsv]

Where:
    bam_file: Path to the input BAM file.
    fasta_file: Path to the reference FASTA file.
    output_tsv: Path for the output TSV file containing the error frequencies.

Example:
    python get_bcerror_freqs.py sample.bam reference.fasta output.tsv

Requires: pysam and pandas python libraries.
"""


def calculate_error_frequencies(bam_file, fasta_file):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    faidx = pysam.FastaFile(fasta_file)

    error_data = []

    for ref in faidx.references:
        ref_len = faidx.get_reference_length(ref)
        # initialize arrays to store this data with as many 0s as the length of the ref
        coverage = [0] * ref_len
        base_counts = {'A': [0] * ref_len, 'T': [0] * ref_len, 
                       'G': [0] * ref_len, 'C': [0] * ref_len, 'N': [0] * ref_len}
        mismatches = [0] * ref_len
        insertions = [0] * ref_len
        deletions = [0] * ref_len
        quality_scores = [0] * ref_len # get the Phreds
        quality_counts = [0] * ref_len  # for averaging later


        for read in samfile.fetch(ref):
            if read.is_unmapped:
                continue

            # skip reads aligned to tRNA antisense
            if read.is_reverse:
                continue

            ref_pos = read.reference_start
            read_pos = 0
            read_seq = read.query_sequence
            read_qual = read.query_qualities 

            for cigar_op, cigar_len in read.cigartuples:
            # integers fed to read.cigartuples are defined by pysam
            # cigar_len is the length of each CIGAR operation
                if cigar_op == 0 or cigar_op == 7 or cigar_op == 8:  # handle all match/mismatch options, including X and =
                    for i in range(cigar_len): # differentiate matches from mismatches
                        coverage[ref_pos + i] += 1
                        ref_base = faidx.fetch(ref, ref_pos + i, ref_pos + i + 1)
                        read_base = read_seq[read_pos + i]
                        qual = read_qual[read_pos + i]

                        # update quality scores
                        quality_scores[ref_pos + i] += qual
                        quality_counts[ref_pos + i] += 1

                        # update base counts
                        if read_base in base_counts:
                            base_counts[read_base][ref_pos + i] += 1

                        # update mismatches
                        if read_base != ref_base or cigar_op == 8:
                            mismatches[ref_pos + i] += 1

                    # increment positions
                    ref_pos += cigar_len
                    read_pos += cigar_len
                elif cigar_op == 1:  # Insertions
                    insertions[ref_pos] += cigar_len
                    read_pos += cigar_len
                elif cigar_op == 2:  # Deletions
                    deletions[ref_pos:ref_pos + cigar_len] = [x + 1 for x in deletions[ref_pos:ref_pos + cigar_len]]
                    ref_pos += cigar_len
                elif cigar_op == 3: # Skips
                    ref_pos += cigar_len
                elif cigar_op == 4:  #  Soft clips
                    read_pos += cigar_len
                else:  # 5 and 6 = Hard clipped bases, padding
                    pass

        # Calculate error frequencies and other metrics
        for pos in range(ref_len):
            pos_data = {
                "Reference": ref,
                "Position": pos + 1,
                "Coverage": coverage[pos],
                "A_Freq": base_counts['A'][pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "T_Freq": base_counts['T'][pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "G_Freq": base_counts['G'][pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "C_Freq": base_counts['C'][pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "N_freq": base_counts['N'][pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "MismatchFreq": mismatches[pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "InsertionFreq": insertions[pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "DeletionFreq": deletions[pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "MeanQual": quality_scores[pos] / quality_counts[pos] if quality_counts[pos] > 0 else 0

            }
            error_data.append(pos_data)

    samfile.close()
    faidx.close()

    # Convert the results to a DataFrame
    err_df = pd.DataFrame(error_data)
    return err_df

if __name__ == "__main__":
    # Setup for command-line argument parsing
    parser = argparse.ArgumentParser(description="Calculate Basecalling Error Frequencies")
    parser.add_argument("bam_file", help="Path to the BAM file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_tsv", help="Path for the output TSV file")
    args = parser.parse_args()

    # Call the function with provided arguments
    error_freq_df = calculate_error_frequencies(args.bam_file, args.fasta_file)

    # Save the DataFrame to a TSV file
    error_freq_df.to_csv(args.output_tsv, sep="\t", index=False)