import argparse
import pysam
import pandas as pd

"""
This script processes a BAM file to calculate per-nucleotide error frequencies, 
considering only alignments on the positive strand. The script computes metrics 
such as mismatch, insertion, and deletion frequencies, along with the mean quality score
for each position in the reference sequence.

Note that this script defines Coverage as number of reads spanning a position.
Consequently, the calculation of basecalling error frequencies are based
on the total number of reads spanning each position, whether there is a base
aligned to that exact position or not.

Mismatch, nt frequencies & mean quality scores are calculated based on bases mapped,
not spanning read coverage.

The script requires a BAM file and a corresponding FASTA file as input.
The results are outputted as a TSV file.

Where:
    bam_file: Path to the input BAM file.
    fasta_file: Path to the reference FASTA file.
    output_tsv: Path for the output TSV file containing the error frequencies.

Example:
    python get_bcerror_freqs.py sample.bam reference.fasta output.tsv
"""

def calculate_error_frequencies(bam_file, fasta_file):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    faidx = pysam.FastaFile(fasta_file)

    error_data = []

    for ref in faidx.references:
        ref_len = faidx.get_reference_length(ref)
        coverage = [0] * ref_len
        base_counts = {'A': [0] * ref_len, 'T': [0] * ref_len, 
                       'G': [0] * ref_len, 'C': [0] * ref_len, 'N': [0] * ref_len}
        mismatches = [0] * ref_len
        insertions = [0] * ref_len
        deletions = [0] * ref_len
        quality_scores = [0] * ref_len
        bases_mapped = [0] * ref_len # for mismatch & insertion frequency calcs

        for read in samfile.fetch(ref):
            if read.is_unmapped or read.is_reverse:
                continue

            ref_pos = read.reference_start
            read_pos = 0
            read_seq = read.query_sequence
            read_qual = read.query_qualities

            for cigar_op, cigar_len in read.cigartuples:
                if cigar_op in [0, 7, 8]:  # Matches and mismatches
                    for i in range(cigar_len):
                        coverage[ref_pos + i] += 1
                        ref_base = faidx.fetch(ref, ref_pos + i, ref_pos + i + 1).upper()
                        read_base = read_seq[read_pos + i].upper()
                        qual = read_qual[read_pos + i]

                        quality_scores[ref_pos + i] += qual
                        bases_mapped[ref_pos + i] += 1

                        if read_base in base_counts:
                            base_counts[read_base][ref_pos + i] += 1

                        if read_base != ref_base or cigar_op == 8:
                            mismatches[ref_pos + i] += 1

                    ref_pos += cigar_len
                    read_pos += cigar_len
                elif cigar_op == 1:  # Insertions should not increment coverage
                    insertions[ref_pos] += 1
                    read_pos += cigar_len
                elif cigar_op == 2:  # Deletions
                    for i in range(cigar_len):
                        coverage[ref_pos + i] += 1
                        deletions[ref_pos + i] += 1
                    ref_pos += cigar_len
                elif cigar_op in [3, 4]:  # Skips and soft clips
                    ref_pos += cigar_len if cigar_op == 3 else 0
                    read_pos += cigar_len if cigar_op == 4 else 0
                # Hard clipped bases and padding (5 and 6) are ignored

        for pos in range(ref_len):
            pos_data = {
                "Reference": ref,
                "Position": pos + 1,
                "Spanning_Reads": coverage[pos],
                "Bases_Mapped": bases_mapped[pos],
                "A_Freq": base_counts['A'][pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "T_Freq": base_counts['T'][pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "G_Freq": base_counts['G'][pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "C_Freq": base_counts['C'][pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "N_freq": base_counts['N'][pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "MismatchFreq": mismatches[pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "InsertionFreq": insertions[pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0,
                "DeletionFreq": deletions[pos] / coverage[pos] if coverage[pos] > 0 else 0,
                "BCErrorFreq": (mismatches[pos] + insertions[pos] + deletions[pos]) / coverage[pos] if coverage[pos] > 0 else 0,
                "MeanQual": quality_scores[pos] / bases_mapped[pos] if bases_mapped[pos] > 0 else 0
            }
            error_data.append(pos_data)

    samfile.close()
    faidx.close()

    return pd.DataFrame(error_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Basecalling Error Frequencies")
    parser.add_argument("bam_file", help="Path to the BAM file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_tsv", help="Path for the output TSV file")
    args = parser.parse_args()

    error_freq_df = calculate_error_frequencies(args.bam_file, args.fasta_file)
    error_freq_df.to_csv(args.output_tsv, sep="\t", index=False)