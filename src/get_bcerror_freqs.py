import argparse
import pysam
import pandas as pd

def calculate_error_frequencies(bam_file, fasta_file):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    faidx = pysam.FastaFile(fasta_file)

    error_data = []

    for ref in faidx.references:
        ref_len = faidx.get_reference_length(ref)
        base_calls = [0] * ref_len
        errors = [0] * ref_len

        for read in samfile.fetch(ref):
            if read.is_unmapped:
                continue

            ref_pos = read.reference_start
            read_pos = 0
            read_seq = read.query_sequence

            for cigar_op, cigar_len in read.cigartuples:
                if cigar_op == 0:  # Matches and mismatches
                    for i in range(cigar_len): # differentiate matches from mismatches
                        ref_base = faidx.fetch(ref, ref_pos + i, ref_pos + i + 1)
                        read_base = read_seq[read_pos + i]
                        base_calls[ref_pos + i] += 1
                        if read_base != ref_base:
                            errors[ref_pos + i] += 1
                    ref_pos += cigar_len
                    read_pos += cigar_len
                elif cigar_op == 1:  # Insertions
                    read_pos += cigar_len
                elif cigar_op == 2:  # Deletions
                    ref_pos += cigar_len
                elif cigar_op == 3:  # Skipped regions
                    ref_pos += cigar_len
                elif cigar_op == 4:  # Soft clipped bases
                    read_pos += cigar_len
                elif cigar_op == 7: # in the event of = used as match
                    for i in range(cigar_len):
                        ref_pos += cigar_len
                        read_pos += cigar_len
                elif cigar_op == 8: # in the event of X used as mismatch
                    for i in range(cigar_len):
                        errors[ref_pos + i] += 1
                        base_calls[ref_pos + i] += 1
                    ref_pos += cigar_len
                    read_pos += cigar_len
                else:  # 5 and 6 = Hard clipped bases, padding
                    pass

        # Calculate error frequencies for this reference
        error_freqs = [errors[i] / base_calls[i] if base_calls[i] > 0 else 0 for i in range(ref_len)]
        
        # Store the error frequencies for this reference
        for pos, freq in enumerate(error_freqs):
            error_data.append((ref, pos + 1, freq))

    samfile.close()
    faidx.close()

    # Convert the results to a DataFrame
    err_df = pd.DataFrame(error_data, columns=["Reference", "Position", "ErrorFrequency"])
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