"""
Clean up tRNA .afa files downloaded from modomics
(e.g., https://genesilico.pl/modomics/rnafamilies/rf00005/)
so that they can be easily cross-referenced to the reference files
we use for aligning nanopore tRNA sequencing reads

example usage:

python clean_afas.py path/to/input.fasta path/to/output.fasta

"""
import argparse

# clean up afa headers so that they match sequence alignment reference formats
def clean_fasta_header(input_file, intermediate_file):
    with open(input_file, 'r') as infile, open(intermediate_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.split('|')
                type_ = parts[1]
                amino_acid = parts[2]
                codon = parts[3]
                location = parts[5].strip()

                # Rename initiator methionines
                if amino_acid == "Ini":
                    amino_acid = "iMet"

                if 'cytosol' in location:
                    new_header = f'>nuc-{type_}-{amino_acid}-{codon}\n'
                elif 'mitochondrion' in location:
                    new_header = f'>mito-{type_}-{amino_acid}-{codon}\n'
                else:
                    new_header = line

                outfile.write(new_header)
            else:
                outfile.write(line)

# remove entries with missing info; add tRNA splint adapter sequences to all entries
def clean_and_modify_sequences(intermediate_file, output_file):
    with open(intermediate_file, 'r') as infile, open(output_file, 'w') as outfile:
        skip_next_line = False
        for line in infile:
            if line.startswith('>'):
                # Check for missing amino acid or anticodon
                parts = line.split('-')
                if parts[2] == '' or parts[2] == 'None' or parts[3] == '\n':
                    skip_next_line = True
                    continue
                else:
                    skip_next_line = False
                    outfile.write(line)
            else:
                if not skip_next_line:
                    # Replace U with T and append/prepend adapter sequences
                    modified_sequence = 'CCTAAGAGCAAGAAGAAGCCTGGN' + line.strip().replace('U', 'T') + 'GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA\n'
                    outfile.write(modified_sequence)

def main():
    parser = argparse.ArgumentParser(description='Clean and modify FASTA file sequences.')
    parser.add_argument('input_file', type=str, help='Input FASTA file path')
    parser.add_argument('output_file', type=str, help='Output FASTA file path')

    args = parser.parse_args()

    intermediate_file = 'intermediate.fasta'
    clean_fasta_header(args.input_file, intermediate_file)
    clean_and_modify_sequences(intermediate_file, args.output_file)

if __name__ == "__main__":
    main()