"""
Clean up the headers of tRNA .afa files downloaded from modomics
(e.g., https://genesilico.pl/modomics/rnafamilies/rf00005/)
so that they can be easily joined to sequence reference files

example usage:

python clean_afas.py path/to/input.fasta path/to/output.fasta

"""
import argparse

def clean_fasta_header(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.split('|')
                type_ = parts[1]
                amino_acid = parts[2]
                codon = parts[3]
                location = parts[5].strip()

                if 'cytosol' in location:
                    new_header = f'>nuc-{type_}-{amino_acid}-{codon}\n'
                elif 'mitochondrion' in location:
                    new_header = f'>mito-{type_}-{amino_acid}-{codon}\n'
                else:
                    new_header = line

                outfile.write(new_header)
            else:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description='Clean FASTA file headers.')
    parser.add_argument('input_file', type=str, help='Input FASTA file path')
    parser.add_argument('output_file', type=str, help='Output FASTA file path')

    args = parser.parse_args()

    clean_fasta_header(args.input_file, args.output_file)

if __name__ == "__main__":
    main()