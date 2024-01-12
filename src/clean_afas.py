"""
Clean up tRNA .afa files downloaded from modomics
(e.g., https://genesilico.pl/modomics/rnafamilies/rf00005/)
so that they can be easily cross-referenced to the reference files
we use for aligning nanopore tRNA sequencing reads

example usage:

python clean_afas.py path/to/input.fasta path/to/output.fasta

The output file will contain reference entries in format:
>nuc-tRNA-Val-TAC-TAC
>mito-tRNA-Trp-TCA-!CA

this anticodon duplication preserves special characters from modomics
that indicate isodecoder-specific anticodon modifications (with unique structures)

currently the script only handles special characters present in yeast tRNA anticodons

"""
import argparse
import tempfile

# convert special characters in anticodon to match ref format but retain modomics anticodon
def convert_anticodon(anticodon):
    # Dictionary for special character conversion (yeast only)
    conversion = {'I': 'A', '3': 'T', '!': 'T', '$': 'T', 'N': 'T'}
    
    # Convert special characters
    converted_anticodon = ''.join(conversion.get(char, char) for char in anticodon)
    
    # Return the duplicated anticodon
    return f'{converted_anticodon}-{anticodon}'

# clean up afa headers so that they match sequence alignment reference formats
def clean_fasta_header(input_file, intermediate_file):
    with open(input_file, 'r') as infile, open(intermediate_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.split('|')
                type_ = parts[1]
                amino_acid = parts[2]
                anticodon = parts[3].replace('U', 'T')  # convert U to T
                location = parts[5].strip()

                # Rename initiator methionines
                if amino_acid == "Ini":
                    amino_acid = "iMet"

                # Use convert_anticodon function to deal with special characters
                anticodon_converted = convert_anticodon(anticodon.strip())

                if 'cytosol' in location:
                    new_header = f'>nuc-{type_}-{amino_acid}-{anticodon_converted}\n'
                elif 'mitochondrion' in location:
                    new_header = f'>mito-{type_}-{amino_acid}-{anticodon_converted}\n'
                else:
                    new_header = line

                outfile.write(new_header)
            else:
                outfile.write(line)

# remove entries with missing info
# add tRNA splint adapter sequences to all entries
# deduplicate entries after cleanup
def clean_and_modify_sequences(intermediate_file, output_file):
    unique_entries = set()  # Set to store unique entries
    current_header = None
    skip_next_line = False  # Flag to skip lines with missing info

    with open(intermediate_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Check for missing amino acid or anticodon
                parts = line.split('-')
                if parts[2] == '' or parts[2] == 'None' or parts[3] == '\n':
                    skip_next_line = True
                else:
                    skip_next_line = False
                    current_header = line.strip()  # Store the header for the current entry
            else:
                if not skip_next_line:
                    # Replace U with T and append/prepend adapter sequences
                    modified_sequence = 'CCTAAGAGCAAGAAGAAGCCTGGN' + line.strip().replace('U', 'T') + 'GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA'
                    entry = (current_header, modified_sequence)  # Create a tuple of header and sequence

                    # Write to file only if this entry is unique
                    if entry not in unique_entries:
                        unique_entries.add(entry)  # Add the entry to the set
                        outfile.write(current_header + '\n')  # Write the header
                        outfile.write(modified_sequence + '\n')  # Write the sequence

def main():
    parser = argparse.ArgumentParser(description='Clean and modify FASTA file sequences.')
    parser.add_argument('input_file', type=str, help='Input FASTA file path')
    parser.add_argument('output_file', type=str, help='Output FASTA file path')

    args = parser.parse_args()

    # Use a temporary file
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as temp_file:
        clean_fasta_header(args.input_file, temp_file.name)
        temp_file.seek(0)  # Reset file pointer to the beginning of the file
        clean_and_modify_sequences(temp_file.name, args.output_file)

if __name__ == "__main__":
    main()