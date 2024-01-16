import csv
import argparse

"""
Test script for generating a TSV that maps positions between a FASTA file
and a structurally annotated FSA (also called AFA) file. This script assumes that multiple FASTA
entries can be mapped to a single AFA entry, allowing us to annotate multiple
tRNA isodecoders with their structurally-indexed positions.

The output TSV file has the following columns:
    tRNA: The name of the tRNA entry in the FSA file
    fsa: The nucleotide at the given position in the FSA file
    (fsa entries can include '-' for gaps or other characters for modified nucleotides)
    structure_pos: The position of the nucleotide in the FSA file
    sequence_pos: The position of the nucleotide in the FASTA file

FASTA reference names and FSA entry names need to be sanity checked 
for this script to work. Recommended format for FSA:
>mito-tRNA-Phe-GAA
>nuc-tRNA-Phe-GAA

If FASTA entries don't have a prefix, all tRNAs will be assumed to be from the nuclear genome.
FASTA file names need to match at the start, but can contain additional characters after the
anticodon sequence (e.g., mito-tRNA-Phe-GAA-1, nuc-tRNA-Phe-GAA-2-1, etc.)

FSA files should be processed with clean_afas.py for naming consistency between the two files.
Afterwards, names will contain duplicated anticodons to retain special characters
(e.g., mito-tRNA-Val-AAC-IAC) but still permit matching to the FASTA file.

Example usage:
python seq_to_structure_mappings.py path/to/reference.fasta path/to/annotated.fsa path/to/output.tsv
"""

def read_fasta(file_name):
    """Read a FASTA or FSA file and return a dictionary of sequences with full headers, and the headers separately."""
    sequences = {}
    headers = []
    with open(file_name, 'r') as file:
        seq_id = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:]  # Use the full header as the key
                headers.append(seq_id)
                sequences[seq_id] = ''
            else:
                sequences[seq_id] += line.upper()
    return sequences, headers

def create_mapping(ref_seq, ann_seq):
    """Create a mapping from reference to annotated sequence."""
    mapping = {}
    ref_index = 0
    for ann_index, ann_nuc in enumerate(ann_seq):
        if ann_nuc != '-':  # If not a gap
            mapping[ann_index + 1] = ref_index + 1
            ref_index += 1
        else:
            # If it's a gap, there's no corresponding sequence position
            mapping[ann_index + 1] = None
    return mapping

def main(fasta_file, fsa_file, output_file):
    # Load sequences and headers from files
    ref_seqs, ref_headers = read_fasta(fasta_file)
    ann_seqs, ann_headers = read_fasta(fsa_file)

    # Write to TSV
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['seq_ref', 'struct_ref', 'struct_nt', 'struct_pos', 'seq_pos'])

        for ann_id in ann_seqs:
            matching_refs = []
            if ann_id.startswith('mito'):
                # Logic for mitochondrial tRNAs
                parts = ann_id.split('-')
                amino_acid = parts[2]
                anticodon_fsa = parts[3].replace('T', 'U')
                matching_refs = [ref_id for ref_id in ref_seqs if ref_id.startswith(f'mito-tRNA-{amino_acid}-{anticodon_fsa}')]
            elif ann_id.startswith('nuc'):
                # Logic for nuclear tRNAs
                parts = ann_id.split('-')
                amino_acid = parts[2]
                anticodon_fsa = parts[3].replace('T', 'U')
                matching_refs = [ref_id for ref_id in ref_seqs if ref_id.startswith(f'nuc-tRNA-{amino_acid}-{anticodon_fsa}')]
            else:
                # Flexible matching for nuclear tRNAs
                # Extract key parts of the header for comparison
                ann_key_parts = ann_id.split('-')[1:4]  # e.g., ['tRNA', 'Tyr', 'GTA']
                for ref_id in ref_seqs:
                    ref_key_parts = ref_id.split('-')[:3]  # e.g., ['tRNA', 'Tyr', 'GTA']
                    if ann_key_parts == ref_key_parts:
                        matching_refs.append(ref_id)

            for ref_id in matching_refs:
                mapping = create_mapping(ref_seqs[ref_id], ann_seqs[ann_id])
                for ann_index, ann_nuc in enumerate(ann_seqs[ann_id]):
                    fasta_pos = mapping.get(ann_index + 1, '')
                    writer.writerow([ref_id, ann_id, ann_nuc, ann_index + 1, fasta_pos])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map positions from a FASTA file to a structurally annotated FSA file.")
    parser.add_argument("fasta_file", help="Path to the FASTA reference file")
    parser.add_argument("fsa_file", help="Path to the annotated structure file")
    parser.add_argument("output_file", help="Path to the output TSV file")
    args = parser.parse_args()

    main(args.fasta_file, args.fsa_file, args.output_file)
