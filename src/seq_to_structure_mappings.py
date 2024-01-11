import csv
import argparse

"""
Test script for generating a TSV that maps positions between a FASTA file
and a structurally annotated FSA file. This script assumes that multiple FASTA
entries can be mapped to a single FSA entry, allowing us to annotate multiple
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

FASTA file names need to match at the start, but can contain additional characters after the
anticodon sequence (e.g., mito-tRNA-Phe-GAA-1, nuc-tRNA-Phe-GAA-2-1, etc.)

Example usage:
python script_name.py path/to/reference.fasta path/to/annotated.fsa path/to/output.tsv
"""

def read_fasta(file_name):
    """Read a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:]  # Remove the '>'
                sequences[seq_id] = ''
            else:
                sequences[seq_id] += line.upper()
    return sequences

def create_mapping(ref_seq, ann_seq):
    """Create a mapping from reference to annotated sequence."""
    mapping = {}
    ref_index = 0
    for ann_index, ann_nuc in enumerate(ann_seq):
        if ann_nuc != '-':  # If not a gap
            mapping[ref_index + 1] = ann_index + 1
            ref_index += 1
    return mapping

def main(fasta_file, fsa_file, output_file):
    # Load sequences from files
    ref_seqs = read_fasta(fasta_file)
    ann_seqs = read_fasta(fsa_file)

    # Create mappings for each annotated sequence
    pos_maps = {}
    for ann_id, ann_seq in ann_seqs.items():
        # Find all matching reference sequences
        matching_refs = [ref_id for ref_id in ref_seqs if ref_id.startswith(ann_id)]
        for ref_id in matching_refs:
            ref_seq = ref_seqs[ref_id]
            pos_maps[ref_id] = create_mapping(ref_seq, ann_seq)

    # Write to TSV
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['tRNA', 'fsa', 'structure_pos', 'sequence_pos'])

        for ann_id, ann_seq in ann_seqs.items():
            matching_refs = [ref_id for ref_id in ref_seqs if ref_id.startswith(ann_id)]
            for ref_id in matching_refs:
                mapping = pos_maps.get(ref_id, {})
                ref_index = 0
                for ann_index, ann_nuc in enumerate(ann_seq):
                    fasta_pos = mapping.get(ref_index + 1, '')
                    writer.writerow([ref_id, ann_nuc, ann_index + 1, fasta_pos])
                    if ann_nuc != '-':
                        ref_index += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map positions from a FASTA file to a structurally annotated FSA file.")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("fsa_file", help="Path to the annotated FSA file")
    parser.add_argument("output_file", help="Path to the output TSV file")
    args = parser.parse_args()

    main(args.fasta_file, args.fsa_file, args.output_file)
