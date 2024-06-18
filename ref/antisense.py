#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def generate_sense_and_antisense_fasta(input_fasta, output_handle):
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Write sense sequence
        sense_record = record
        sense_record.id = "sense-" + record.id
        sense_record.description = "sense-" + record.description
        SeqIO.write(sense_record, output_handle, "fasta")

        # Write antisense sequence
        antisense_seq = record.seq.reverse_complement()
        antisense_record = record
        antisense_record.id = "antisense-" + record.id
        antisense_record.description = "antisense-" + record.description
        antisense_record.seq = antisense_seq
        SeqIO.write(antisense_record, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("Usage: python antisense.py <input_fasta> [<output_fasta>]")
        sys.exit(1)

    input_fasta = sys.argv[1]

    if len(sys.argv) == 3:
        output_fasta = sys.argv[2]
        with open(output_fasta, 'w') as output_handle:
            generate_sense_and_antisense_fasta(input_fasta, output_handle)
    else:
        generate_sense_and_antisense_fasta(input_fasta, sys.stdout)

