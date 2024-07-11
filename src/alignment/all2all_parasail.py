import parasail as ps

def read_sequences(filepath):
    with open(filepath, 'r') as file:
        header = None
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if header:
                    yield header, sequence
                header = line.strip()
                sequence = ''
            else:
                sequence += line.strip()
        if header:
            yield header, sequence

def get_alignment_score(sequences, output_filepath):
    # Open output file for all alignments
    with open(output_filepath, 'w') as output_file:
        # Compare each sequence against each other
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                seq1_header, seq1 = sequences[i]
                seq2_header, seq2 = sequences[j]
                
                # parasail alignment
                alignment = ps.sw_trace_striped_32(seq1, seq2, 10, 1, ps.blosum62)
                stats = ps.sw_stats_striped_32(seq1, seq2, 10, 1, ps.blosum62)
                
                # Calculate sequence percent identity
                matches = stats.matches
                length = stats.length
                identity = (matches / length) * 100 if length > 0 else 0

                # write alignment to txt file
                output_file.write(f"Comparing {seq1_header} to {seq2_header}\n")
                output_file.write("Query Alignment: " + alignment.traceback.query + "\n")
                output_file.write("Match Alignment: " + alignment.traceback.comp + "\n")
                output_file.write("Reference Alignment: " + alignment.traceback.ref + "\n")
                output_file.write("Alignment score: " + str(alignment.score) + "\n")
                output_file.write("Sequence Identity: " + str(identity) + "\n")
                output_file.write("-------------------------------------------\n")

# File path
seqfile = "/Users/jillbilodeaux/devel/projects/parasail_alignments/ref/noAdapter_t5_ecoli.fa"

# Read sequences into a list
sequences = list(read_sequences(seqfile))

# Define output file path
output_filepath = 'ecoli-phage-alignment.txt'

# Call the function to get alignment scores and write to output file
get_alignment_score(sequences, output_filepath)

# Inform user of completion and file writing
print(f"tRNA alignment output has been written to '{output_filepath}'")
