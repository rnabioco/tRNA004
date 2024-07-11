import pysam
import pandas as pd
import matplotlib.pyplot as plt

def extract_ml_values_with_positions(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")

    ml_values = []

    for read in bam:
        if read.has_tag('MM') and read.has_tag('ML'):
            mm_tag = read.get_tag('MM')
            ml_tag = read.get_tag('ML')
            seq = read.query_sequence
            start_pos = read.reference_start

            # Debugging output
            print(f"\nRead start position: {start_pos}")
            print(f"Sequence: {seq}")
            print(f"MM tag: {mm_tag}")
            print(f"ML tag: {ml_tag}")

            mm_parts = mm_tag.split(';')
            ml_parts = list(ml_tag)

            ml_index = 0  # Index to keep track of ML values

            for mm_part in mm_parts:
                if not mm_part:
                    continue

                mod_info = mm_part.split(',')
                base_and_mod = mod_info[0]

                if len(base_and_mod) < 3:
                    # Debugging output
                    print(f"Skipping invalid base_and_mod: '{base_and_mod}'")
                    continue

                base = base_and_mod[0]
                mod_type = base_and_mod[2:]

                if mod_type == '17802':
                    deltas = list(map(int, mod_info[1:]))
                    print(f"Deltas: {deltas}")

                    seq_position = 0
                    base_count = 0

                    for delta in deltas:
                        # Move to the next position by counting the bases
                        base_count = 0
                        while base_count <= delta and seq_position < len(seq):
                            if seq[seq_position] == base:
                                base_count += 1
                            seq_position += 1

                        current_pos = start_pos + seq_position - 1
                        try:
                            if ml_index < len(ml_parts):
                                ml_value = ml_parts[ml_index]
                                ml_values.append({'Position': current_pos, 'ML_value': ml_value})
                                # Debugging output
                                print(f"Appended: Position: {current_pos}, ML Value: {ml_value}")
                                ml_index += 1
                            else:
                                print(f"ml_index ({ml_index}) out of range for ml_parts length {len(ml_parts)}")
                        except IndexError as e:
                            print(f"IndexError: {e}, Length of ml_values: {len(ml_values)}, ML parts length: {len(ml_parts)}")
                            break

    # Debugging output
    print(f"\nML Values Extracted: {ml_values}")
    df = pd.DataFrame(ml_values)
    print(f"\nExtracted DataFrame:\n{df.head()}")  # Debugging print statement
    return df

# Path to your subsampled BAM file
subsampled_bam_file = "/Users/laurakwhite/Library/CloudStorage/Dropbox/Laura - UC Denver/2023/tRNAworkshop/alignments/004supv5/subsampled.yeast.10k.bam"

# Extract ML values within the sequence
ml_df = extract_ml_values_with_positions(subsampled_bam_file)

# Check if 'Position' column exists and if DataFrame is not empty
if not ml_df.empty and 'Position' in ml_df.columns:
    # Group by position and calculate the median ML value
    ml_median_df = ml_df.groupby('Position')['ML_value'].median().reset_index()

    # Plot the median ML values from 5' to 3'
    plt.figure(figsize=(10, 5))
    plt.plot(ml_median_df['Position'], ml_median_df['ML_value'], color='blue')
    plt.title("Median ML Values from 5' to 3'")
    plt.xlabel('Position')
    plt.ylabel('Median ML Value')
    plt.show()
else:
    print("No valid positions extracted or 'Position' column missing.")
