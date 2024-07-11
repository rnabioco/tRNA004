#!/usr/bin/env python3

import pysam
import os
import csv
import sys

def calculate_precision_recall(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    total_reads = bam.mapped + bam.unmapped  # total number of reads in the BAM file
    tp = 0  # true positives (forward alignments)
    fp = 0  # false positives (reverse alignments)
    fn = 0  # false negatives (unmapped reads)

    for read in bam:
        if read.is_unmapped:
            # Unmapped reads are counted as false negatives
            fn += 1
        else:
            if read.flag & 16:  # read is mapped to the reverse strand
                # Reads mapped to the reverse strand are false positives
                fp += 1
            else:
                # Reads mapped to the forward strand are true positives
                tp += 1

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0

    return [total_reads, tp, fp, fn, precision, recall]

def write_aggregate_csv(output_csv, results):
    with open(output_csv, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['File', 'Total Reads', 'TP', 'FP', 'FN', 'Precision', 'Recall'])
        for result in results:
            writer.writerow(result)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate precision and recall for BAM files.')
    parser.add_argument('input_dir', type=str, help='Directory containing input BAM files.')
    parser.add_argument('output_csv', type=str, help='Output CSV file to save the aggregated results.')

    args = parser.parse_args()

    input_dir = args.input_dir
    output_csv = args.output_csv

    results = []

    # Process each BAM file in the input directory
    for bam_file in os.listdir(input_dir):
        if bam_file.endswith('.bam'):
            bam_path = os.path.join(input_dir, bam_file)
            metrics = calculate_precision_recall(bam_path)
            base_name = os.path.basename(bam_file).replace('.bam', '')
            results.append([base_name] + metrics)
    
    write_aggregate_csv(output_csv, results)

