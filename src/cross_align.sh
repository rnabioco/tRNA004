#! /usr/bin/env bash

#BSUB -J bwa[1-6]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams/crossalignments"

set -o nounset -o pipefail -o errexit -x

samples=(
WTyeast004_20231111_1104_P2S-00519-B_PAQ47538_fa3726ec.rna004_130bps_sup@v3.0.1
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_sup@v3.0.1
Drosophila004_20231208_1307_MN31004_FAX73799_72c200e4.rna004_130bps_sup@v3.0.1
HumanFibroblast004_20231208_1351_MN42516_FAX71838_b8f7b990.rna004_130bps_sup@v3.0.1
Tetrahymena004_20231208_0920_MN42516_FAX71838_a244513d.rna004_130bps_sup@v3.0.1
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_sup@v3.0.1
)

refloc="/beevol/home/whitel/tRNAworkshop/ref"

references=(
ecoliK12MG1655-mature-tRNAs.fa
sacCer3-mature-tRNAs.fa
tthermophila-numbered-mature-tRNAs.fa
dm6-mature-tRNAs.fa
danRer11-mature-tRNAs.fa
hg38-mature-tRNAs.fa
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

# Loop over each reference
for ref in "${references[@]}"; do
    # Extract a simple name for the reference to use in naming the output files
    ref_name=$(basename "$ref" .fa)

    ### align to tRNA reference - adapter anchored reference alignment
    bwa mem -W 13 -k 6 -T 20 -x ont2d "$refloc/$ref" "$in/${u}.fastq" \
    | samtools view -F4 -hu - \
    | samtools sort -o "${dest}/${u}.${ref_name}.bwa.bam"

    samtools index "${dest}/${u}.${ref_name}.bwa.bam"
done
