#! /usr/bin/env bash

#BSUB -J bwa[1-8]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/ref/antisense-Ecoli.fa"
in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams/antisense"

set -o nounset -o pipefail -o errexit -x

samples=(
Ecoli002_20231208_1312_P2S-00519-A_PAS44628_9e5214a7.rna002_70bps_fast@v3
Ecoli002_20231208_1312_P2S-00519-A_PAS44628_9e5214a7.rna002_70bps_hac@v3
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_fast@v3.0.1
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_fast@v5.0.0
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_hac@v3.0.1
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_hac@v5.0.0
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_sup@v3.0.1
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_sup@v5.0.0
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam
