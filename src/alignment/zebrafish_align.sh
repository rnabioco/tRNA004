#! /usr/bin/env bash

#BSUB -J bwa[1-3]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/ref/danRer11-mature-tRNAs.fa"
in="$HOME/tRNAworkshop/rebasecalled/supv5_fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams"

set -o nounset -o pipefail -o errexit -x

samples=(
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_fast@v5.0.0
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_hac@v5.0.0
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_sup@v5.0.0
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -C -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam
