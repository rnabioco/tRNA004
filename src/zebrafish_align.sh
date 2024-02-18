#! /usr/bin/env bash

#BSUB -J bwa[1-5]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/refs_with_adapters/danRer11-mature-tRNAs.fa"
in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams"

set -o nounset -o pipefail -o errexit -x

samples=(
Zebrafish002_20231207_1528_P2S-00519-A_PAS44221_7e583c4c.rna002_70bps_fast@v3
Zebrafish002_20231207_1528_P2S-00519-A_PAS44221_7e583c4c.rna002_70bps_hac@v3
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_fast@v3.0.1
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_hac@v3.0.1
Zebrafish004_20231208_0919_MN31004_FAX73799_f3dea9f8.rna004_130bps_sup@v3.0.1
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam
