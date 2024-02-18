#! /usr/bin/env bash

#BSUB -J bwa[1-5]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/refs_with_adapters/dm6-mature-tRNAs.fa"
in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams"

set -o nounset -o pipefail -o errexit -x

samples=(
Drosophila002_20231207_1500_MN35252_FAX30455_48242de2.rna002_70bps_fast@v3
Drosophila002_20231207_1500_MN35252_FAX30455_48242de2.rna002_70bps_hac@v3
Drosophila004_20231208_1307_MN31004_FAX73799_72c200e4.rna004_130bps_fast@v3.0.1
Drosophila004_20231208_1307_MN31004_FAX73799_72c200e4.rna004_130bps_hac@v3.0.1
Drosophila004_20231208_1307_MN31004_FAX73799_72c200e4.rna004_130bps_sup@v3.0.1
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam

