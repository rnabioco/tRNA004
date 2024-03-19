#! /usr/bin/env bash

#BSUB -J bwa[1-2]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/refs_with_adapters/T5_infected_ecoli.fa"
in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams/allreads"

set -o nounset -o pipefail -o errexit -x

samples=(
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_sup@v3.0.1
T5infectedecoli.20240202_1332_P2S-01618-B_PAS25444_31820489.rna004_130bps_sup@v3.0.1
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.T5ref.bwa.bam

samtools index ${dest}/${u}.T5ref.bwa.bam
