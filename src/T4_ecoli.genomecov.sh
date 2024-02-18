#! /usr/bin/env bash

#BSUB -J genomecov[1-2]
#BSUB -o logs/genomecov.%J.out
#BSUB -e logs/genomecov.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 1

in="$HOME/tRNAworkshop/rebasecalled/alignedbams/full_length"
dest="$HOME/tRNAworkshop/rebasecalled/bedgraphs/full_length"

set -o nounset -o pipefail -o errexit -x

samples=(
Ecoli004_20240105_1224_MN31004_FAX73799_1fe84761.rna004_130bps_sup@v3.0.1.T4ref
T4infectedecoli_20240202_1332_P2S-01618-A_PAU05281_fd3805fc.rna004_130bps_sup@v3.0.1.T4ref
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

# normalize by CPM - counts per million reads per lib
# determined from line counting the length of deduplicated bam files, dividing by 10^6
# genomecov scale factor is a multiplier, so inverse is calculated in tmpscale

tmpscale=$(samtools view -c ${in}/${u}.bwa_full_length.bam | awk '{print 1000000/$1}')

bedtools genomecov -ibam ${in}/${u}.bwa_full_length.bam -scale ${tmpscale} -bg > ${dest}/${u}.bg
