#! /usr/bin/env bash

#BSUB -J bwa[1-11]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/refs_with_adapters/hg38-mature-tRNAs.fa"
in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams"

set -o nounset -o pipefail -o errexit -x

samples=(
Trub1doxH295R004_20231207_1527_MN42516_FAX73799_8b069090.rna004_130bps_fast@v3.0.1
Trub1doxH295R004_20231207_1527_MN42516_FAX73799_8b069090.rna004_130bps_hac@v3.0.1
Trub1doxH295R004_20231207_1527_MN42516_FAX73799_8b069090.rna004_130bps_sup@v3.0.1
Trub1untrH295R004_20231207_1522_MN31004_FAX71838_4571d21c.rna004_130bps_fast@v3.0.1
Trub1untrH295R004_20231207_1522_MN31004_FAX71838_4571d21c.rna004_130bps_hac@v3.0.1
Trub1untrH295R004_20231207_1522_MN31004_FAX71838_4571d21c.rna004_130bps_sup@v3.0.1
HumanFibroblast002_20231207_1455_P2S-00519-B_PAS29437_b7f6c4bd.rna002_70bps_fast@v3
HumanFibroblast002_20231207_1455_P2S-00519-B_PAS29437_b7f6c4bd.rna002_70bps_hac@v3
HumanFibroblast004_20231208_1351_MN42516_FAX71838_b8f7b990.rna004_130bps_fast@v3.0.1
HumanFibroblast004_20231208_1351_MN42516_FAX71838_b8f7b990.rna004_130bps_hac@v3.0.1
HumanFibroblast004_20231208_1351_MN42516_FAX71838_b8f7b990.rna004_130bps_sup@v3.0.1
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam
