#! /usr/bin/env bash

#BSUB -J bwa[1-12]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 18

bwaidx="$HOME/tRNAworkshop/ref/sacCer-mito-and-nuclear-tRNAs.fa"
in="$HOME/tRNAworkshop/rebasecalled/fastqs"
dest="$HOME/tRNAworkshop/rebasecalled/alignedbams/mito_analysis"

set -o nounset -o pipefail -o errexit -x

samples=(
WTyeast004_20231111_1104_P2S-00519-B_PAQ47538_fa3726ec.rna004_130bps_sup@v3.0.1
WTyeast004_20231111_1104_P2S-00519-B_PAQ47538_fa3726ec.rna004_130bps_sup@v5.0.0
Petityeast004_20240105_1222_P2S-00519-A_PAQ47538_49891fce.rna004_130bps_sup@v3.0.1
Petityeast004_20240105_1222_P2S-00519-A_PAQ47538_49891fce.rna004_130bps_sup@v5.0.0
Pus1yeast004_20240105_1330_P2S-00519-A_PAS92267_56de71b8.rna004_130bps_sup@v3.0.1
Pus1yeast004_20240105_1330_P2S-00519-A_PAS92267_56de71b8.rna004_130bps_sup@v5.0.0
Pus2yeast004_20240105_1330_P2S-00519-B_PAS98920_19ba906c.rna004_130bps_sup@v3.0.1
Pus2yeast004_20240105_1330_P2S-00519-B_PAS98920_19ba906c.rna004_130bps_sup@v5.0.0
Pus4yeast004_20231208_1152_MN42516_FAX71838_e2d70aac.rna004_130bps_sup@v3.0.1
Pus4yeast004_20231208_1601_MN31004_FAX30455_c72fd68f.rna002_130bps_sup@v5.0.0
Grandeyeast004_20240105_1222_P2S-00519-B_PAS98845_231d7d7e.rna004_130bps_sup@v3.0.1
Grandeyeast004_20240105_1222_P2S-00519-B_PAS98845_231d7d7e.rna004_130bps_sup@v5.0.0
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u}.fastq \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam
