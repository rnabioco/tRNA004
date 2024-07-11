#! /usr/bin/env bash

#BSUB -J bwa[1-16]
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -q rna
#BSUB -n 6

bwaidx="$HOME/tRNAworkshop/ref/sacCer-mito-and-nuclear-tRNAs.fa"
in="$HOME/tRNAworkshop/alignment_evals/syntheticreads"
dest="$HOME/tRNAworkshop/alignment_evals/bams"

set -o nounset -o pipefail -o errexit -x

samples=(
Q5synthetic.fasta
Q6synthetic.fasta
Q7synthetic.fasta
Q8synthetic.fasta
Q9synthetic.fasta
Q10synthetic.fasta
Q11synthetic.fasta
Q12synthetic.fasta
Q13synthetic.fasta
Q14synthetic.fasta
Q15synthetic.fasta
Q16synthetic.fasta
Q17synthetic.fasta
Q18synthetic.fasta
Q19synthetic.fasta
Q20synthetic.fasta
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

### align to tRNA reference - adapter anchored reference alignment
bwa mem -W 13 -k 6 -T 20 -x ont2d $bwaidx $in/${u} \
| samtools view -F4 -hu - \
| samtools sort -o ${dest}/${u}.bwa.bam

samtools index ${dest}/${u}.bwa.bam
