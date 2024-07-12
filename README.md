# tRNA004
Nanopore direct tRNA sequencing method refinement and benchmarking of RNA modification detection between old (RNA002) and new (RNA004) chemistry.

Accompanying manuscript describes:
* our method improvements to nanopore tRNA sequencing protocols (use of beads, optimizing molar ratios of ligations to reduce free adapter)
* the yield and basecalling accuracy improvements that result from adapting the protocol to RNA004 chemistry
* ability to capture mitochondrial tRNA 
* leveraging tRNA as a way to systematically benchmark differences in basecalling error signals produced by old vs. new nanopore dRNAseq chemistry
* comparison of select tRNA modification signals across ecoli, budding yeast, drosophila, tetrahymena, T4 phage, human tRNAs

The directory `modplots` contains additional plots for all 43 RNA modifications evaluated in the manuscript, at all sequence contexts >29X coverage in RNA004 supv5 basecalled data.
