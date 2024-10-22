## Introduction

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12802127.svg)](https://doi.org/10.5281/zenodo.12802127)

This project focuses on nanopore direct tRNA sequencing method refinement and benchmarking of RNA modification detection between Oxford Nanopore's previous direct RNA sequencing chemistry, RNA002, and the RNA004 chemistry released in November 2023. The majority of the data was collected during a December 2023 tRNA sequencing workshop sponsored by the Hesselberth lab at the University of Colorado, in which participants prepared matched sequencing libraries using both chemistries from tRNA isolated from six different species.

## Comparative analysis of 43 distinct RNA modifications by nanopore tRNA sequencing
The accompanying manuscript describes how we have leveraged the many well characterized modifications on tRNA (cataloged in [MODOMICS](https://genesilico.pl/modomics/)) to evaluate the signals produced by nanopore basecallers at more than 40 unique types of RNA modifications across both the previous and current nanopore direct RNA sequencing chemistry. In the process, we also quantified the yield and basecalling accuracy improvements that result from switching to RNA004 chemistry for tRNA sequencing, and take advantage of these these higher yields for mitochondrial and bacteriophage tRNA sequence analysis, including RNA modification detection on low abundance tRNAs.

## Data analysis
This repository contains an R project and associated R markdown documents to generate the plots throughout the manuscript (in `rmd`), as well as custom scripts used for alignment, read filtering, quantification of abundance, calculation of basecalling errors and annotation of tRNA references with modification location in `src`. The subdirectory `bcerror` contains tsvs with outputs of some of these for modification induced error analysis, which are directly called by the R markdown documents to generate plots. The directory `modplots` contains additional plots for all 43 RNA modifications evaluated in the manuscript, at all sequence contexts >29X coverage in RNA004 supv5 basecalled data. These plots are not included in the manuscript, but are generated using the function `plot_modification_data` in `rmd/fig2_RNA_mods.Rmd`.

## Data access
Raw sequencing data (pod5 and fastq files) associated with this manuscript have been uploaded to SRA and Zenodo for private review and will be made publically available upon upon final publication.

## tRNA sequencing logistics
Our manuscript also details a few methodological improvements to nanopore tRNA sequencing, including the use of tRNA-specific SPRI beads and optimized molar ratios for all ligations, as well as updated sequencing settings for runs using the RNA004 chemistry.

### Nanopore tRNA sequencing protocols
* [Condensed protocol](https://docs.google.com/document/d/1gpVMYq6TfkbFxCs395O-HhrZjDWMoeaVCyw4LQSlWLU/edit?usp=sharing) (Google doc)
* [Long form protocol](https://benchling.com/s/prt-hGDLQuwsw0qKaWrtc70u?m=slm-kowx2TQzGlx2LAxpZwXk) with dynamic tables for calculating reactions with varying inputs (Benchling, thanks to Kezia Dobson)

### Sequencing settings
As of July 2024, the default MinKNOW settings for RNA004 libraries recapitulate the solution to the read length problem identified in [Lucas et al. 2023](https://www.nature.com/articles/s41587-023-01743-6) and described in more detail for RNA002 libraries [here](https://lkwhite.github.io/tRNAseq/sequencing-settings/), enabling real time capture of short reads without the need for saving bulk files, running sequencing simulations or any editing of config files. Our matched libraries in this dataset were sequenced using the same millisecond-based threshold for defining reads vs. adapters for both the RNA004 and RNA002 chemistry.

NB: an upcoming "RNA short mode" announced by ONT in May 2024 has the potential to alter these settings; you can check them yourself in your MinKNOW directory under `/Applications/MinKNOW.app/Contents/Resources/conf/package/sequencing` on a Mac, or `/opt/ont/minknow/conf/package/` on a Linux machine by examine `sequencing_PRO004RA_RNA.toml` and/or `sequencing_MIN106_RNA.toml`, depending on which instrument you use.
