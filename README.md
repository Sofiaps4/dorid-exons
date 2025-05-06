# TARGET CAPTURE ANALYSES
## ABOUT
This repository contains scripts and instructions for analyzing exons obtained through target capture sequencing. Nudibranch mollusc species were sequenced using the bait set designed by XXX for the suborder Doridina. To evaluate the efficiency of this bait set across broader phylogenetic clades, we sequenced both Doridina taxa (within the bait set’s target scope) and more divergent taxa from the suborder Cladobranchia.
We partially followed the methodology previously applied to nudibranchs by Layton et al. (2020) and XXXX. The software Hybpiper, designed for targeted sequence capture, was applyied for assembly and extract target genes, and a phylogenetic analysis was conducted. Further analyses for efficiency test were carried out in R, studing the correlation between bait-to-target DNA distances and Hybpiper recovery statistics.

Citation: Paz-Sedano S., Valdes A., Stout C.C., Feliciano K., Muro S., Wilson N., Layton K., Goodheart J.A. 202X. Assessing the limits of exon capture efficiency for phylogenomics, using nudibranch gastropods (Mollusca: Heterobranchia) as a case study. 

## REQUIREMENTS
This workflow requires the following software: 

- [Trimmomatic](https://github.com/timflutre/trimmomatic) – for trimming raw reads

- [HybPiper](https://github.com/mossmatters/HybPiper) – for read assembly and target gene recovery

- [trimAl](https://vicfero.github.io/trimal/) – for trimming alignments

- [IQ-TREE](http://www.iqtree.org) – for phylogenetic inference

- [R](https://rstudio-education.github.io/hopr/starting.html) – for statistical analyses

We recommend installing and managing software using Conda. 

      conda install -c bioconda trimmomatic
      conda install -c bioconda hybpiper
      conda install -c bioconda trimal
      conda install -c bioconda iqtree


Installing HybPiper using conda with a new environment will install HybPiper, all required Python packages, and all required external programs. If you have conda installed and the channels bioconda and conda-forge have already been added, this can be done using the commands:

# Workflow overview

- [Assemble and extracting genes](https://github.com/Sofiaps4/dorid-exons/tree/main/assemble%20and%20extracting%20genes): for trimming raw sequences using Trimmomatic, assemble clean reads and extracting genes using Hybpiper.

- [Phylogeny](https://github.com/Sofiaps4/dorid-exons/tree/main/phylogeny): for phylogenetics analyses using extracted genes as input, including aligning with MAFFT, trimming alignments with TrimAl and infer Maximum Likelihood with IQtree.

- [Efficiency_analyses](https://github.com/Sofiaps4/dorid-exons/tree/main/efficiency_analyses): for statistics analyses of capture efficiency using R.
