# EFFICIENCY ANALYSES
To investigate the dataset's capture efficiency across divergent species, we observed variations in different parameters depending on the suborder and superfamily to which they belonged, including variations in:

- sensitivity (percentage of genes recovered)
- specificity (percentage of reads mapping to targets)
- bait-to-target DNA distances
- sequence length of genes relative to the bait set.
- number of reads

## Visualize stats
To observe if the capture efficiency differs among taxa, descriptive statistics were performed, [descriptive_statistic.R](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/descriptive_statistic.R)

The analyses focused on sensitivity and specificity (expressed as PctOnTarget), based on outputs from the hybpiper stats analysis, specifically the hybpiper_stats.tsv file.

Hybpiper stats output seq_lengths.tsv was used for plot the percentage of genes recovered at different proportions of sequence length relative to the bait set, classified by suborder [barplot_percents_lengths.py](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/barplot_percents_lengths.py)

    python barplot_percents_lengths.py

## Correlation analyses

To investigated the hypotheses of whether capture sensitivity or capture specificity are correlated to the bait-to-target DNA distance we calculate the genetic distances.

### Calculate bait-to-target DNA distance

The bait set includes sequences for the same gene from different species. However, the specific sequence in the bait set used for recovery of target genes can be seen in the Exonerate output during the extraction and assembly of exons in Hybpiper. 

First, we observed which of the baits comprising the bait set was used to extract each gene for each species. A list was created that associated each bait with the target [list_target_genes.sh](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/list_target_genes.sh)

The output list species_gene_headers.txt are the input for [calculate_distances.sh](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/calculate_distances.sh) script.

This script generates a temporal alignment using MAFFT, comparing each target gene retrieved for a species (via HybPiper Retrieve) with its corresponding gene from the bait set. Uncorrected pairwise distances were calculated, excluding gaps and missing data, are then calculated. The average of these distances is taken as the genetic distance between each species and the bait set.

### Correlation tests.

The hypotheses of whether capture sensitivity or capture specificity are correlated to the bait-to-target DNA distance were investigated using R: [correlation_analyses](https://github.com/Sofiaps4/dorid-exons/tree/main/efficiency_analyses)

The workflow performs the following steps:

- Computes the mean of the genetic distances obtained

- Assesses the normality of parameters using the Shapiro-Wilk test

- Evaluates the correlation between genetic distance and sensitivity using Kendall’s rank correlation test

- Applies piecewise (segmented) regression to model the relationship between distance and sensitivity

- Assesses the correlation between genetic distance and specificity using Kendall’s rank correlation test using the non-parametric Kendall rank correlation coefficient test.

This workflow also assess the correlation between number of initial reads included in the analyses and percentage of genes extracted, to investigate if the sensitivity could be related with input sequences more than genetic distances.

