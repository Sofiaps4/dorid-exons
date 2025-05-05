# EFFICIENCY ANALYSES

To observe if the capture efficiency differs among taxa, descriptive statistics were performed, [descriptive_statistic.R](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/descriptive_statistic.R), exploring whether sensitivity (percentage of genes) and specificity (percentage of reads that map to targets) vary depending on the taxa. 
The number of genes for each taxa, as well as percentage of reads that map to targets (expressed as PctOnTarget) are obtained during hybpiper stats analyses, file hybpiper_stats.tsv

To investigated the hypotheses of whether capture sensitivity or capture specificity are correlated to the bait-to-target DNA distance we calculate the genetic distances.
The bait set includes sequences for the same gene from different species. However, the specific sequence in the bait set used for recovery of target genes can be seen in the Exonerate output during the extraction and assembly of exons in Hybpiper. 
To first obtain the bait-to-target DNA distance, we observed which of the baits comprising the bait set was used to extract each gene for each species. A list was created that associated each bait with the target [list_target_genes.sh](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/list_target_genes.sh)

The output list species_gene_headers.txt are the input of the following [calculate_distances.sh](https://github.com/Sofiaps4/dorid-exons/blob/main/efficiency_analyses/calculate_distances.sh)
This script generates a temporal alignment using MAFFT of each target gene extracted for each species, resulting from Hybpiper Retrieve, with the gene from the baitset used for its retrieve. 
Uncorrected pairwise distances were calculated, excluding gaps and missing data, and the average of the distances obtained was obtained, establishing it as the genetic distance of the species with the bait set. 

The hypotheses of whether capture sensitivity or capture specificity are correlated to the bait-to-target DNA distance were investigated using R, using the non-parametric Kendall rank correlation coefficient test.
