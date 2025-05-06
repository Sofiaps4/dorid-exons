# PHYLOGENETIC ANALYSES

The retrieved sequences were aligned and trimmed, and phylogenetic analyses were carried out, following the tutorial [GeneTrees](https://hackmd.io/@mossmatters/SkXRVOdSc)

## MAFFT
   
The resulting gene files were aligned using MAFFT v7 (Katoh & Standley, 2013). To be able to infer a phylogeny, we first need to align the sequences. MAFFT takes unaligned raw sequences and creates multiple sequence alignments of nucleotide sequences.

      #!/bin/bash
      for i in retrieve/*.FNA; do   
      gene=$(basename $i | cut -f 1 -d .)

      mafft --preservecase --maxiterate 1000 --localpair retrieve/$gene.*.FNA > MAFFT/$gene.mafft.fasta;
      done


## Trimal
After aligning sequences, there will be spaces in between base pairs that need to be removed before building the tree. Trimal will remove any illegitimate or poorly aligned sections of the sequences.
First, make a new directory named TRIMAL, mkdir TRIMAL for your trimal output files to be directed to.

      #!/bin/bash
      for i in MAFFT/*.mafft.fasta; do   
      gene=$(basename $i | cut -f 1 -d .)
      trimal -in MAFFT/$gene.mafft.fasta -out trimal/$gene.trimal.fasta;
      done 
      
|options|function|
|:------|:-----|
|-in|Input fasta file|
-out|Output trimmed fasta file|

### Occupancy matrix
The resulted aligned and trimmed files are fasta files for each gene, inside which there are the sequences for each sample that has that gene. 
We wanted to run phylogenetic analyses in different matrices with different occupancies, including in the matrix only those genes for which at least 25%, 50% and 75% of the taxa contained that gene. 
Gene alignments were copied in different folders, Matrix_all, Matrix_25, Matrix_50 and Matrix_75. Then, the genes fasta files with a lower percentage of species of interest in the occupation were eliminated

      for i in *.fasta; do 
      gene=$(basename "$i") 
      count=$(grep -o ">" "$i" | wc -l)
       echo "$gene: $count" 
      if [ "$count" -lt #number-occupancy ]; then #change the number depending on the % of occupancy you are interested in
      rm "$i" 
      fi
      done

As a result, only the genes of interest remained, which were concatenated using IQtree.


## IQTREE
Individual gene alignments were concatenated for model testing and maximum Likelihood analyses in IQtree v2.3.5 (Nguyen et al. 2015) were performed with 1,500 ultrafast bootstrap (UFBoot) replicates, using MoldelFinder (Nguyen et al. 2015; Kalyaanamoorthy et al. 2017) to find the best-fitting model for each concatenated matrix and matrix subset. Nodes with bootstrap support less than 50% were collapsed. 

IQ tree command will take your input of multiple sequence alignment and will reconstruct a phylogeny that is best explained by your input data. From version 2 you can input a directory of alignment files. IQ-TREE 2 will load and concatenate all alignments within the directory, eliminating the need for users to manually perform this step.

      iqtree2 -p FOLDER_NAME --out-aln OUTFILE_NAME

This will produce:

•	OUTFILE_NAME: the concatenated alignment (default is phylip format)

•	OUTFILE_NAME.nex: the partition file in nexus format

•	OUTFILE_NAME.partitions: the partition file in RAxML format

Optionally, you can add --out-format FASTA|NEXUS option to specify concatenated alignment format.
Once you have the concatenated nexus file:

      #!/bin/bash
      iqtree -s alignment.nex -m MFP+MERGE -bb 1500 -czb -nt AUTO
      
|options|function|
|:------|:-----|
|-s|imput alignment file|
|-m MFP+MERGE| Tells IQ-TREE to perform ModelFinder to find the best-fitting model and merge similar partitions based on a statistical test.|
|-bb 1000|Indicates that 1000 ultrafast bootstrap replicates should be performed.|
|-czb|Collapses nodes with bootstrap support less than 50%.|

