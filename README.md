# TARGET CAPTURE ANALYSES
## ABOUT
The following script aims to analyze exons obtained by target capture sequencing. Nudibranch mollusc species were sequenced using the bait set designed by XXX. We evaluated its effectiveness in gene recovery across broader phylogenetic clades by sequencing both taxa within the bait set's target scope and more divergent taxa. For this purpose, we partially followed the same methodology previously applied to nudibranchs by Layton et al. (2020) and XXXX. The software Hybpiper, designed for targeted sequence capture, was applyied for assembly and extract target genes and and a phylogenetic analysis was conducted. Further analyses for efficiency test were carried out in R using bait-to-target DNA distances and Hybpiper statistics.

Citation: Paz-Sedano S., Valdes A., Stout C.C., Feliciano K., Muro S., Wilson N., Layton K., Goodheart J.A. 202X. Assessing the limits of exon capture efficiency for phylogenomics, using nudibranch gastropods (Mollusca: Heterobranchia) as a case study. 

## TRIMMOMATIC
Trimmomatic v0.36 (Bolger, Lohse, & Usadel, 2014) was used to remove adapter sequences, exon capture reads with a quality score below 15 in a 4-bp sliding window, and reads shorter than 26 bp.

      #!/bin/bash
      for folder_name in */; do   
      # Remove the trailing slash
      folder_name=${folder_name%/};   
      # Ensure it's a directory 
      if [ -d "$folder_name" ]; then     
      cd "$folder_name";     

      trimmomatic PE -threads 4 -phred33 ${folder_name}_R1_001.fastq.gz ${folder_name}_R2_001.fastq.gz ${folder_name}_R1_paired.fq.gz ${folder_name}_R1_unpaired.fq.gz ${folder_name}_R2_paired.fq.gz ${folder_name}_R2_unpaired.fq.gz ILLUMINACLIP:../ALL_adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:26;
      cd ..;
      fi;
      done
|options|function|
|:------|:-----|
|-threads|indicates the number of threads to use, which improves performance on multi-core computers. If not specified, it will be chosen automatically.|
|-phred33 or -phred64|specifies the base quality encoding. If no quality encoding is specified, it will be determined automatically (since version 0.32). The prior default was -phred64.| 
|ILLUMINACLIP| includes fasta with adapters, seed mismatches, palindrome clip threshold, simple clip threshold|
|SLIDINGWINDOW| including windowSize and required quality| 
|windowSize|specifies the number of bases to average across|
|requiredQuality|specifies the average quality required.|
|MINLEN|This module removes reads that fall below the specified minimal length.|

## HYBPIPER
Clean reads were analyzed using HybPiper v2.2.0 (Johnson et al. 2016). The full tutorial can be found at: [hybpiper](https://github.com/mossmatters/HybPiper)

### hybpiper assemble

Assemble contigs and extract sequences. To organize everything like the tutorial, create a folder called hybpiper_exon with all the clean reads.

Create a namelist.txt with the names of all the samples

      #!/bin/bash
      while read name; 
      do 
          hybpiper assemble -t_dna LN03_Ktar_chr_concatExonsDeduped_DNA.fasta -r ${name}_R1.fq.gz ${name}_R2.fq.gz --prefix ${name} --bwa; #change the name of your fasta file depending on the bait set name
      done < namelist.txt

|options|function|
|:------|:-----|
|-t|Target sequences file, which can be -t_dna for nucleotides and -t_aa for amino acids|
|-r|cleaned read fileas as input|
|--prefix|Prefix for sample name|
|--BWA| BWA for assembly. It works with a NUCLEOTIDE target file. If you supply a nucleotide target file but omit the flag --bwa, HybPiper will translate your target file and map reads to targets using BLAST.|

The BLASTx version of the pipeline (default) will map the reads to amino-acid target sequences. Although it is slower than the BWA version, it may have higher specificity. Reads may not align to divergent nucleotide target sequences, which are required for the BWA version. 

The reads are distributed to separate directories, where they are assembled separately using SPAdes. The main output is a FASTA file of the (in frame) CDS portion of the sample for each target region, and a separate file with the translated protein sequence.

HybPiper will:
1.	map the reads to the target sequences,
2.	 sort the reads by gene, 
3.	assemble the reads for each gene separately,
4.	align the contigs to the target sequence, 
5.	extract a coding sequence from each gene.
   
Output from each of these phases is saved in a standardized directory hierarchy, making it easy for post-processing commands to summarize information across many samples. It creates one additional folder for each gene extracted for that specimen, where you can find fasta files of the outputs generated.

### hybpiper stats

The parent directory contains one or more Base directories corresponding to the output of hybiper assemble for each sample. The descriptions below assume that the command hybpiper stats have been run from the parent directory.

      #!/bin/bash
      hybpiper stats -t_dna LN03_Ktar_chr_concatExonsDeduped_DNA.fasta gene namelist.txt #change the name of your fasta file depending on the bait set name

output:

•	seq_lengths.tsv. A table in tab-separated-values format, containing the lengths of each recovered gene sequence for each sample, along with the mean sequence length for each gene within the target file. The name of this file can be changed using the parameter --seq_lengths_filename <filename>.

•	hybpiper_stats.tsv. A table in tab-separated-values format, containing statistics on the HybPiper run. The name of this file can be changed using the parameter --stats_filename <filename>.

Output files of hybpiper stats were used for further investigation on the target capture efficiency CAN I INCLUDE HERE A LINK TO R ANALYSES?


### hybpiper retrieve_sequences
This command fetches either:
1.	The sequences recovered from the same gene for all samples; generates an unaligned multi-FASTA file for each gene.
2.	The sequences recovered from each gene for a single sample; generates a single multi-FASTA file.

            #!/bin/bash
            hybpiper retrieve_sequences -t_dna LN03_Ktar_chr_concatExonsDeduped_DNA.fasta dna --sample_names namelist.txt #change the name of your fasta file depending on the bait set name



## REFERENCES
(Bolger, Lohse, & Usadel, 2014) 
(Johnson et al. 2016)
Layton et al 2020
