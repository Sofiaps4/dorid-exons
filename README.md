# SCRIPTS FOR TARGET CAPTURE ANALYSES

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

-phred33 or -phred64 specifies the base quality encoding. If no quality encoding is specified, it will be determined automatically (since version 0.32). The prior default was -phred64. 

-threads indicates the number of threads to use, which improves performance on multi-core computers. If not specified, it will be chosen automatically.

ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>

LEADING:<quality> Remove low quality bases from the beginning. As long as a base has a value below this threshold the base is removed and the next base will be investigated

quality: Specifies the minimum quality required to keep a base.

TRAILING:<quality> Remove low quality bases from the end

quality: Specifies the minimum quality required to keep a base.

SLIDINGWINDOW:<windowSize>:<requiredQuality> 

windowSize: specifies the number of bases to average across 

requiredQuality: specifies the average quality required.

MINLEN:<length> This module removes reads that fall below the specified minimal length.

## HYBPIPER
HybPiper v1.3.1 (Johnson et al., 2016) was employed to assemble the cleaned reads into contigs of the targeted regions of the genes. Briefly, the reads were mapped to a reference file of con-catenated bait sequences using Burrows-Wheeler Aligner (BWA) (Li& Durbin, 2009).
HybPiper produced an unaligned fasta file for each gene, containing a DNA sequence for each sample, and a series of summary statistics. For example, Hybpiper uses BWA to map the reads to contigs to present a value for percent reads on target.Genes that did not enrich or enriched poorly (genes whose contigswere <50% of the reference) were removed (n = 125)
FULL TUTORIAL: [hybpiper](https://github.com/mossmatters/HybPiper)

### hybpiper assemble

Assemble contigs and extract sequences. 
To organize everything like the tutorial, create a folder called hybpiper_exon with all the clean reads
Create a namelist.txt with the names of all the samples

      while read name; 
      do 
          hybpiper assemble -t_dna LN03_Ktar_chr_concatExonsDeduped_DNA.fasta -r ${name}_R1.fq.gz ${name}_R2.fq.gz --prefix ${name} --bwa; #change the name of your fasta file depending on the bait set name
      done < namelist.txt

--BWA, you need a NUCLEOTIDE target file. If you supply a nucleotide target file but omit the flag --bwa, HybPiper will translate your target file and map reads to targets using BLAST.

The reads are distributed to separate directories, where they are assembled separately using SPAdes. The main output is a FASTA file of the (in frame) CDS portion of the sample for each target region, and a separate file with the translated protein sequence.

HybPiper will:
1.	map the reads to the target sequences,
2.	 sort the reads by gene, 
3.	assemble the reads for each gene separately,
4.	align the contigs to the target sequence, 
5.	extract a coding sequence from each gene.
   
Output from each of these phases is saved in a standardized directory hierarchy, making it easy for post-processing commands to summarize information across many samples.

It creates one additional folder for each gene extracted for that specimen, where you can find fasta files of the outputs generated.

### hybpiper stats

The parent directory contains one or more Base directories corresponding to the output of hybiper assemble for each sample. The descriptions below assume that the command hybpiper stats have been run from the parent directory.

      hybpiper stats -t_dna LN03_Ktar_chr_concatExonsDeduped_DNA.fasta gene namelist.txt #change the name of your fasta file depending on the bait set name



