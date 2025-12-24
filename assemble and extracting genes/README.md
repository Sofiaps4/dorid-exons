 ## TRIMMOMATIC
Trimmomatic v0.36 was used to remove adapter sequences, trim exon capture reads with a quality score below 15 in a 4-bp sliding window, and reads shorter than 26 bp.

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
|Options|Description|
|:------|:-----|
|-threads|Number of threads to use, which improves performance on multi-core computers. If not specified, it will be chosen automatically.|
|-phred33|Specifies the base quality encoding. If no quality encoding is specified, it will be determined automatically. The prior default was -phred64.| 
|ILLUMINACLIP| Specifies adapters file, seed mismatches, palindrome clip threshold, simple clip threshold|
|SLIDINGWINDOW|Performs quality trimming using a sliding window (window size: 4, quality threshold: 15)| 
|MINLEN|Discards reads that fall below the specified minimal length.|

## HYBPIPER
Clean reads were analyzed using HybPiper v2.2.0. The full tutorial can be found at: [hybpiper](https://github.com/mossmatters/HybPiper)

### Hybpiper assemble

Assemble contigs and extract sequences from exon-capture reads. To organize everything like the tutorial, 
- Create a folder called hybpiper_exon with all the clean reads.
- Create a namelist.txt with the names of all the samples.

       #!/bin/bash
       while read name; 
       do 
       hybpiper assemble -t_dna BaitSetName.fasta -r ${name}_R1.fq.gz ${name}_R2.fq.gz --prefix ${name} --bwa; 
       done < namelist.txt
    

|options|function|
|:------|:-----|
|-t|Target sequences file, which can be -t_dna for nucleotides and -t_aa for amino acids|
|-r|cleaned read as input|
|--prefix|Prefix for each sample's output files|
|--BWA| BWA for assembly. It works with a NUCLEOTIDE target file. If you supply a nucleotide target file but omit the flag --bwa, HybPiper will translate your target file and map reads to targets using BLAST.|

The BLASTx version of the pipeline (default) will map the reads to amino-acid target sequences. Although it is slower than the BWA version, it may have higher specificity. Reads may not align to divergent nucleotide target sequences, which are required for the BWA version. 

The reads are distributed to separate directories, where they are assembled separately using SPAdes. The main output is a FASTA file of the (in frame) CDS portion of the sample for each target region, and a separate file with the translated protein sequence.

HybPiper will:
1.	map the reads to the target sequences,
2.	 sort the reads by gene, 
3.	assemble the reads for each gene separately,
4.	align the contigs to the target sequence, 
5.	extract a coding sequence from each gene.

Outputs are saved in a standardized directory hierarchy, making it easy for post-processing commands to summarize information across many samples. It creates one additional folder for each gene extracted for that specimen, where you can find fasta files of the outputs generated.

### Hybpiper stats

The parent directory contains one or more Base directories corresponding to the output of hybiper assemble for each sample. The descriptions below assume that the command hybpiper stats have been run from the parent directory.

      #!/bin/bash
      hybpiper stats -t_dna BaitSetName.fasta gene namelist.txt #change the name of your fasta file depending on the bait set name

Output files:

•	seq_lengths.tsv. A table in tab-separated-values format, containing the lengths of each recovered gene sequence for each sample, along with the mean sequence length for each gene within the target file. The name of this file can be changed using the parameter --seq_lengths_filename <filename>.

•	hybpiper_stats.tsv. A table in tab-separated-values format, containing statistics on the HybPiper run. The name of this file can be changed using the parameter --stats_filename <filename>.

Output files of hybpiper stats were used for further investigation on the target capture efficiency [efficiency_analyses](https://github.com/Sofiaps4/dorid-exons/tree/main/efficiency_analyses)


### Hybpiper retrieve_sequences
This command fetches either:
1.	The sequences recovered from the same gene for all samples; generates an unaligned multi-FASTA file for each gene.
2.	The sequences recovered from each gene for a single sample; generates a single multi-FASTA file.

            #!/bin/bash
            hybpiper retrieve_sequences -t_dna BaitSetName.fasta dna --sample_names namelist.txt #change the name of your fasta file depending on the bait set name

The retrieved sequences were used as input files for further [phylogenetic analyses](https://github.com/Sofiaps4/dorid-exons/tree/main/phylogeny) as well as calculate bait-to-target DNA distances for [efficiency_analyses](https://github.com/Sofiaps4/dorid-exons/tree/main/efficiency_analyses)



### FILTERING OUTLIERS

To filter possible contaminations, paralogs and general outlier genes obtained during the assembly and retrieve steps, we run [phylTER](https://github.com/damiendevienne/phylter) tool in R. phylTER use a collection of gene trees to look for outliers, selecting specific gene in a specific taxon. It use a collection of gene trees as an input data. Therefore, we first aligned retrieved sequences using MAFFT v7 (Katoh & Standley, 2013) and
run IQtree v2.3.5 (Nguyen et al. 2015) for gene trees, using MoldelFinder (Kalyaanamoorthy et al. 2017) to find the best-fitting model for each gene.  


## MAFFT

         #!/bin/bash
         for gene in MAFFT/*.fasta.gz; do
          base=$(basename "$gene" .FNA.gz)
          gunzip -c "$gene" > tmp_fna/$base.FNA
         iqtree2 \
          -s tmp_fna/$base.FNA \
          -m MFP \
          -nt AUTO \
          -pre gene-trees/$base/$base
        done

## IQtree

         #!/bin/bash
         for gene in MAFFT/*.fasta.gz; do
         # Clean base name
           base=$(basename "$gene" .fasta.gz)
         # Decompress to temporary folder
           gunzip -c "$gene" > tmp_fna/$base.fasta
         # Run IQ-TREE
            mkdir -p gene-trees/$base  # ensure folder exists
           iqtree2 \
           -s tmp_fna/$base.fasta \
           -m MFP \
           -nt AUTO \
           -pre gene-trees/$base/$base
      done

The resulting gene trees were used as input for removing outliers in [phylTER](https://github.com/Sofiaps4/dorid-exons/blob/main/assemble%20and%20extracting%20genes/phylTER.R), run in R. Resulted outliers were removed from the initial retrieved data using the following code [filger_outliers](https://github.com/Sofiaps4/dorid-exons/blob/main/assemble%20and%20extracting%20genes/filter_outliers.py). These filtered genes were used for further [phylogenetic](https://github.com/Sofiaps4/dorid-exons/tree/main/phylogeny) and [capture efficiency](https://github.com/Sofiaps4/dorid-exons/tree/main/efficiency_analyses) analyses. 
