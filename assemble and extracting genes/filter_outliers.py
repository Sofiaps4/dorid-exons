#!/usr/bin/env python3
import gzip
from Bio import SeqIO
import os
import shutil

# ----------------------------
# Input / Output paths
# ----------------------------
outlier_file = "phylter_outliers_result.txt"  # txt with gene + taxa
fasta_folder = "retrieve"                     # folder with all .FNA.gz files
output_folder = "retrieve_filtered"           # folder to store filtered FASTAs
os.makedirs(output_folder, exist_ok=True)

# ----------------------------
# Read outlier table
# ----------------------------
# Dictionary: gene -> set of taxa to remove
outliers = {}
with open(outlier_file, "r") as f:
    next(f)  # skip header
    for line in f:
        # Strip whitespace AND remove Windows \r
        line = line.strip().replace('\r', '')
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        gene = parts[0].strip()
        taxon = parts[1].strip()
        # Store lower-case version for case-insensitive matching
        outliers.setdefault(gene, set()).add(taxon.lower())

# ----------------------------
# Process all FASTA files
# ----------------------------
for fasta_file in os.listdir(fasta_folder):
    if not fasta_file.endswith(".FNA.gz"):
        continue

    gene = fasta_file.replace(".FNA.gz", "")
    input_path = os.path.join(fasta_folder, fasta_file)
    output_path = os.path.join(output_folder, fasta_file)

    if gene in outliers:
        taxa_to_remove = outliers[gene]

        with gzip.open(input_path, "rt") as infile:
            seqs = list(SeqIO.parse(infile, "fasta"))

            # Filter sequences: check first word in header, lower-case
            seqs_filtered = [
                seq for seq in seqs
                if seq.id.split()[0].strip().lower() not in taxa_to_remove
            ]

            if len(seqs_filtered) == 0:
                # If all sequences are outliers, delete the gene (do not write file)
                print(f"All sequences removed for gene {gene}. Skipping this gene.")
                continue

            # Otherwise, write filtered sequences
            with gzip.open(output_path, "wt") as outfile:
                SeqIO.write(seqs_filtered, outfile, "fasta")

        print(f"Processed {gene}: kept {len(seqs_filtered)}/{len(seqs)} sequences")
    
    else:
        # Copy unmodified gene if not in txt
        shutil.copy2(input_path, output_path)
        print(f"Copied {gene} without changes")

