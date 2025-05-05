#!/bin/bash

# ---------------------
# CONFIGURATION SECTION
# ---------------------
retrieve_dir="retrieve"                # Path to the 'retrieve' folder containing the FASTA files
bait_fasta="Nudibaitfile.fas"          # Path to the bait FASTA file
header_map="species_gene_headers.txt"  # Path to the header map file for species and gene mapping
long_csv="species_gene_distances_long.csv"   # Output CSV file (long format)
wide_csv="p_distance_matrix.csv"       # Output CSV file (wide format)

# ---------------------
# SETUP
# ---------------------
mkdir -p temp_alignments           # Create a temporary directory for alignments
echo "Species,Gene,Distance" > "$long_csv"  # Write header to the long CSV file

# ---------------------
# PROCESS EACH GENE
# ---------------------
for gene_file in "$retrieve_dir"/*.FNA.gz; do
  gene_name=$(basename "$gene_file" .FNA.gz)  # Extract gene name from file

  # Unwrap FASTA: Combine multiline sequences into single-line format
  zcat "$gene_file" | awk '
    /^>/ {if (seq) print seq; seq=""; print; next}
    {seq = seq $0}
    END {if (seq) print seq}
  ' > temp.fa

  # Iterate through each species in the file
  awk '
    NR%2==1 {header=$0}
    NR%2==0 {print header "\t" $0}
  ' temp.fa | while IFS=$'\t' read -r header sequence; do
    species_name=$(echo "$header" | cut -d' ' -f1 | sed 's/^>//')  # Extract species name

    # Find corresponding bait header from map
    target_name=$(awk -v s="$species_name" -v g="$gene_name" '$1 == s && $2 == g {print $3}' "$header_map")

    if [[ -z "$target_name" ]]; then
      echo "⚠️  No bait match for $species_name $gene_name, skipping"
      continue
    fi

    # Extract full bait sequence from Nudibaitfile.fas
    awk -v id="$target_name" '
      BEGIN { found=0 }
      $0 ~ "^>" {
        if (found) exit;
        found=($0 == id)
      }
      found { print }
    ' "$bait_fasta" > temp_ref.fa

    if [[ ! -s temp_ref.fa ]]; then
      echo "⚠️  No bait sequence found for $target_name, skipping"
      continue
    fi

    # Write query sequence to temp_query.fa
    echo ">$species_name" > temp_query.fa
    echo "$sequence" >> temp_query.fa
# Align both sequences using MAFFT
    cat temp_ref.fa temp_query.fa > temp_pair.fa
    mafft --quiet --auto temp_pair.fa > temp_aligned.fa

    # Check the number of sequences in the alignment
    seq_count=$(grep -c "^>" temp_aligned.fa)
    if [[ "$seq_count" -ne 2 ]]; then
      echo "⚠️  MAFFT alignment failed for $species_name $gene_name"
      continue
    fi

    # Calculate p-distance using Python (excludes gaps)
    distance=$(python3 - <<EOF
from Bio import AlignIO
import numpy as np

def calculate_p_distance(alignment):
    total_positions = 0
    differing_positions = 0

    # Iterate over aligned sequences, ignoring gaps
    for i in range(len(alignment[0].seq)):  # Iterate over the positions in the alignment
        residues = [record.seq[i].upper() for record in alignment]

        # Skip if there's a gap in either sequence
        if any(base in ['-', 'N'] for base in residues):
            continue

        # Count differing positions
        if residues[0] != residues[1]:
            differing_positions += 1

        total_positions += 1

    if total_positions > 0:
        p_distance = differing_positions / total_positions
    else:
        p_distance = 0  # If no positions to compare, set p-distance to 0

    return p_distance

# Read the alignment file
alignment = AlignIO.read("temp_aligned.fa", "fasta")

# Calculate the p-distance
p_distance = calculate_p_distance(alignment)

# Output the p-distance as a percentage difference
print(f"{p_distance:.5f}")
EOF
)

    # Log the species, gene, and distance
    echo "$species_name,$gene_name,$distance" >> "$long_csv"

  done
done

# ---------------------
# CONVERT TO MATRIX FORMAT
# ---------------------
python3 - <<EOF
import pandas as pd
df = pd.read_csv("$long_csv")
df_wide = df.pivot(index="Gene", columns="Species", values="Distance")
df_wide.to_csv("$wide_csv", float_format="%.5f", na_rep="NA")
EOF

# ---------------------
# CLEANUP
# ---------------------
rm -f temp.fa temp_ref.fa temp_query.fa temp_pair.fa temp_aligned.fa  # Clean up temporary files

echo "✅ Done! Distance matrix written to $wide_csv"
