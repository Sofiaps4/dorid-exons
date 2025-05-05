#!/bin/bash
#The code copies the target gene into the bait set that has been used for the extraction of each of the genes for each species.

# Output file to store species, gene, and FASTA header
output="species_gene_headers.txt"

> "$output"  # Clear the output file

# Loop through each species folder
for species_dir in */; do
  cd "$species_dir" || continue
  species_name="${species_dir%%/}"

  # Loop through each gene folder
  for gene_dir in */; do
    cd "$gene_dir" || continue
    gene_name="${gene_dir%%/}"

    # Look for the _target.fasta file
    target_fasta=$(ls *_target.fasta 2>/dev/null)

    if [ -n "$target_fasta" ]; then
      # Extract the first header line
      header=$(grep '^>' "$target_fasta" | head -n 1)
      echo -e "${species_name}\t${gene_name}\t${header}" >> ../../"$output"


    fi

    cd ..
  done

  cd ..
done
