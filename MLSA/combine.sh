#!/bin/bash

# Loop through each locus references for mlsa

for file in mlsa_loci/*.fas; do
	locus=$(basename "$file" .fas)
	cat ./seq_extract/*$locus.tsv.fasta > ./combined_fasta/${locus}_combined.fasta
done
