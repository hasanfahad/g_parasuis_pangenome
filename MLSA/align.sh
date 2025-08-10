#!/bin/bash

for file in combined_fasta/*.fasta; do
	locus=$(basename "$file" .fasta)
	mafft --auto "$file" > aligned_fasta/${locus}_aligned.fasta
done
