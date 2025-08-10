#!/bin/bash

# Running BLAST for each locus against each assembly
for fasta in mlsa_loci/*.fas
do
	locus=$(basename "$fasta" .fas)
	for asm in assembly_list/*.fna
	do
		sample=$(basename "$asm" .fna)
		# BLAST allele fasta vs assembly
		blastn -query "$fasta" \
		-db "$asm" \
		-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
		-max_target_seqs 1 \
		-num_threads 8 \
		> blast_results/${sample}_${locus}.tsv
	done
done
