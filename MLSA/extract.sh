#!/bin/bash

blast_results="./blast_results"
assemblies="./assembly_list"
loci="./mlsa_loci"
outdir="./seq_extract"

for tsv in blast_results/*.tsv; do
	fname=$(basename "$tsv")
	sample="${fname%_*.*}"
	locus="${fname##*_}"
	asm="$assemblies/$sample.fna"

	# Read the first BLAST hit (assume tab-separated columns)
	read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore < <(head -n 1 "$tsv")

	# Debug printout:
	# echo "DEBUG -- sample=$sample locus=$locus sseqid=$sseqid sstart=$sstart send=$send"

	# Skip if no valid hit information
	if [[ -z "$sseqid" || -z "$sstart" || -z "$send" ]]; then
		echo "No hit or incomplete BLAST output for $fname"
		continue
	fi

	# Ensure coordinates are integers
	if [[ ! "$sstart" =~ ^[0-9]+$ || ! "$send" =~ ^[0-9]+$ ]]; then
		echo "BLAST coordinates are not integers for $fname (sstart=$sstart, send=$send)"
		continue
	fi

	# Choose region, and reverse if needed
	if [ "$sstart" -le "$send" ]; then
		region="$sstart:$send"
		seqkit subseq --chr "$sseqid" -r "$region" "$asm" > "$outdir/${sample}_${locus}.fasta"
	else
		region="$send:$sstart"
		seqkit subseq --chr "$sseqid" -r "$region" "$asm" | seqkit seq -r > "$outdir/${sample}_${locus}.fasta"
	fi

done
