# Pan-genomic Analysis of 121 _Glaesserella parasuis_ Strains
## Overview
This repository contains the computational workflow and analysis methods used to characterize 121 Glaesserella parasuis strains through comprehensive genomic approaches.

## Results

### Multilocus Sequence Analysis (MLSA) Phylogeny tree
<img width="5432" height="5432" alt="mlsa_tree2" src="https://github.com/user-attachments/assets/8c134eb0-5b85-4543-a9cd-13894490919f" />

### Average Neucleotide Identity (ANI) Heatmap
<img width="2400" height="1800" alt="ani_plot" src="https://github.com/user-attachments/assets/67be1ea3-4a2d-4fdb-b9be-2e4898d0fa6e" />
<img width="2400" height="1800" alt="ani_plot_2" src="https://github.com/user-attachments/assets/af42d427-c64d-4e72-9023-c38f92f5ba43" />



## Dataset
* **Species**: *Glaesserella parasuis*
* **Number of Strains**: 121 (105 from NCBI Genbank, 16 from China National Genbank)
* **Data type:** Whole genome assemblies with GFF3 Annotations
* **Source 1:** https://www.ncbi.nlm.nih.gov/genbank/
* **Source 2:** https://db.cngb.org/

# Methods
## MLSA Phylogeny
### Housekeeping gene selection
The following seven housekeeping genes were selected based on:
* Previous G. parasuis MLSA studies
* Single-copy nature across bacterial genomes
* Adequate phylogenetic signal
* Primer availability for experimental validation

| Gene  | Product | Function |
| ------------- | ------------- | -------- |
|atpD|ATP synthase subunit beta|Energy metabolism|
|gyrB|DNA gyrase subunit B|DNA replication|
|recA|Recombinase A|DNA repair|
|rpoB|RNA polymerase subunit beta|Transcription|
|dnaJ|Molecular chaperone DnaJ|Protein folding|
|mdh|Malate dehydrogenase|Central metabolism|
|pgi|Glucose-6-phosphate isomerase|Glycolysis|

### Analysis Strategy
* Gene extraction from genome assemblies
* Sequence alignment for each locus
* Alignment concatenation into supermatrix
* Phylogenetic reconstruction with statistical support
* Tree visualization and interpretation

### Tools used
```
bakta >= 1.11.3           # Genome annotation
blast+ >= 2.12.0          # Sequence search
seqkit >= 2.3.0           # Sequence manipulation
mafft >= 7.490            # Multiple sequence alignment
amas >= 1.0.0             # Alignment manipulation
iqtree >= 2.1.0           # Phylogenetic inference
```

### Command
#### Bakta
```
bakta --db bakta_db/db assembly_list/GCA_000444565.1_ASM44456v1_genomic.fna \
	 --output annotations/GCA_000444565.1_ASM44456v1_genomic

```
#### BLAST
```
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
```

#### Combine each assembly with the housekeeping genes
```
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
```

#### Loop through each locus reference for MLSA
```
#!/bin/bash
# Loop through each locus references for mlsa
for file in mlsa_loci/*.fas; do
	locus=$(basename "$file" .fas)
	cat ./seq_extract/*$locus.tsv.fasta > ./combined_fasta/${locus}_combined.fasta
done
```
#### Align sequence with mafft
```
#!/bin/bash

for file in combined_fasta/*.fasta; do
	locus=$(basename "$file" .fasta)
	mafft --auto "$file" > aligned_fasta/${locus}_aligned.fasta
done
```

#### Create MLSA supermatrix using FASconCAT-G_v1.06.1.pl

#### IQTree for calculating MLSA
```
iqtree -s mlsa_supermatrix.fas -m MFP -bb 1000 -nt AUTO
```
#### Visualization using iTOL (https://itol.embl.de/)

## Average Nucleotide Identity
ANI analysis was performed to assess overall genomic similarity between strains and establish phylogenetic relationships.
### Tools Used
* FastANI v1.34
### Command
```
fastANI --ql ani_list.txt --rl ani_list.txt -t 32 -o all_vs_all_fastani_output.txt
```

## Pangenome construction using Roary
### Overview
Roary was used to perform pan-genome analysis on 121 Glaesserella parasuis strains to:

* Identify core, accessory, and unique genes
* Generate gene presence/absence matrices
* Perform phylogenetic analysis
* Enable pan-GWAS studies for virulence traits

### Tools used
* Rorry v3.13.0

### Command
```
roary -e -n -v -p 32 gff3/*.gff3
```
## Contact
Created by Fahad Hasan <br>
**PhD Student, Department of Genomics, Phenomics, and Bioinformatics** <br>
North Dakota State University

For questions about this pipeline, please contact [mdfahad.hasan@ndsu.edu]
