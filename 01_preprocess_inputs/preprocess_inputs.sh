#!/bin/bash
module load samtools #used for retrieving only chromosomes 
module load bedops/2.4.39 #check the version before downloading; used for converting gff3 to bed


#1) Process reference and query genome
#1.1) Remove fragments and retain chromosomes only
cat refgenome.fa | grep ">chr" | sed s/">"//g > refgenome_chr_names.txt
samtools faidx refgenome.fa -r refgenome_chr_names.txt > refgenome_chr.fa
cat qrygenome.fa | grep ">chr" | sed s/">"//g > qrygenome_chr_names.txt
samtools faidx qrygenome.fa -r qrygenome_chr_names.txt > qrygenome_chr.fa
#1.1.2) If the fasta files have messy headline, remove spacer and pick a unique column (often the case for genomes donwnloaded directly from NCBI)
cat refgenome.fa | cut -f 1 -d " " > ref.fa #This separates by a space, retains the first column.
cat refgenome.fa | grep ">CM" | cut -f 1 -d " " | sed s/">"//g > ref_chr_names.txt

#2) Process gene annotation file 
## For converting gff to bed files, the format has to be:
##column 1: name of chromosome
##column 2: start of chromosome
##column 3: end of chromosome
## and no header. 
#2.1) Cut only CDS rows 
awk '$3 ~ /CDS/' genome.gff > genome_cds.gff
#2.2) Convert gff3 into bed
gff2bed < genome_cds.gff > genome_cds.bed
