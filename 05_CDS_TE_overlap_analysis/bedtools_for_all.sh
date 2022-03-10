#!/bin/bash
# bedtools for getting stats around overlapping sequences between two files
module load bedtools/2.30.0 #current version (Feb 16, 2022)
##For each feature defined in "A", how many features in "B" are overlapping? -c
##For each feature defined in "A", how many basepairs are overlapping? -wo
#For the basepair count, we should merge the CDS first so that there is no duplicates/isoforms that overrepresent those CDS.

##First make bed file out of inversion/syntenic table for each genus
cut -f 1-3 invOut_table.txt | uniq | sed '1d' > inversion_regions.bed
cut -f 1-3 synOut_table.txt | uniq | sed '1d' > syntenic_regions.bed
cut -f 1-3 *_refgenome_chr.fa.mod.EDTA.TEanno.bed > genome_{}_TE.bed

#genome_pair.txt lists my parent directory (ie. Citrus_genome, etc.) in rows
for genus in `cat genome_pair.txt`;
do
  sort -k1,1 -k2,2n ${genus}_genome/*_reduced_cds.bed > ${genus}_genome/sorted_genome_cds.bed
  bedtools merge -i ${genus}_genome/sorted_genome_cds.bed > ${genus}_genome/merged_cds.bed
  #1. Inversion
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/*_reduced_cds.bed -c > ${genus}_genome/bedtools_count_cds_inv.txt
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_inv.txt
  #2. Syntenic region
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/*_reduced_cds.bed -c > ${genus}_genome/bedtools_count_cds_syn.txt
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_syn.txt
done
#For the basepair count, we should merge the TE first so that there is no overlaps and overrepresentation.
for genus in `cat genome_pair.txt`;
do
  sort -k1,1 -k2,2n ${genus}_genome/chr_by_chr_fasta/All_chr_EDTA.TEanno.bed > ${genus}_genome/sorted_${genus}_TE.bed
  bedtools merge -i ${genus}_genome/sorted_${genus}_TE.bed > ${genus}_genome/merged_TE.bed
  #1. Inversion (TE)
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/*EDTA.TEanno.bed -c > ${genus}_genome/bedtools_count_TE_inv.txt
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/merged*_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_inv.txt
  #2. Syntenic region (TE)
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/*EDTA.TEanno.bed -c > ${genus}_genome/bedtools_count_TE_syn.txt
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/merged*_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_syn.txt
done
