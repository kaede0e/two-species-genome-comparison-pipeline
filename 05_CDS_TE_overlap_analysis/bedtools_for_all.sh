#!/bin/bash
# bedtools for getting stats around overlapping sequences between two files
module load bedtools/2.30.0 #current version (Feb 16, 2022)
##For each feature defined in "A", how many features in "B" are overlapping? -c
##For each feature defined in "A", how many basepairs are overlapping? -wo
#For the basepair count, we should merge the CDS first so that there is no duplicates/isoforms that overrepresent those CDS.

#1) Comparison between inverted region vs. syntenic region
## Make bed file out of inversion/syntenic table for each genus
cut -f 1-3 invOut_table.txt | uniq | sed '1d' > inversion_regions.bed
cut -f 1-3 synOut_table.txt | uniq | sed '1d' > syntenic_regions.bed
cut -f 1-3 *_refgenome_chr.fa.mod.EDTA.TEanno.bed > genome_{}_TE.bed

for genus in `cat genome_pair.txt`; #genome_pair.txt lists my parent directory (ie. Citrus_genome, etc.) in rows
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
## For the basepair count, we should merge the TE first so that there is no overlaps and overrepresentation.
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


#2) Comparison between inversion breakpoints vs. randomized regions in the reference genome 
### Create a randomized regions file ###
##1 merge all aligned regions between the genome pair into full length (syri_aligned_regions.bed)
for genus in `cat genome_pair.txt`;
do
  #remove unaligned ("NOTAL") regions, indels ("INS", "DEL") that are identified in syri.out
  cat ${genus}_genome/SyRI_output_${genus}_genome/syri.out | grep -v "NOTAL" | grep -v "SNP" | grep -v "DEL" | grep -v "INS" > ${genus}_genome/SyRI_output_${genus}_genome/aligned_syri.out
  cut -f 1-3 ${genus}_genome/SyRI_output_${genus}_genome/aligned_syri.out > ${genus}_genome/SyRI_output_${genus}_genome/syri_aligned_regions.bed
  sort -k1,1 -k2,2n ${genus}_genome/SyRI_output_${genus}_genome/syri_aligned_regions.bed > ${genus}_genome/SyRI_output_${genus}_genome/sorted_syri_aligned_regions.bed
  bedtools merge -i ${genus}_genome/SyRI_output_${genus}_genome/sorted_syri_aligned_regions.bed > ${genus}_genome/SyRI_output_${genus}_genome/merged_syri_aligned_regions.bed
  cp ${genus}_genome/SyRI_output_${genus}_genome/merged_syri_aligned_regions.bed /home/kaedeh/scratch/paired_genome_done_EDTA/${genus}_genome/
  rm */SyRI_output_*/sorted_syri_aligned_regions.bed */SyRI_output_*/syri_aligned_regions.bed */SyRI_output_*/aligned_syri.out
done
##2 randomly pick regions from SyRI aligned regions (merged_syri_aligned_regions.bed)
#-n adjust to 1000 x number of chromosomes
bedtools random -n 19000 -l 4000 -g merged_syri_aligned_regions.bed | cut -f 1-3 | sort -k1,1 -k2,2n > bedtools_random_4k_aligned_regions.bed

##3 extract inversion breakpoints by R and bedtools intersect
for genus in `cat invalid_genome_pair.txt`;
do
  #0. Convert breakpoints_regions.csv to .bed
  cat ${genus}_genome/breakpoints_regions_${genus}.csv | tail -n +2 | cut -f 3-5 -d "," | tr -d '"' | sed s/","/"\t"/g > ${genus}_genome/breakpoints_regions.bed
  #1. Inversion breakpoints
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b TE_files/sorted_${genus}_TE.bed -c > ${genus}_genome/bedtools_count_TE_inv_breakpoints.txt
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b ${genus}_genome/merged_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_inv_breakpoints.txt
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b ${genus}_genome/sorted_genome_cds.bed -c > ${genus}_genome/bedtools_count_cds_inv_breakpoints.txt
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_inv_breakpoints.txt
  #2. Randomized 4kbp alignable regions in the genome
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b TE_files/sorted_${genus}_TE.bed -c > ${genus}_genome/bedtools_count_TE_random_4k.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b ${genus}_genome/merged_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_random_4k.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b ${genus}_genome/sorted_genome_cds.bed -c > ${genus}_genome/bedtools_count_cds_random_4k.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_random_4k.txt
done
