#!/bin/bash
# bedtools for getting stats around overlapping sequences between two files (bed)
module load bedtools/2.30.0 #current version (Feb 16, 2022)
##For each feature defined in "A", how many features in "B" are overlapping? -c
##For each feature defined in "A", how many basepairs are overlapping? -wo

#NOTE: 
#Merge CDS or TE features before doing bedtools intersect -wo so that there is no duplicates/isoforms that could lead to overrepresentation.
#bedtools intersect -wo does not output 0 bp overlap regions. 

#1) Comparison between inverted region vs. syntenic region
## Make bed file out of inversion/syntenic table for each genus
cut -f 1-3 invOut_table.txt | uniq | sed '1d' > inversion_regions.bed
cut -f 1-3 synOut_table.txt | uniq | sed '1d' > syntenic_regions.bed
## Reduce the bed file for CDS and TE to the first three column only
cut -f 1-3 ${genus}_genome_cds.bed > ${genus}_reduced_cds.bed
cut -f 1-3 */chr_by_chr_fasta/All_chr_EDTA.TEanno.bed > */All_chr_EDTA.TEanno.bed

for genus in `cat genome_pair.txt`; #genome_pair.txt lists my parent directory (ie. Acer, Actinidia, Acer, etc.) in rows
do
  sort -k1,1 -k2,2n ${genus}_genome/*_reduced_cds.bed > ${genus}_genome/sorted_genome_cds.bed
  bedtools merge -i ${genus}_genome/sorted_genome_cds.bed > ${genus}_genome/merged_cds.bed
  #1. Inversion (CDS)
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/*_reduced_cds.bed -c > ${genus}_genome/bedtools_count_cds_inv.txt
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_inv.txt
  #2. Syntenic region (CDS)
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/*_reduced_cds.bed -c > ${genus}_genome/bedtools_count_cds_syn.txt
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_syn.txt
  
  sort -k1,1 -k2,2n ${genus}_genome/All_chr_EDTA.TEanno.bed > ${genus}_genome/sorted_${genus}_TE.bed
  bedtools merge -i ${genus}_genome/sorted_${genus}_TE.bed > ${genus}_genome/merged_TE.bed
  #1. Inversion (TE)
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/*EDTA.TEanno.bed -c > ${genus}_genome/bedtools_count_TE_inv.txt
  bedtools intersect -a ${genus}_genome/inversion_regions.bed -b ${genus}_genome/merged*_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_inv.txt
  #2. Syntenic region (TE)
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/*EDTA.TEanno.bed -c > ${genus}_genome/bedtools_count_TE_syn.txt
  bedtools intersect -a ${genus}_genome/syntenic_regions.bed -b ${genus}_genome/merged*_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_syn.txt
done


#2) Comparison between inversion breakpoint regions vs. randomized regions in the reference genome 
##Create a randomized regions file
---------------------------------------
##Option a: randomization using SyRI output
#1. Merge all aligned regions identified by SyRI between the genome pair into full length (syri_aligned_regions.bed)
for genus in `cat genome_pair.txt`;
do
  ## remove unaligned ("NOTAL") regions and inversions ("INV") that are identified in syri.out
  cat ${genus}_genome/SyRI_output_${genus}_genome/syri.out | grep -v "NOTAL" | grep -v "INV" > ${genus}_genome/SyRI_output_${genus}_genome/aligned_syri.out
  cut -f 1-3 ${genus}_genome/SyRI_output_${genus}_genome/aligned_syri.out > ${genus}_genome/SyRI_output_${genus}_genome/syri_aligned_regions.bed
  cat ${genus}_genome/SyRI_output_${genus}_genome/syri_aligned_regions.bed | awk '{if ($2 <= $3){print}}' > ${genus}_genome/SyRI_output_${genus}_genome/valid_syri_aligned_regions.bed
  sort -k1,1 -k2,2n ${genus}_genome/SyRI_output_${genus}_genome/valid_syri_aligned_regions.bed > ${genus}_genome/SyRI_output_${genus}_genome/sorted_syri_aligned_regions.bed
  bedtools merge -i ${genus}_genome/SyRI_output_${genus}_genome/sorted_syri_aligned_regions.bed > ${genus}_genome/SyRI_output_${genus}_genome/merged_syri_aligned_regions.bed
  cp ${genus}_genome/SyRI_output_${genus}_genome/merged_syri_aligned_regions.bed /home/kaedeh/scratch/paired_genome_done_EDTA/${genus}_genome/
  rm */SyRI_output_*/sorted_syri_aligned_regions.bed */SyRI_output_*/syri_aligned_regions.bed */SyRI_output_*/aligned_syri.out
done
#2. Randomly pick regions from SyRI aligned regions (merged_syri_aligned_regions.bed)
##  -n adjust to 1000 x number of chromosomes
for genus in `cat genome_pair.txt`;
do
  cat ${genus}_genome/merged_syri_aligned_regions.bed | cut -f 1 | uniq -d | wc -l > ${genus}_genome/n.txt
  bedtools random -n `cat ${genus}_genome/n.txt`000 -l 4000 -g ${genus}_genome/merged_syri_aligned_regions.bed | cut -f 1-3 | sort -k1,1 -k2,2n > ${genus}_genome/bedtools_random_4k_aligned_regions.bed
done

##Option b: randomization using whole chromosome 
#1. Prepare chromosome size file from .fai
samtools faidx refgenome.fa
#2. Randomly pick regions from chromosome file
##  -n adjust to 1000 x number of chromosomes
for genus in `cat genome_pair.txt`;
do
  cat ${genus}_genome/*chr.txt | tail -n +2 | tr ' ' '  ' | awk '{print $1,$3}' > ${genus}_genome/chr_size.txt
  cat ${genus}_genome/chr_size.txt | wc -l > n.txt
  cat ${genus}_genome/chr_size.txt | tr ' ' '  ' > ${genus}_genome/input_chr.txt
  bedtools random -n `cat ${genus}_genome/n.txt`000 -l 4000 -g ${genus}_genome/input_chr.txt | cut -f 1-3 | sort -k1,1 -k2,2n > ${genus}_genome/bedtools_random_4k_regions.bed
done
----------------------------

## Extract inversion breakpoints by R (refer to extract_inv_breakpoints.R) and perform overlap analysis using bedtools intersect
for genus in `cat genome_pair.txt`;
do
  #0. Convert breakpoints_regions.csv to .bed
  cat ${genus}_genome/breakpoints_regions_${genus}.csv | tail -n +2 | cut -f 3-5 -d "," | tr -d '"' | sed s/","/"\t"/g > ${genus}_genome/breakpoints_regions.bed
  #1. Inversion breakpoints
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b TE_files/sorted_${genus}_TE.bed -c > ${genus}_genome/bedtools_count_TE_inv_breakpoints.txt
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b ${genus}_genome/merged_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_inv_breakpoints.txt
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b ${genus}_genome/sorted_genome_cds.bed -c > ${genus}_genome/bedtools_count_cds_inv_breakpoints.txt
  bedtools intersect -a ${genus}_genome/breakpoints_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_inv_breakpoints.txt
  #2a. Randomized 4kbp alignable regions in the genome
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b TE_files/sorted_${genus}_TE.bed -c > ${genus}_genome/bedtools_count_TE_random_4k.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b ${genus}_genome/merged_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_random_4k.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b ${genus}_genome/sorted_genome_cds.bed -c > ${genus}_genome/bedtools_count_cds_random_4k.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_aligned_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_random_4k.txt
  #2b. Randomized 4kbp regions in the genome
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_regions.bed -b TE_files/sorted_${genus}_TE.bed -c > ${genus}_genome/bedtools_count_TE_random_4k_wholechr.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_regions.bed -b ${genus}_genome/merged_TE.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_TE_random_4k_wholechr.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_regions.bed -b ${genus}_genome/sorted_genome_cds.bed -c > ${genus}_genome/bedtools_count_cds_random_4k_wholechr.txt
  bedtools intersect -a ${genus}_genome/bedtools_random_4k_regions.bed -b ${genus}_genome/merged_cds.bed -wo > ${genus}_genome/bedtools_bpoverlap_merged_cds_random_4k_wholechr.txt
done

#3) Comparison between breakpoints vs. random points in the genome for frequency of containing a gene
## Extract breakpoints coordinates in .bed format
for genus in `cat genome_pair.txt`;
do
  cat ${genus}_genome/inversion_regions.bed | cut -f -2 > ${genus}_genome/inversion_start_regions.txt
  cat ${genus}_genome/inversion_regions.bed | awk '{print $1,$3}' | tr ' ' '      ' > ${genus}_genome/inversion_end_regions.txt
  cat ${genus}_genome/inversion_start_regions.txt ${genus}_genome/inversion_end_regions.txt >> ${genus}_genome/inversion_breakpoints.txt
  awk '{print $0,$NF}' ${genus}_genome/inversion_breakpoints.txt | tr ' ' '       ' > ${genus}_genome/inversion_breakpoints.bed
done

## Use bedtools intersect -c (count) to see if breakpoints overlaps with any genes
for genus in `cat genome_pair.txt`;
do
  cat ${genus}_genome/*genes.bed | cut -f 1-3 | sort -k1,1 -k2,2n > ${genus}_genome/sorted_genes.bed
  bedtools merge -i ${genus}_genome/sorted_genes.bed > ${genus}_genome/merged_genes.bed
  bedtools intersect -a ${genus}_genome/inversion_breakpoints.bed -b ${genus}_genome/merged_genes.bed -c > ${genus}_genome/bedtools_count_gene_at_breakpoints.txt
done

## To have a comparison to the above stats, extract random points from the genome within aligned regions
## Extract random point coordinates in .bed format
for genus in `cat genome_pair.txt`;
do
  cat ${genus}_genome/bedtools_random_4k_aligned_regions.bed | cut -f -2 > ${genus}_genome/random_aligned_start_regions.txt
  cat ${genus}_genome/bedtools_random_4k_aligned_regions.bed | awk '{print $1,$3}' | tr ' ' '      ' > ${genus}_genome/random_aligned_end_regions.txt
  cat ${genus}_genome/random_aligned_start_regions.txt ${genus}_genome/random_aligned_end_regions.txt >> ${genus}_genome/random_aligned_points.txt
  awk '{print $0,$NF}' ${genus}_genome/random_aligned_points.txt | tr ' ' '        ' > ${genus}_genome/random_aligned_points.bed
done

## Use bedtools intersect -c (count) to see if random points overlap with any genes
for genus in `cat genome_pair.txt`;
do
  bedtools intersect -a ${genus}_genome/random_aligned_points.bed -b ${genus}_genome/merged_genes.bed -c > ${genus}_genome/bedtools_count_gene_at_random_points.txt
done


