#!/bin/bash

#Finding segmental duplications or inverted segments in breakpoint regions upstream/downstream of inversions
module load bedtools
module load samtools
module load StdEnv/2020 gcc/9.3.0 blast+/2.12.0

#This can extract sequences specified in bed file
bedtools getfasta -fi ../*refgenome_chr.fa -bed ../breakpoints_regions_onesided_10k.bed > ../refgenome_breakpoints_regions_onesided_10k.fa

cat ../breakpoints_regions_onesided_10k.bed | head -79 > breakpoints_regions_onesided_start10k.bed #adjust the number to number_of_filtered_inversions = 1/2 of the row number of this breakpoints_regions_onesided_10k.bed file
cat ../breakpoints_regions_onesided_10k.bed | tail -79 > breakpoints_regions_onesided_end10k.bed
split -l 1 -a 3 -d --additional-suffix=.start.txt breakpoints_regions_onesided_start10k.bed
split -l 1 -a 3 -d --additional-suffix=.end.txt breakpoints_regions_onesided_end10k.bed

for i in {000..78};#change the numbers to match the number of inversions you have -1
do
  bedtools getfasta -fi ../*refgenome_chr.fa -bed x$i.start.txt > refgenome_inv_x$i.start.fa
  bedtools getfasta -fi ../*refgenome_chr.fa -bed x$i.end.txt > refgenome_inv_x$i.end.fa
  makeblastdb -in refgenome_inv_x$i.start.fa -dbtype nucl
  blastn -query refgenome_inv_x$i.end.fa -db refgenome_inv_x$i.start.fa -out refgenome_inv_10k_startx${i}_to_endx${i}.txt -outfmt 6
  rm refgenome_inv_x$i.start.fa*
  rm refgenome_inv_x$i.end.fa*
done
cat refgenome_inv_10k_startx*_to_endx*.txt > refgenome_inv_10k_blast_score_aln.txt

for genus in `cat genome_pair.txt`;
do
  #cat ${genus}_genome/breakpoints_regions_onesided_10k.bed | wc | awk '{print $1*0.5}' > ${genus}_genome/number_of_filtered_inversions.txt
  cat ${genus}_genome/fixed_breakpoints_regions_onesided.bed | tail -`cat ${genus}_genome/number_of_filtered_inversions.txt` > ${genus}_genome/start_breakpoints_regions_onesided.bed
  cat ${genus}_genome/fixed_breakpoints_regions_onesided_10k.bed | tail -`cat ${genus}_genome/number_of_filtered_inversions.txt` > ${genus}_genome/start_breakpoints_regions_onesided_10k.bed
done
