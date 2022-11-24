#!/bin/bash

#Finding segmental duplications in inversion breakpoint regions
module load bedtools
module load samtools
module load StdEnv/2020 gcc/9.3.0 blast+/2.12.0

##### Duplication detection using BLAST with smaller pieces (1kbp query, to 10kbp reference on each side)
# this finds duplicates within the 10kbp breakpoint region. 
for genus in `cat genome_pair.txt`;
do
  #0. Convert 1kbp_coordinates.csv to .bed
  mkdir ${genus}_genome/blast_1kbp_start
  mkdir ${genus}_genome/blast_1kbp_end
  mkdir ${genus}_genome/blast_1kbp_start/individual_1kbp_start
  mkdir ${genus}_genome/blast_1kbp_end/individual_1kbp_end
  mkdir ${genus}_genome/blast_1kbp_start/individual_1kbp_start/BLAST_database
  mkdir ${genus}_genome/blast_1kbp_end/individual_1kbp_end/BLAST_database
  cat ${genus}_genome/breakpoints_start_1k_coordinates_${genus}.csv | tail -n +2 | cut -f 3-5 -d "," | tr -d '"' | sed s/","/"\t"/g > ${genus}_genome/blast_1kbp_start/breakpoints_start_1k_coordinates.bed
  cat ${genus}_genome/breakpoints_end_1k_coordinates_${genus}.csv | tail -n +2 | cut -f 3-5 -d "," | tr -d '"' | sed s/","/"\t"/g > ${genus}_genome/blast_1kbp_end/breakpoints_end_1k_coordinates.bed
  #1. Make separate files for every inversion. 10 x 1kbp in each coordinate file
  cat ${genus}_genome/blast_1kbp_start/breakpoints_start_1k_coordinates.bed | sort > ${genus}_genome/blast_1kbp_start/sorted_breakpoints_start_1k_coordinates.bed
  split -l 10 ${genus}_genome/blast_1kbp_start/sorted_breakpoints_start_1k_coordinates_fixed.bed ${genus}_genome/blast_1kbp_start/1kbp_coordinate_ --additional-suffix=.bed -d -a 3
  cat ${genus}_genome/blast_1kbp_end/breakpoints_end_1k_coordinates.bed | sort > ${genus}_genome/blast_1kbp_end/sorted_breakpoints_end_1k_coordinates.bed
  split -l 10 ${genus}_genome/blast_1kbp_end/sorted_breakpoints_end_1k_coordinates_fixed.bed ${genus}_genome/blast_1kbp_end/1kbp_coordinate_ --additional-suffix=.bed -d -a 3
done

###### From here, submit the following script from each Genus directory. #####
#2. Make individual 1kbp coordinate file. (do it in each genus and double check numbers are right)
  awk '{print > ("blast_1kbp_start/individual_1kbp_start/seq_" NR ".bed")}' blast_1kbp_start/1kbp_coordinate_*.bed
  awk '{print > ("blast_1kbp_end/individual_1kbp_end/seq_" NR ".bed")}' blast_1kbp_end/1kbp_coordinate_*.bed

  #3. Extract sequences for each of those positions.
  for i in {000..225};#change the numbers to match the number of inversions you have
  do
    bedtools getfasta -fi *refgenome_chr.fa -bed blast_1kbp_start/1kbp_coordinate_$i.bed > blast_1kbp_start/1kbp_coordinate_$i.fa #start breakpoints; all 10 x 1kbp sequences
    bedtools getfasta -fi *refgenome_chr.fa -bed blast_1kbp_end/1kbp_coordinate_$i.bed > blast_1kbp_end/1kbp_coordinate_$i.fa #end breakpoints
  done
  for j in {1..2260}; #change to 10x the total inv number
  do
    bedtools getfasta -fi *refgenome_chr.fa -bed blast_1kbp_start/individual_1kbp_start/seq_$j.bed > blast_1kbp_start/individual_1kbp_start/seq_$j.fa #start breakpoints
    bedtools getfasta -fi *refgenome_chr.fa -bed blast_1kbp_end/individual_1kbp_end/seq_$j.bed > blast_1kbp_end/individual_1kbp_end/seq_$j.fa #end breakpoints
  done
  #4. Run BLAST for each individual 1kbp to the corresponding x10 fasta file.
  for h in {000..225};
  do
    makeblastdb -in blast_1kbp_start/1kbp_coordinate_${h}.fa -dbtype nucl -out blast_1kbp_start/individual_1kbp_start/BLAST_database_start_${h}
    makeblastdb -in blast_1kbp_end/1kbp_coordinate_${h}.fa -dbtype nucl -out blast_1kbp_end/individual_1kbp_end/BLAST_database_end_${h}
    for times in {1..10};
    do
      blastn -query blast_1kbp_start/individual_1kbp_start/seq_`expr ${h} \* 10 + ${times}`.fa -db blast_1kbp_start/individual_1kbp_start/BLAST_database_start_${h} -out blast_1kbp_start/individual_1kbp_start/seq_`expr ${h} \* 10 + ${times}`_to_1kbp_start_${h}.txt -outfmt 6
      blastn -query blast_1kbp_end/individual_1kbp_end/seq_`expr ${h} \* 10 + ${times}`.fa -db blast_1kbp_end/individual_1kbp_end/BLAST_database_end_${h} -out blast_1kbp_start/individual_1kbp_start/seq_`expr ${h} \* 10 + ${times}`_to_1kbp_end_${h}.txt -outfmt 6
    done
  done
   #this is run from the _start directory so I need to change that. if it was successful!

#5. Delete self-alignment from each blast.
for files in `ls seq*_to_1kbp_*.txt`; #This must be called from blast_1kbp_start/individual* directory.
do
  cat $files | tail -n +2 > self_aln_filtered_${files}
done

#6. Concatenate results by genus for R import.
cat blast_1kbp_*/individual*/self_aln_filtered_* > self_aln_filtered_1kbp_blast_score.txt

for genus in `cat genome_pair.txt`;
do
  awk '{print var, "  ",$0}' var="${genus}" ${genus}_genome/self_aln_filtered_1kbp_blast_score.txt > ${genus}_self_aln_filtered_1kbp_blast_score.txt;
done

cat *_self_aln_filtered_1kbp_blast_score.txt > Master_inv_breakpoints_within_aln_blast_output_1kb_window.txt
