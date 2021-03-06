#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=YOURACCOUNT
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G

genus=$1
echo $genus

#0 Give paths to the command and load modules
export PATH=$PATH:/PATH/TO/anchorwave/minimap2/
export PATH=$PATH:/PATH/TO/bin
export PATH=$PATH:/PATH/TO/bin/syri
export PATH=$PATH:/PATH/TO/YOURDATA

#1 Whole genome alignment round 1
minimap2 -t 10 -ax asm5 --eqx ${genus}_genome/*refgenome_chr.fa ${genus}_genome/*qrygenome_chr.fa > ${genus}_genome/minimap2_ref_qrygenome_chr_alignment.sam

#2 Run SyRI round 1
module load scipy-stack
source /PATH/TO/python_env/bin/activate
PATH_TO_SYRI="/PATH/TO/bin/syri/syri/bin/syri"
PATH_TO_PLOTSR="/PATH/TO/bin/syri/syri/bin/plotsr"

  python3 $PATH_TO_SYRI \
-c ${genus}_genome/minimap2_ref_qrygenome_chr_alignment.sam \
-r ${genus}_genome/*refgenome_chr.fa \
-q ${genus}_genome/*qrygenome_chr.fa \
--dir ${genus}_genome \
-k -F S 2> ${genus}_genome/syri_error_out.txt ;
  cat ${genus}_genome/syri_error_out.txt | grep WARN  | grep "high fraction of inverted" | cut -f 23 -d " " | sed s/\(//g | sed s/\)\.//g > ${genus}_genome/chr_to_rev.txt  ;
  echo "Done first $genus comparison" ;
  echo "Done first round of SyRI" ;
##when first round of SyRI fails, it's due to complementary sequences of chromosomes comapared. The following complements those problematic chromosomes identified in "chr_to_rev.txt". 
  rm -r ${genus}_genome/qrygenome_chr ;
  mkdir ${genus}_genome/qrygenome_chr ;
  cat ${genus}_genome/*qrychr_names.txt | grep -v -f ${genus}_genome/chr_to_rev.txt > ${genus}_genome/chr_to_keep.txt
  for chr in `cat ${genus}_genome/chr_to_keep.txt` ;
  do
    echo $chr > ${genus}_genome/qrygenome_chr/$chr.txt ;
    samtools faidx ${genus}_genome/*qrygenome_chr.fa -r ${genus}_genome/qrygenome_chr/$chr.txt > ${genus}_genome/qrygenome_chr/$chr.fwd.fa ;
  done #forward strand fasta for each chromosome
  for chr in `cat ${genus}_genome/chr_to_rev.txt` ;
  do
    echo $chr > ${genus}_genome/qrygenome_chr/rev_$chr.txt ;
    samtools faidx ${genus}_genome/*qrygenome_chr.fa -i -r ${genus}_genome/qrygenome_chr/rev_$chr.txt > ${genus}_genome/qrygenome_chr/$chr.rev.fa ;
  done #reverse strand fasta for each chromosome
  cat ${genus}_genome/qrygenome_chr/* > ${genus}_genome/qrygenome_chr_rev.fa ;
  echo "Done first round of SyRI"
  
#3 Redo Whole Genome Alignment on reversed genome
minimap2 -t 10 -ax asm5 --eqx ${genus}_genome/*refgenome_chr.fa ${genus}_genome/qrygenome_chr_rev.fa > ${genus}_genome/minimap2_ref_qrygenome_chr_rev_alignment.sam

#4 Run SyRI round 2
    python3 $PATH_TO_SYRI \
  -c ${genus}_genome/minimap2_ref_qrygenome_chr_rev_alignment.sam \
  -r ${genus}_genome/*refgenome_chr.fa \
  -q ${genus}_genome/qrygenome_chr_rev.fa \
  --dir ${genus}_genome \
  -k -F S 2> ${genus}_genome/syri_error2_out.txt
  echo "Done $genus second SyRI comparison" ;
  echo "Done second round of SyRI"
  
 #5 Format SyRI output into tidy table
cat ${genus}_genome/invOut.txt | perl /PATH/TO/syriout2table.pl > ${genus}_genome/invOut_table.txt
cat ${genus}_genome/synOut.txt | perl /PATH/TO/syriout2table.pl > ${genus}_genome/synOut_table.txt
 
