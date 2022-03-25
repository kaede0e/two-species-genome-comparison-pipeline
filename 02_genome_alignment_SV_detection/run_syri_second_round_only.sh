#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-gowens
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G

genus=$1
echo $genus

#0 Give paths to the command and load modules
export PATH=$PATH:/PATH/TO/anchorwave/minimap2/
export PATH=$PATH:/PATH/TO/bin
export PATH=$PATH:/PATH/TO/bin/syri
export PATH=$PATH:/PATH/TO/YOURDATA

module load scipy-stack
source /PATH/TO/python_env/bin/activate
PATH_TO_SYRI="/PATH/TO/bin/syri/syri/bin/syri"
PATH_TO_PLOTSR="/PATH/TO/bin/syri/syri/bin/plotsr"

#1 reverse chromosomes from chr_to_rev file and redo minimap
  rm -r ${genus}_genome/qrygenome_chr ;
  mkdir ${genus}_genome/qrygenome_chr ;
  cat ${genus}_genome/*qrychr_names.txt | grep -v -f ${genus}_genome/chr_to_rev.txt > ${genus}_genome/chr_to_keep.txt
  # Reverse chromosomes that cuased problem
  for chr in `cat ${genus}_genome/chr_to_keep.txt` ;
  do
    echo $chr > ${genus}_genome/qrygenome_chr/$chr.txt ;
    samtools faidx ${genus}_genome/*qrygenome_chr.fa -r ${genus}_genome/*qrygenome_chr/$chr.txt > ${genus}_genome/qrygenome_chr/$chr.fwd.fa ;
  done #forward strand fasta for each chromosome
  for chr in `cat ${genus}_genome/chr_to_rev.txt` ;
  do
    echo $chr > ${genus}_genome/qrygenome_chr/rev_$chr.txt ;
    samtools faidx ${genus}_genome/*qrygenome_chr.fa -i -r ${genus}_genome/*qrygenome_chr/rev_$chr.txt > ${genus}_genome/qrygenome_chr/$chr.rev.fa ;
  done #reverse strand fasta for each chromosome
  cat ${genus}_genome/qrygenome_chr/* > ${genus}_genome/qrygenome_chr_rev.fa ;
  minimap2 -ax asm5 --eqx ${genus}_genome/*refgenome_chr.fa ${genus}_genome/qrygenome_chr_rev.fa > ${genus}_genome/minimap2_ref_qrygenome_chr_rev_alignment.sam

#2 second SyRI
        python3 $PATH_TO_SYRI \
  -c ${genus}_genome/minimap2_ref_qrygenome_chr_rev_alignment.sam \
  -r ${genus}_genome/*refgenome_chr.fa \
  -q ${genus}_genome/qrygenome_chr_rev.fa \
  --dir ${genus}_genome \
  -k -F S 2> ${genus}_genome/syri_error2_out.txt
echo "Done $genus second SyRI comparison"
