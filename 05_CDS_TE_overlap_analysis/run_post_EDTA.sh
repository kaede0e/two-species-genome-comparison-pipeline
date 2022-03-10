#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --account=def-gowens
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=2000M

export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin
module load bedops/2.4.39
module load samtools

# 4) After completion of EDTA
#Loops for getting TE.bed and archiving and deleting unnecessary files. Call this script in your genome chr_by_chr_fasta/:
mkdir completed_EDTA
for file in `cat chrnames.txt`;
do
  gff2bed < ${file}.fwd.fa.mod.EDTA.TEanno.gff3 > completed_EDTA/${file}.fwd.fa.mod.EDTA.TEanno.bed
  tar -cvf ${file}.fwd.fa.mod.EDTA.final.tar ${file}.fwd.fa.mod.EDTA.final/ > ${file}.fwd.fa.mod.EDTA.final_zipped_files.txt
  rm -r ${file}.fwd.fa.mod.EDTA.final/
done
#Concatenate all chromosomes TE results
echo "Concatenate all chromosomes TE results"
cat completed_EDTA/*.mod.EDTA.TEanno.bed >> All_chr_EDTA.TEanno.bed
