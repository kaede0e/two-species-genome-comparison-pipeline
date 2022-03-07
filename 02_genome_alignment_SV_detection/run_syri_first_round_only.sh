#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-gowens
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G

genus=$1
echo $genus

#0 Give paths to the command and load modules
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin2/anchorwave/minimap2/
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/syri
export PATH=$PATH:/home/kaedeh/scratch/paired_genome_for_syri

module load scipy-stack
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/python_env/bin/activate
PATH_TO_SYRI="/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/syri/syri/bin/syri"
PATH_TO_PLOTSR="/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/syri/syri/bin/plotsr"

python3 $PATH_TO_SYRI \
-c ${genus}_genome/minimap2_ref_qrygenome_chr_alignment.sam \
-r ${genus}_genome/*refgenome_chr.fa \
-q ${genus}_genome/qrygenome_chr_rev.fa \
--dir ${genus}_genome \
-k -F S 2> ${genus}_genome/syri_error_out.txt;
echo "Done $genus first SyRI comparison"
