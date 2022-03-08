#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --account=def-gowens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000M

singularity exec -B /home -B /project -B /scratch -B /localscratch /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA-2.0.0.sif perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA_raw.pl --genome Pd01.fwd.fa --cds Prudul26A.cds.fa --exclude Prunus_dulcis_cds.bed --overwrite 1 --sensitive 1 --anno 1 --threads 16
singularity exec -B /home -B /project -B /scratch -B /localscratch /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA-2.0.0.sif perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA.pl --genome Pd01.fwd.fa --cds Prudul26A.cds.fa --exclude Prunus_dulcis_cds.bed --overwrite 0 --sensitive 1 --anno 1 --threads 16
