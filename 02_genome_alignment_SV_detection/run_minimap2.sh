#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=YOURACCOUNT
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=6G

export PATH=$PATH:/PATH/TO/anchorwave/minimap2/

minimap2 -ax asm5 --eqx refgenome.fa qrygenome.fa > minimap2_whole_genome_alignment.sam
