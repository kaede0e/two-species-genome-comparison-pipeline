#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --account=def-gowens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000M

#BEFORE DOING THIS, MAKE SURE I HAVE SPACE IN MY ACCOUNT
for FOLDER in $(find . -maxdepth 1 -type d | tail -n +2); do
  echo -ne "$FOLDER:\t"
  find $FOLDER -type f | wc -l
done #This prints how many files there are in your current directory, by subdirectory.
diskusage_report #check how much more your account can support.

singularity exec -B /home -B /project -B /scratch -B /localscratch /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA-2.0.0.sif \
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA.pl \
--genome refgenome_chr.fa \
--cds genome_cds.fasta \
--exclude genome_cds.bed \
--overwrite 1 --sensitive 1 --anno 1 --threads ${SLURM_CPUS_PER_TASK}
  #Singularity exec basically lets the node enter this virtual shell (Singularity img. file),
  #runs one command, and exits the shell automatically once the command is completed.
  #You can change --overwrite 0 to allow running from where you left off - when a job timed out.
