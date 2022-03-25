 #!/bin/bash

while IFS= read -r file
do
  echo -e "#!/bin/bash" >> run_EDTA_by_chr_${file}.sh
  echo -e "#SBATCH --time=1-00:00:00" >> run_EDTA_by_chr_${file}.sh
  echo -e "#SBATCH --account=YOURACCOUNT" >> run_EDTA_by_chr_${file}.sh
  echo -e "#SBATCH --ntasks=1" >> run_EDTA_by_chr_${file}.sh
  echo -e "#SBATCH --cpus-per-task=16" >> run_EDTA_by_chr_${file}.sh
  echo -e "#SBATCH --mem-per-cpu=2000M" >> run_EDTA_by_chr_${file}.sh
  echo -e "singularity exec -B /home -B /project -B /scratch -B /localscratch /PATH/TO/bin/EDTA/EDTA-2.0.0.sif perl /PATH/TO/bin/EDTA/EDTA_raw.pl --genome ${file}.fwd.fa --cds Prudul26A.cds.fa --exclude Prunus_dulcis_cds.bed --overwrite 1 --sensitive 1 --anno 1 --threads 16" >> run_EDTA_by_chr_${file}.sh
  echo -e "singularity exec -B /home -B /project -B /scratch -B /localscratch /PATH/TO/bin/EDTA/EDTA-2.0.0.sif perl /PATH/TO/bin/EDTA/EDTA.pl --genome ${file}.fwd.fa --cds Prudul26A.cds.fa --exclude Prunus_dulcis_cds.bed --overwrite 0 --sensitive 1 --anno 1 --threads 16" >> run_EDTA_by_chr_${file}.sh
done < "chrnames.txt"
