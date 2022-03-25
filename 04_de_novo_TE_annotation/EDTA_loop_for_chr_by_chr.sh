#!/bin/bash
#!/usr/bin/env python3
export PATH=$PATH:/PATH/TO/bin/EDTA
module load bedops/2.4.39
module load samtools

# 1) Extract individual chromosomes and names file
#This loop creates fasta files for individual chromosomes in *_genome/chr_by_chr_fasta/, chrnames.txt in *_genome/
genome_pair="/home/kaedeh/scratch/paired_genome_done_syri/xxxx_genome"
for genus in $genome_pair
do
  cat ${genus}/*refgenome_chr.fa | grep ">" | sed s/">"//g > ${genus}/chrnames.txt
  mkdir ${genus}/chr_by_chr_fasta
  for chr in `cat ${genus}/chrnames.txt`;
  do
    echo ${genus}/${chr} > ${genus}/chr_by_chr_fasta/${chr}.txt;
    cat ${genus}/chr_by_chr_fasta/${chr}.txt | cut -d "/" -f 7 > ${genus}/chr_by_chr_fasta/${chr}_chrname.txt
    rm ${genus}/chr_by_chr_fasta/${chr}.txt ;
    samtools faidx ${genus}/*refgenome_chr.fa -r ${genus}/chr_by_chr_fasta/${chr}_chrname.txt > ${genus}/chr_by_chr_fasta/${chr}.fwd.fa;
    rm ${genus}/chr_by_chr_fasta/${chr}_chrname.txt ;
  done
  mv ${genus}/chrnames.txt ${genus}/chr_by_chr_fasta/ ;
done

# 2) Make batch of scripts to run EDTA individually
chmod +x run_EDTA_*.sh
source run_EDTA_*.sh
bash ../../while_loop_for_EDTA.sh #Call this from *_genome/chr_by_chr_fasta/ where run_EDTA_*.sh is located in.
                                  #This submits all individual chromosome scripts.

# 3) Script to run EDTA in Singularity
for file in `cat chrnames.txt`;
do
  singularity exec -B /home -B /project -B /scratch -B /localscratch /PATH/TO/bin/EDTA/EDTA-2.0.0.sif \
perl /PATH/TO/bin/EDTA/EDTA_raw.pl \
--genome chr_by_chr_fasta/${file}.fwd.fa \
--cds aradu.V14167.gnm1.ann1.cxSM.cds.fna \
--exclude Arachis_duranensis_chr_renamed_cds.bed \
--overwrite 1 --sensitive 1 --anno 1 --threads ${SLURM_CPUS_PER_TASK}

  singularity exec -B /home -B /project -B /scratch -B /localscratch /PATH/TO/bin/EDTA/EDTA-2.0.0.sif \
perl /PATH/TO/EDTA/EDTA.pl \
--genome chr_by_chr_fasta/${file}.fwd.fa \
--cds aradu.V14167.gnm1.ann1.cxSM.cds.fna \
--exclude Arachis_duranensis_chr_renamed_cds.bed \
--overwrite 0 --sensitive 1 --anno 1 --threads ${SLURM_CPUS_PER_TASK}
done

# 4) After completion of EDTA
#Loops for getting TE.bed and archiving and deleting unnecessary files:
mkdir completed_EDTA
for file in `cat chrnames.txt`;
do
  gff2bed < ${file}.fwd.fa.mod.EDTA.TEanno.gff3 > completed_EDTA/${file}.fwd.fa.mod.EDTA.TEanno.bed
  tar -cvf ${file}.fwd.fa.mod.EDTA.final.tar ${file}.fwd.fa.mod.EDTA.final/ > ${file}.fwd.fa.mod.EDTA.final_zipped_files.txt
  rm -r ${file}.fwd.fa.mod.EDTA.final/
done
#Concatenate all chromosomes TE results
cat completed_EDTA/*.mod.EDTA.TEanno.bed >> All_chr_EDTA.TEanno.bed
#If it requires renaming of chromosomes
cat All_chr_EDTA.TEanno.bed | perl /PATH/TO/replace_contig_names.pl *_chrnames_to_replace.txt > All_chr_names_fixed_EDTA.TEanno.bed
