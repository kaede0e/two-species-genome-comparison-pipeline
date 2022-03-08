# Case study: Genome comparison between almond (Prunus dulcis) and apricot (Prunus armeniaca)
<img width="534" alt="Screen Shot 2022-03-06 at 7 06 24 PM" src="https://user-images.githubusercontent.com/91504464/156960974-db044e5d-5b7d-40a0-bb57-5745855a5566.png">
As it was mentioned in the introduction, this pipeline is not limited to one study system. For example, you could apply this pipeline to your animal genomes, microorganisms, insects, or plants as long as high-quality genome is available for both species you’re comparing, and they are reasonably close in their phylogeny. Here, we will walk through the pipeline using Prunus species in plant kingdom to illustrate what each command/code does and outputs, and let you adjust to your needs. 

Data collected from: 
- Alioto T, Alexiou KG, Bardil A, Barteri F, Castanera R, Cruz F, Dhingra A, Duval H, Fernández i Martí A, Frias L, Galán B, García JL, Howad W, Gómez‐Garrido J, Gut M, Julca I, Morata J, Puigdomènech P, Ribeca P, Rubio Cabetas MJ, Vlasova A, Wirthensohn M, Garcia‐Mas J, Gabaldón T, Casacuberta JM, Arús P. Transposons played a major role in the diversification between the closely related almond and peach genomes: results from the almond genome sequence. plant journal. 2020; 101(2):455-472.
- Groppi A, Liu S, Cornille A, Decroocq S, Bui QT, Tricon D, Cruaud C, Arribat S, Belser C, Marande W, Salse J, Huneau C, Rodde N, Rhalloussi W, Cauet S, Istace B, Denis E, Carrère S, Audergon JM, Roch G, Lambert P, Zhebentyayeva T, Liu WS, Bouchez O, Lopez-Roques C, Serre RF, Debuchy R, Tran J, Wincker P, Chen X, Pétriacq P, Barre A, Nikolski M, Aury JM, Abbott AG, Giraud T, Decroocq V. Population genomics of apricots unravels domestication history and adaptive events.. Nature communications. 2021 06 25; 12(1):3956.

## 0. Preprocess input data

### 0.1 Download genomes and gene annotation file and unzip all. 
The default scripts are written in a way so that all your genomes and other data are deposited under a parent directory named "Genus_genome". 
In this tutorial, we are analyzing the Prunus genome so we are going to make a directory named "Prunus_genome". 
```
mkdir Prunus_genome
cd Prunus_genome
```
Then we will download the following raw genome assembly into the directory we just created. 

Prunus dulcis genome (refgenome): 
```
wget https://www.rosaceae.org/rosaceae_downloads/Prunus_dulcis/pdulcis_v2.0/assembly/pdulcis26.chromosomes.fasta.gz
wget https://www.rosaceae.org/rosaceae_downloads/Prunus_dulcis/pdulcis_v2.0/genes/Prudul26A.cds.fa.gz
wget https://www.rosaceae.org/rosaceae_downloads/Prunus_dulcis/pdulcis_v2.0/genes/Prudul26A.chromosomes.gff3.gz
gunzip *
```
Rename to generalize scripts
```
mv pdulcis26.chromosomes.fasta refgenome.fa
mv Prudul26A.cds.fa refgenome_cds.fa
mv Prudul26A.chromosomes.gff3 refgenome.gff3
```
Prunus armeniaca genome (qrygenome): 
```
wget https://www.rosaceae.org/rosaceae_downloads/Prunus_armeniaca/parmeniaca_cv.Stella_v1/assembly/Prunus_armeniaca_cv_Stella.fasta.gz
gunzip Prunus_armeniaca_cv_Stella.fasta.gz
```
Rename to generalize scripts
```
mv Prunus_armeniaca_cv_Stella.fasta qrygenome.fa
```
### 0.2 Retain only chromosomes and remove the scaffolds in genome files. 
We use samtools faidx tool to create an index file listing names of chromosomes and extract them from the raw assembly fasta file. 
Depending on the chromosome names, fix the ">xxxx" name. 
```
module load samtools
cat refgenome.fa | grep ">Pd" | sed s/">"//g > refgenome_chr_names.txt
samtools faidx refgenome.fa -r refgenome_chr_names.txt > refgenome_chr.fa
cat qrygenome.fa | grep ">chr" | sed s/">"//g > qrygenome_chr_names.txt
samtools faidx qrygenome.fa -r qrygenome_chr_names.txt > qrygenome_chr.fa
```
### 0.3 Process gene annotation file. 
The gff3 is the common format for gene annotation. Annotation file describes what feature of the genome is found in where in the genome. It generally contains many information including the locations of coding sequence (CDS), genes, exons, mRNA, etc. in the genome. However, our downstream analysis requires the use of bed formatted files (another way of describing locations of features) and we're only interested in CDS. The reason for choosing only CDS is that we want to look at regions that could make a significant phenotypic effect to the organism. CDS is where the sequences are actively translated to proteins, performing functional roles in an organism and it's where the nonsynonymous mutations arise that lead to evolutionary implications. 

Extract CDS from gff3 (or any gff formatted gene annotation file)
```
awk '$3 ~ /CDS/' refgenome.gff3 > genome_cds.gff
```
Convert gff3 into bed format
```
module load bedops/2.4.39
gff2bed < genome_cds.gff > genome_cds.bed
```
To fit with the scripts, you'd need to 
Now you have all the input data you need for the pipeline to run! 

## 1. Genome alignment and structural variant detection 

The script to run this section in one step is found in 02_genome_alignment_SV_detection/run_syri_automation.sh. If you are to reproduce the whole pipeline, you would simply run: 
```
cd ..
sbatch run_syri_automation.sh Prunus
```
which will submit the job script that is coded to run minimap2 and subsequent structural variant detection software (SyRI) in "Prunus_genome" directory. 

When the job is done successfully, following outputs are produced in "Prunus_genome" folder: 
```
minimap2_ref_qrygenome_chr_alignment.sam
syri.out
syri.vcf
syri.summary 
invOut_table.txt 
synOut_table.txt
```
"minimap2_ref_qrygenome_chr_alignment.sam" is your alignment result in SAM format, "syri.out", "syri.vcf" are the genomic structural differences identified by SyRI in tabular format (TSV and VCF format respectively) which includes syntenic region, inverted region, translocation, duplication, and small variants (SNPs and InDels). "syri.summary" summarises the stats related to how many of those features were found in the genome alignment. "invOut_table.txt" and "synOut_table.txt" provide us the useful % identity scores for each aligned region. Because we are interested in comparing syntenic region vs inverted region and computing species divergence score from sequence alignment, these score tables are crucial in the subsequent analyses. 

Additionally, for computing species divergence score more accurately, we need to consider the length and % identity of each small aligned blocks that constitute those larger syntenic and invereted regions. This score is not found directly in any of the SyRI output files and we need to extract these using the command 02_genome_alignment_SV_detection/samtocoords.py or 02_genome_alignment_SV_detection/samtocoords.sh using minimap2 result (SAM). 

At this point, it is useful to visually see what the alignment and structural variants look like between the two genomes you've compared. Once we have "table_original.txt", we can plot the alignment using 02_genome_alignment_SV_detection/visualize_genome_alignment.R which creates the following figure in R: 
![coords_file_for_visualization_original v0](https://user-images.githubusercontent.com/91504464/157148323-ef152261-3dc7-4be3-9236-dbd4bb84db3c.png)

Each box is representative of the aligned chromosome from refgenome on x-axis, qrygenome on y-axis. When the two genome sequences match perfectly (syntenic), then the straight right-side-up diagonal line will appear. What it means is that the one location of the refgenome chromosome aligns with the same location on the qrygneome chromosome. Different colours represent the direction of the alignment (+ or - based on which way the sequence should be read). 

For visualizing the structural variants identified by SyRI, fortunately, SyRI provides this straightforward script to produce the figure: 
```
source /PATH/TO/python_env
pip install matplotlib
python3 $PATH_TO_PLOTSR syri.out refgenome_chr.fa qrygenome_chr.fa -H 8 -W 5
```
The script produces "syri.pdf" figure in your curret directory which looks like this: 
<img width="578" alt="rotated" src="https://user-images.githubusercontent.com/91504464/157101596-a56ae954-b854-4a1c-9953-8eb49f2e7f84.png">
and you can see interesting genomic structural variants present between apricot and almond genomes. Especially the large inversion in the end of chromosome 4 (Pd04) stands out to me. 
This plot gives you a rough idea of how well the two genomes align at nucleotide sequence level. 


Another common way of visualizing alignment result that does not depend on SyRI software is to simply plot on R. Using the 02_genome_alignment_SV_detection/visualize_genome_alignment.R code in R studio, you can straightaway visualize the minimap2 result (SAM file) 

Genome alignment is performed by minimap2 to align qrygenome onto refgenome. It uses seed-chain-align procedure in which . 
To run the program, you submit the job script 
