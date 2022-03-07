# Case study: Genome comparison between almond (Prunus dulcis) and apricot (Prunus armeniaca)
<img width="534" alt="Screen Shot 2022-03-06 at 7 06 24 PM" src="https://user-images.githubusercontent.com/91504464/156960974-db044e5d-5b7d-40a0-bb57-5745855a5566.png">
As it was mentioned in the introduction, this pipeline is not limited to one study system. For example, you could apply this pipeline to your animal genomes, microorganisms, insects, or plants as long as high-quality genome is available for both species you’re comparing, and they are reasonably close in their phylogeny. Here, we will walk through the pipeline using Prunus species in plant kingdom to illustrate what each command/code does and outputs, and let you adjust to your needs. 

Data collected from: 
- Alioto T, Alexiou KG, Bardil A, Barteri F, Castanera R, Cruz F, Dhingra A, Duval H, Fernández i Martí A, Frias L, Galán B, García JL, Howad W, Gómez‐Garrido J, Gut M, Julca I, Morata J, Puigdomènech P, Ribeca P, Rubio Cabetas MJ, Vlasova A, Wirthensohn M, Garcia‐Mas J, Gabaldón T, Casacuberta JM, Arús P. Transposons played a major role in the diversification between the closely related almond and peach genomes: results from the almond genome sequence. plant journal. 2020; 101(2):455-472.
- Groppi A, Liu S, Cornille A, Decroocq S, Bui QT, Tricon D, Cruaud C, Arribat S, Belser C, Marande W, Salse J, Huneau C, Rodde N, Rhalloussi W, Cauet S, Istace B, Denis E, Carrère S, Audergon JM, Roch G, Lambert P, Zhebentyayeva T, Liu WS, Bouchez O, Lopez-Roques C, Serre RF, Debuchy R, Tran J, Wincker P, Chen X, Pétriacq P, Barre A, Nikolski M, Aury JM, Abbott AG, Giraud T, Decroocq V. Population genomics of apricots unravels domestication history and adaptive events.. Nature communications. 2021 06 25; 12(1):3956.

## 0. Preprocess input data
The following uses bash script. 
``` 
#!/bin/bash
```
### 0.1 Download genomes and gene annotation file and unzip all. 
Prunus dulcis genome (refgenome) 
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
Prunus armeniaca genome (qrygenome) 
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

Now you have all the input data you need for the pipeline to run! 

## 1. 
