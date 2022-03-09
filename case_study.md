# Case study: Genome comparison between almond (Prunus dulcis) and apricot (Prunus armeniaca)
<img width="534" alt="Screen Shot 2022-03-06 at 7 06 24 PM" src="https://user-images.githubusercontent.com/91504464/156960974-db044e5d-5b7d-40a0-bb57-5745855a5566.png">
As it was mentioned in the introduction, this pipeline is not limited to one study system. For example, you could apply this pipeline to your animal genomes, microorganisms, insects, or plants as long as high-quality genome is available for both species you’re comparing, and they are reasonably close in their phylogeny. Here, we will walk through the pipeline using Prunus species in plant kingdom to illustrate what each command/code does and outputs, and let you adjust to your needs. 

Data collected from: 
- Alioto T, Alexiou KG, Bardil A, Barteri F, Castanera R, Cruz F, Dhingra A, Duval H, Fernández i Martí A, Frias L, Galán B, García JL, Howad W, Gómez‐Garrido J, Gut M, Julca I, Morata J, Puigdomènech P, Ribeca P, Rubio Cabetas MJ, Vlasova A, Wirthensohn M, Garcia‐Mas J, Gabaldón T, Casacuberta JM, Arús P. Transposons played a major role in the diversification between the closely related almond and peach genomes: results from the almond genome sequence. plant journal. 2020; 101(2):455-472.
- Groppi A, Liu S, Cornille A, Decroocq S, Bui QT, Tricon D, Cruaud C, Arribat S, Belser C, Marande W, Salse J, Huneau C, Rodde N, Rhalloussi W, Cauet S, Istace B, Denis E, Carrère S, Audergon JM, Roch G, Lambert P, Zhebentyayeva T, Liu WS, Bouchez O, Lopez-Roques C, Serre RF, Debuchy R, Tran J, Wincker P, Chen X, Pétriacq P, Barre A, Nikolski M, Aury JM, Abbott AG, Giraud T, Decroocq V. Population genomics of apricots unravels domestication history and adaptive events.. Nature communications. 2021 06 25; 12(1):3956.

## 0. Download genomes and gene annotation file and unzip all
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
## 1. Preprocess input data
### 1.1 Retain only chromosomes and remove the scaffolds in genome files
We use samtools faidx tool to create an index file listing names of chromosomes and extract them from the raw assembly fasta file. 
Depending on the chromosome names, fix the ">xxxx" name. 
```
module load samtools
cat refgenome.fa | grep ">Pd" | sed s/">"//g > refgenome_chr_names.txt
samtools faidx refgenome.fa -r refgenome_chr_names.txt > refgenome_chr.fa
cat qrygenome.fa | grep ">chr" | sed s/">"//g > qrygenome_chr_names.txt
samtools faidx qrygenome.fa -r qrygenome_chr_names.txt > qrygenome_chr.fa
```
### 1.2 Process gene annotation file
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

## 2. Genome alignment and structural variant detection 

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

Additionally, to compute species divergence score more accurately, we need to consider the length and % identity of each small aligned blocks that constitute those larger syntenic and invereted regions. This score is not found directly in any of the SyRI output files and we need to extract these using the command 02_genome_alignment_SV_detection/samtocoords.py or 02_genome_alignment_SV_detection/samtocoords.sh using minimap2 result (SAM). 
```
python3       #the three ">" appears which tells you that you're writing python command
>>> from syri.pyxFiles.synsearchFunctions import samtocoords
>>> table=samtocoords("/PATH/TO/minimap2_ref_qrygenome_chr_alignment.sam")
>>> table.to_csv(r'table.txt', index = False, sep = '\t') #makes table.txt file in CSV format in current directory
exit()
```
So now, you should have the following output files in your working directory: 
```
minimap2_ref_qrygenome_chr_alignment.sam
syri.out
syri.vcf
syri.summary 
invOut_table.txt 
synOut_table.txt
table.txt
```

At this point, it is useful to visually see what the alignment and structural variants look like between the two genomes you've compared. 

Once we have "table.txt", we can plot the alignment using 02_genome_alignment_SV_detection/visualize_genome_alignment.R which creates the following figure in R: 
![coords_file_for_visualization_original v0](https://user-images.githubusercontent.com/91504464/157148323-ef152261-3dc7-4be3-9236-dbd4bb84db3c.png)
Each box is representative of the aligned chromosome from refgenome on x-axis and qrygenome on y-axis. When the two genome sequences match perfectly (syntenic), then the straight right-side-up diagonal line will appear. What it means is that the position of the reference chromosome aligns with the same position on the query chromosome. So when large structural variants such as inversions are found, the position on query chromosome slides by some sequences/position, and aligns to the reference chromosome in the opposite direction resulting in left-side-up diagonal. Different colours represent the direction of the alignment (+ in blue or - in black based on which way the alignd sequence was read). 

For visualizing the structural variants identified by SyRI, fortunately, SyRI provides this straightforward script to produce the figure: 
```
source /PATH/TO/python_env
pip install matplotlib
python3 $PATH_TO_PLOTSR syri.out refgenome_chr.fa qrygenome_chr.fa -H 8 -W 5
```
The script makes "syri.pdf" figure in your curret directory which looks like this: 
<img width="578" alt="rotated" src="https://user-images.githubusercontent.com/91504464/157101596-a56ae954-b854-4a1c-9953-8eb49f2e7f84.png">

and you can see interesting genomic structural variants present between apricot and almond genomes more intuitively. It is always a good idea to double-check that both figures are roughly matching in alignment. Here, especially the large inversion in the end of chromosome 4 (Pd04) stands out in both plots. 

### Troubleshooting tips ###
If you are stuck on somewhere in the pipline and doesn't proceed to completion, it is a good idea to see which step is causing the problem. Break-up scripts are uploaded for you to do this. 
### 2.1 Genome alignment with minimap2
You can perform only the genome alignment step using minimap by submitting the following job: 
```
sbatch run_minimap2.sh #change the input refgenome and qrygenome as needed. 
```
This program is fairly computaitonally expensive, requiring lots of memory and CPU time to complete a complicated genome comparison. It may help if you run minimap2 separately from SyRI. 
### 2.2 SyRI in two rounds 
This may not be obvious from this example Prunus genomes but SyRI crashes and cannot finish the job if your chromosomes are not in the correct strand. DNA is double-stranded, and thus it has plus and minus strand. The program runs properly only if there is a sufficient amount of syntenic regions present between the two genomes. When one or more chromosomes of your query genome are somehow sequences of the opposite strands, then it fails to finish the job. This is why in some cases it requires you to run SyRI twice. 
First with the non-complementary raw genome and then identify the problematic strands: 
```
sbatch run_syri_first_round_only.sh {your Genus name}
```
If SyRI fails, the text file called "chr_to_rev.txt" will have a list of chromosome names that need flipping. Using this information, you proceed with reverse-complementing chromosomes by samtools, then rerun SyRI. 
```
sbatch run_syri_second_round_only.sh {your Genus name}
```
There are occasionally cases where this does not solve the issue and requires you to manually determine which chromosome needs reverse-complementing. For instance, when a genome consists of many small inversions, this can mess up the program to determine the overall pattern of syntenic alignment. The solution to this issue is to visualize the SAM file as mentioned above, then update the "chr_to_rev.txt" file, and rerun second round of SyRI with the same command as above. 


## 3. Species divergence analysis 
Now that we have the results of structural variants and the individual % identity scores from inversion and syntenic regions, we will calculate the divergence score between the two genomes of interest. 
We define: 
large blocks of inversion (inv) or syntenic region (syn) as inv region or syn region, respectively. 
small blocks of inversion (inv) or syntenic region (syn) as inv block or syn block, respectively. 
The formula to calculate the value termed "region_percent_identity" is below: 

```
- We know the score (= % identity) and length of individual syn/inv block 
- We want to know the score of each syn/inv region

For region A consisting of block a and b, 
 ---------   -------------------
    a                  b
|-------------------------------| region A 
where a = 10bp, 90% 
      b = 90bp, 70% 

Then, region A_percent_identity = length-normalized % identity for region A = {(10x0.9) + (90x0.7)} / (10+90) = 0.72
      region_percent_identity = length-normalized % identity for region x = {(length a x %a) + (length b x %b) + ...} / (length a + length b + ...) 
```

We do this calculation for every syn/inv region identified. Then we create a histogram of this score distribtion for syntenic region and inversion region using ggplot in R. The codes are found in 03_divergence_analysis/syriout_divergence_calc_plot.R and the output figures should look like these: 

syn_plot 
![Prunus_syn_plot](https://user-images.githubusercontent.com/91504464/157343849-21feb567-3232-4f2a-a0a5-b198fc77aaa5.png)

You should expect a normal distribution but some may appear skewed. The vertical line on syn_plot is the mean of the length-normalized region % identity for this dataset corresponding to species divergence score (100% = 100% identical region in sequence, same species/genomes). Because random mutations occur and fix over time in the evolutionary timescale, the less % identity in sequence indicates the more divergent species relationship. Thus, we can interpret the x-axis as genomic divergence between the two species. 

inv_plot
![Prunus_inv_plot2](https://user-images.githubusercontent.com/91504464/157343746-76a3c79f-1e2a-48a9-bf6a-34e85c213e9b.png)
For inversions, the idea of sequence identity is similar to how we would interpret the syn_plot. Although inversions do reduce recombination rate within the region compared to other parts of the genome, they can still undergo mutations over time. As a result, the lower % identity or the smaller value in x-axis on this plot indicate the presence of potentially old, ancient inversions occurred in the ancestor of the two species. 
!Beware of the big difference in counts as there are much fewer inversions than syntenic regions between closely related species genomes. 
