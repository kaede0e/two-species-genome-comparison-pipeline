# Case study: Genome comparison between almond (Prunus dulcis) and apricot (Prunus armeniaca)
<img width="534" alt="Screen Shot 2022-03-06 at 7 06 24 PM" src="https://user-images.githubusercontent.com/91504464/156960974-db044e5d-5b7d-40a0-bb57-5745855a5566.png">
As it was mentioned in the introduction, this pipeline is not limited to one study system. For example, you could apply this pipeline to your animal genomes, microorganisms, insects, or plants as long as high-quality genome is available for both species you’re comparing, and they are reasonably close in their phylogeny. Here, we will walk through the pipeline using Prunus species in plant kingdom to illustrate what each command/code does and outputs, and discuss what they mean or might indicate in terms of molecular evolution. 


Data for this tutorial collected from: 
- Alioto T, Alexiou KG, Bardil A, Barteri F, Castanera R, Cruz F, Dhingra A, Duval H, Fernández i Martí A, Frias L, Galán B, García JL, Howad W, Gómez‐Garrido J, Gut M, Julca I, Morata J, Puigdomènech P, Ribeca P, Rubio Cabetas MJ, Vlasova A, Wirthensohn M, Garcia‐Mas J, Gabaldón T, Casacuberta JM, Arús P. Transposons played a major role in the diversification between the closely related almond and peach genomes: results from the almond genome sequence. plant journal. 2020; 101(2):455-472.
- Groppi A, Liu S, Cornille A, Decroocq S, Bui QT, Tricon D, Cruaud C, Arribat S, Belser C, Marande W, Salse J, Huneau C, Rodde N, Rhalloussi W, Cauet S, Istace B, Denis E, Carrère S, Audergon JM, Roch G, Lambert P, Zhebentyayeva T, Liu WS, Bouchez O, Lopez-Roques C, Serre RF, Debuchy R, Tran J, Wincker P, Chen X, Pétriacq P, Barre A, Nikolski M, Aury JM, Abbott AG, Giraud T, Decroocq V. Population genomics of apricots unravels domestication history and adaptive events. Nature communications. 2021 06 25; 12(1):3956.

## 0. Download genomes and gene annotation file and unzip all
The default scripts are written in a way so that all your genomes and other data are deposited under a parent directory named "xxxx_genome" where xxxx = Genus name. 
In this tutorial, we are analyzing the Prunus genome so we are going to make a directory named "Prunus_genome". 
```
mkdir Prunus_genome
cd Prunus_genome
```
Then we will download the following genome assembly into the directory we just created. 

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
The gff3 is the common format for gene annotation. Annotation file describes what feature of the genome is found in where in the genome. It generally contains many information including the locations of coding sequence (CDS), genes, exons, mRNA, etc. in the genome. However, our downstream analysis requires the use of bed formatted files (another way of describing locations of features; column 1 - chromosome name, column 2 - start of chromosome, column 3 - end of chromosome) and we're only interested in CDS. The reason for choosing only CDS is that we want to look at regions that could make a significant phenotypic effect to the organism. CDS is where the sequences are actively translated to proteins, more likely involved with functional roles in an organism and it's where the nonsynonymous mutations arise that lead to evolutionary implications. 

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

This is a multi-step section as outlined in SyRI pipeline for two species genomes.md. The script to run this section in one-step is found in 02_genome_alignment_SV_detection/run_syri_automation.sh. If you are to reproduce the whole pipeline, you would simply run: 
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
"minimap2_ref_qrygenome_chr_alignment.sam" is your alignment result in SAM format, "syri.out" and "syri.vcf" are the genomic structural differences identified by SyRI in tabular format (TSV and VCF format respectively) which includes syntenic region, inverted region, translocation, duplication, and small variants (SNPs and InDels). "syri.summary" summarises the stats related to how many of those features were found in the genome alignment. "invOut_table.txt" and "synOut_table.txt" provide us the useful % identity scores for each aligned region. Because we are interested in comparing syntenic region vs inverted region and computing species divergence score from sequence alignment, these score tables are crucial in the subsequent analyses. 

Additionally, to compute species divergence score more accurately, we need to consider the length and % identity of each small aligned blocks that constitute those larger syntenic and invereted regions. This score is not found directly in any of the SyRI output files and we need to extract these using the command 02_genome_alignment_SV_detection/samtocoords.py or 02_genome_alignment_SV_detection/samtocoords.sh using minimap2 result (SAM). On my Linux command line system, it looks like this: 
```
python3       #the three ">" appears which tells you that you're writing python command
>>> from syri.pyxFiles.synsearchFunctions import samtocoords
>>> table=samtocoords("/PATH/TO/minimap2_ref_qrygenome_chr_alignment.sam")
>>> table.to_csv(r'table.txt', index = False, sep = '\t') #makes table.txt file in CSV format in current directory
exit()
```
Once you get the score table, you should have the following output files in your working directory: 
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
Overall, this plot shows a good genome alignment between _Prunus dulcis_ and _Prunus armeniaca_ which assures that the program ran properly. However, it's worth mentioning that some sections on chromosome 7 and 8 have poor alignment (where the two parallel diagonal lines appear), which may be indicative of partial genome duplications or genome misassemblies. As it was mentioned in the introduction of the module, the pipeline is limited by the quality of genome assemblies. The more accurate assemblies are, the more accurate our interpretation of the results. 

For visualizing the structural variants identified by SyRI, you can use their built-in script to produce the figure: 
```
source /PATH/TO/python_env
pip install matplotlib
python3 $PATH_TO_PLOTSR syri.out refgenome_chr.fa qrygenome_chr.fa -H 8 -W 5
```
The script above makes "syri.pdf" figure in your curret directory which looks like this: 
<img width="578" alt="rotated" src="https://user-images.githubusercontent.com/91504464/157101596-a56ae954-b854-4a1c-9953-8eb49f2e7f84.png">

Now you can see interesting genomic structural variants present between the two genomes more intuitively. It is always a good idea to double-check that both figures are roughly matching in alignment. Here, especially the large inversion in the end of chromosome 4 (Pd04) stands out in both plots. Rich source of small-size variations is detected in chromosome 2 and 3 (Pd02, Pd03) while chromosome 1 and 5 (Pd01, Pd05) are highly conserved. These are qualitative observations that can serve as a good starting point to dig further into what evolutionary importance these divergent regions may have to the species. 

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

Synteny plpt
![Prunus_syn_plot](https://user-images.githubusercontent.com/91504464/157343849-21feb567-3232-4f2a-a0a5-b198fc77aaa5.png)

You should expect a normal distribution but some may appear skewed. The vertical line on syn_plot is the mean of the length-normalized region % identity for this dataset corresponding to species divergence score (100% = 100% identical region in sequence, same species/genomes). Because random mutations occur and fix over time in the evolutionary timescale, the less % identity in sequence indicates the more divergent species relationship. Thus, we can interpret the x-axis as genomic divergence between the two species. 

Inversion plot
![Prunus_inv_plot2](https://user-images.githubusercontent.com/91504464/157343746-76a3c79f-1e2a-48a9-bf6a-34e85c213e9b.png)

For inversions, the idea of sequence identity is similar to how we would interpret the syn_plot. Although inversions do reduce recombination rate within the region compared to other parts of the genome, they can still undergo mutations over time with baseline recombinations. As a result, the lower % identity or the smaller value in x-axis on this plot indicate the presence of potentially old, ancient inversions occurred in the ancestor of the two species which could be interesting to look into. These old inversions mean that they got fixed and maintained in the population, could have had neutral impact to the species or beneficial, but most likely not deleterious. (Beware of the big difference in counts as there are much fewer inversions than syntenic regions between closely related species genomes. Don't overinterpret things based on the instant appearances.) 

### Checkpoint ### 
Taking these results up to here, you can plan a further investigation into certain structural variants detected, or maybe do a BLAST search on what genes are found in these regions that might have an evolutionary significance. Performing the same steps on multiple species of the same genus might also give you an interesting insight about species relationships and genomic architecture for your study organisms. 


If you are keen to investigate the difference in some genetic content between syntenic region vs. inverted region (coding sequence (CDS) content or trnasposable element (TE)), then you may wish to continue with the pipeline below. 


## 4. _de novo_ Transposable Element annotation pipeline 
The next few sections are driven by some hypotheses around genomic inversions.  
1. Like most mutations are deleterious, and disrupting genes by inversions is bad. Therefore, we would expect fewer genes in inversions than syntenic regions. 
2. Like most mutations are deleterious, and breaking up genes by inversions is bad. Therefore, we would expect fewer genes in inversion breakpoints compared to the rest of the genome. 
3. Transposable elements (TEs) are thought to be the driver of inversions. Therefore, we would expect more TEs in inversions than syntenic regions. 
4. Transposable elements (TEs) are thought to be the driver of inversions. Therefore, we would expect enriched presence of TEs in inversion breakpoints compared to the rest of the genome. 

To test the above hypotheses, first, we need to obtain TE library for the refgenome. 
We will use the extensive de novo TE annotator (EDTA) pipeline which works by making educational guesses on sequences that have structural characteristics that define TE classes and superfamilies. The advantages of de novo method include: it does not rely on high copy number assumption, does not rely on sequence homology from reference (model organism) genome sequences, prior knowledge, or databases. 

The program EDTA depends on many third-party software packages as it is a pipeline that combines multiple separate TE finders (LTR-FINDER, LTRharvest, generic repeat finder, HelitronScanner, etc. - refer to their Github page for details) and then filters them for redundancy, and additionally perform de novo TE annotator (RepeatModeler) to produce the curated nonredundant TE library. As a result, the pipeline comes in a package that contains everything you need which use you can download and access via either Conda or Singularity. 
Here, we will use Singularity module to run the pipeline. 
```
singularity shell -B /home -B /PATH/TO/EDTA/EDTA-2.0.0.sif \
perl /PATH/TO/EDTA/EDTA.pl \
--genome refgenome_chr.fa \
--cds genome_cds.fasta \
--exclude genome_cds.bed \
--overwrite 1 --sensitive 1 --anno 1 --threads ${SLURM_CPUS_PER_TASK}
```
This command starts an interactive job that goes all the steps in EDTA pipeline and runs until completion. However, as it is a resource-heavy program and could run for days and weeks depending on the size of your genome, it may be more feasible to use scripts like 04_de_novo_TE_annotation/example_run_EDTA_by_chr.sh which divides the genome into chromosomes and performs EDTA pipeline separately. The biggest limitation for doing this is potentially missing TEs as it cannot look for TEs found across chromosomes, but it saves a lot of computational resource (and your time to wait for the job to complete) this way. 

The program produces many output files, but in the end we curate the results from separate chromosomes into a combined file named 
```
All_chr_EDTA.TEanno.bed
```
which is the final list of TEs annotated to your refgenome with information about what type of TEs are found for each feature in bed format. This will be used in the next step in parallel with the CDS data to compare their content in syntenic region vs inverted region. 

## 5. Coding sequence and transposable element in syntenic region vs. inversion 
For figuring out how many CDS and TEs are found within syn/inv regions, we will use the software called bedtools. To run the program, you run the script found in 05_CDS_TE_overlap_analysis/bedtools_for_all.sh. The script is written in a way that lets you perform the overlap analyses across multiple genera listed in a text file called "genome_pair.txt" at once if you desire. So, the first thing you need to do is to create a list of genera you're running the analysis for. In this example module, our analysis has only one genus under the folder named "Prunus_genome" so we create a text file "genome_pair.txt" that looks like: 
```
Prunus
```

Then you call the script: 
```
bash bedtools_for_all.sh
```

The script performs the basepair overlap analysis between the two bed files you input. In this case, we are detecting the overlap between syn/inv region and CDS/TE features respectively. At the end, you should have these eight result files created: 
```
bedtools_count_cds_inv.txt
bedtools_bpoverlap_merged_cds_inv.txt
bedtools_count_cds_syn.txt
bedtools_bpoverlap_merged_cds_syn.txt
bedtools_count_TE_inv.txt
bedtools_bpoverlap_merged_TE_inv.txt
bedtools_count_TE_syn.txt
bedtools_bpoverlap_merged_TE_syn.txt
```
They are all in tabular format. Note that the basepair overlap analyses was performed on the 'merged' CDS/TE file which means any overlapping features (eg. some gene annotation file includes essentially the same gene multiple times but identifies them as isoforms, which can lead to falsely overrepresentated regions) are combined and treated as one big section of features. This conveys more accurate information we're looking for (ie. how much proportion of inv/syn region contains CDS/TE). 

The count files give you how many CDS/TE are found within the syn/inv regions on the last column. 
The bpoverlap files give you how many basepairs of CDS/TE are found to overlap with syn/inv regions on the last column. 

Using these data, we can visualize and analyze trends in inversion vs syntenic region with respect to CDS and TE content. For instance, we can ask questions like 
"Are there more CDS in the syntenic region than inverted regions?" 
"What type of TEs are most prevaent in the genome? Are there differences in syn/inv regions?" 
"Does the trend you see differ between different genera?" 

### Some ideas for further analysis in R 
I have curated some example data analysis and plots you can make from bedtools analysis in 05_CDS_TE_overlap_analysis/example_Prunus_data_analysis.R

1) How many CDS/TE are there? Is there a difference in syntenic region vs inversion? 
```
# A tibble: 2 × 3
  type  total_number_of_CDS total_number_of_TE
  <chr>               <int>              <int>
1 inv                  4838               5482
2 syn                142922              51133
```
If you consider the ratio between CDS and TE (remember, there are much fewer inv than syn regions analyzed), this roughly supports our first hypothesis that more TEs are found in inversions than in syntenic regions. In fact, it is assuring that there are so much more CDS found in syntenic regions (x3 the number of TEs) compared to inversions where there are actually more TEs than CDS! 

2) How much proportion of inv/syn region contains CDS or TE sequences? 

CDS: 
<img width="384" alt="Screen Shot 2022-03-11 at 10 19 09 AM" src="https://user-images.githubusercontent.com/91504464/157926987-25ca0a42-25c5-4e1d-962a-5aa5f16a0fbe.png">

This violin plot shows the distribution of CDS/TE proportions per inv/syn region. The area of the violin corresponds to the number of data points. This was unexpected and rather a weird looking distribution as we expected far higher CDS proportion in syntenic region than inversions while the results are showing more proportion of CDS in inversions than syntenic regions. I would suspect that this is due to very small inversion regions identified compared to syntenic regions, making the comparison rather unequally weighted. There is no further explanations I can think of from this plot, unfortunately. 


Okay, what about TEs? 

TE: 
<img width="380" alt="Screen Shot 2022-03-11 at 10 17 30 AM" src="https://user-images.githubusercontent.com/91504464/157926482-40f6e2a9-6b86-4604-a3bc-15f45e0cf469.png">

The same violin plot for TE seems to provide no support for our TE driving inversion hypothesis. It's good to keep in mind whenever you encounter a non-significant result, remember the null hypothesis and importance of null result (neutral theory). 

But this is not fun to end this tutorial. Let's see if we can analyze something interesting about TEs in this genome. Recall from EDTA pipeline, it made a lot of output files and one of them includes stats summary file called "xxxx.mod.EDTA.TEanno.sum" where xxxx is your input genome/chromosome (if you did by chromosome-by-chromosome). The file shows you what type of TEs were found and how much they occupy relative to the whole input genome/chromosome length. 

This summary stats from EDTA tell us: 
```
# A tibble: 8 × 2
  chr   total_percent_of_TE
  <chr>               <dbl>
1 Pd01                 22.1
2 Pd02                 22.5
3 Pd03                 21.8
4 Pd04                 22.0
5 Pd05                 16.6
6 Pd06                 22.8
7 Pd07                 19.3
8 Pd08                 18.5
```
Which seems to suggest about 20% of the almond genome is composed of transposable elements. Compared to other plant species or animal genomes, this seems a bit low. However, as we discussed in 'eukaryotic genome' lecture, there is a strong positive correlation between genome size and TE content in eukaryotic genomes. Taking into consideration that _Prunus dulcis_ genome is relatively small (228Mbp), this TE content can be said average. As a comparison (https://doi.org/10.1155/2011/893546), Arabidopsis genome is ~15% TEs (125Mbp), rice genome is ~35% TEs (472Mbp), maize genome is ~85% TEs (2.3Gbp). TEs have recently gained a lot of attention and have shifted the concept of being 'junk' DNA to a key regulator that holds evolutionary roles especially in plant genomes. Plants are sessile. It is stuck in where it grows. But at the same time, it needs to fight against both biotic and abiotic stresses. At times of changing climate, for instance, it cannot escape by moving to a different habitat but it must adapt by changing its own abiotic stress-coping mechanisms. TEs can be their great molecular tools because of its 'mobile' nature which allows to proliferate in the genome, allowing a rapid physiological change which may be more suited to the changing climate. 

Additionally, I'm interested in asking: What type of TEs are the most prevalent in _Prunus dulcis_? 
Because we've run the pipeline chromosome-by-chromosome, we can additionally ask: "Does the trend differ by chromosomes?"
To answer these questions, we can do some stats in R to get summary like this based on the TE category available in EDTA: 
```
Prunus_TE_table %>%
group_by(Class)%>%
summarise(mean_percent = mean(percent_masked))

# A tibble: 11 × 2
   Class               mean_percent
   <chr>                      <dbl>
 1 LTR_Copia                 3.18  
 2 LTR_Gypsy                 2.35  
 3 LTR_unknown               3.31  
 4 nonLTR_LINE_element       0.0512
 5 nonTIR_helitron           1.12  
 6 repeat_region             4.07  
 7 TIR_CACTA                 1.21  
 8 TIR_hAT                   1.16  
 9 TIR_mutator               2.87  
10 TIR_PIF_Harbinger         1.37  
11 TIR_Tc1_Mariner           0.0312
```
The "repeat_region" is not very informative as it doesn't fall under any families of TEs, but we can infer that LTR or long terminal repeat is the most common class as consistent with other plant genomes. In _Prunus dulcis_ genome, the TIR or terminal inverted repeat seems also common. 


I hope the pipeline gave you a hands-on experience to perform genomic analysis and reveal something interesting about two species genomes of your interest! 
