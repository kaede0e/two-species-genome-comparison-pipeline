# Comparative genomics: Two closely related species genome alignment, structural variant detection & TE annotation pipeline 
Two-species genome comparison pipeline (combines minimap2/SyRI, EDTA, bedtools) to identify structural variants (inv/syn), annotate coding sequence and TE, and plot results in R. 

### Genome alignment and variant detection
-	Genome alignment is useful because it tells you how the two species are related at genomic level. 
-	Relevant for phylogenetics, evolutionary analysis, population genetics etc. 
-	Genomic variants include large-scale structural variations (inversions, translocations, etc.) that are not captured by genome-wide genetic polymorphism scans like SNPs and InDels identifications. Whole genome alignment has advantages to capture synonymous changes which are more common than nonsynonymous, interesting genomic variants that affect the genome evolution in a non-traditional central dogma (Transposable elements or selfish elements that don’t code for particular genes but affects the organism’s phenotype, etc.) 
-	Each step in this pipeline can be performed individually as long as you have the proper input data. But as a whole pipeline, it helps you to produce a readily useful statistics that indicates interesting information about your species of interest at genomic level. 
-	Widely applicable to any organisms you’re interested in. No upper or lower limit in terms of genome complexity (theoretically), though the more complicated the genome the more computational resource it requires to run. 


### Limitation of this pipeline
-	relies on the quality of genome you have. Only chromosome-resolved genomes can be used as the input. 
-	The two species of interest must contain sufficient genomic region that is syntenic 
-	requires appropriate amount of computational resource as shown: 

| Input genome (size)  | Step  | time required (hrs) | CPU memory (GB) |
| :------------------- |:-----:| :------------------:| :-------------: |
| 1 Prunus (202Mbp)    | minimap2 | 3               | 3              |
|                      | SyRI     | 1               | 3              |
|                      | EDTA     | 36*             | 5.5            |
| 2 Arachis (1.01Gbp)  | minimap2+SyRI|    8.7      | 72             |
|                      | EDTA     | 142*             | 18.9          |

*Wall-clock time using 16 CPU, run by chromosome-by-chromosome. 

## General pipeline structure: 
![Picture1](https://user-images.githubusercontent.com/91504464/156959886-30dc48cc-4eba-41bd-8e9d-13a731b5e7d6.png)

## Dependencies/Installation
The pipeline requires following software packages to run: 
- samtools http://www.htslib.org/
- minimap2 https://github.com/lh3/minimap2
- Synteny and Rearrangement Identifier (SyRI) https://schneebergerlab.github.io/syri/
- Extensive de-novo TE Annotator (EDTA) https://github.com/oushujun/EDTA
- bedtools https://bedtools.readthedocs.io/en/latest/#

Other dependencies (in Compute Canada clusters): 
- python_env or virtual environment to install python packages https://docs.computecanada.ca/wiki/Python
- Singularity container https://docs.computecanada.ca/wiki/Singularity/en

## Run the pipeline
0. Collect genomes (fasta) of two species (referred to as reference and query genome = refgenome and qrygenome) and gene annotation file (gff3 or bed) for reference genome
1. **01_preprocess_input**: Pre-process genomes to remove scaffolds from genome assembly and convert gene annotation file into bed format containing only coding sequence (CDS)
2. **02_genome_alignment_SV_detection**: Aligns two refgenome and qrygenome by minimap2 and performs structural variant detection using SyRI. Then the outputs are visualized in R Studio. 
3. 
