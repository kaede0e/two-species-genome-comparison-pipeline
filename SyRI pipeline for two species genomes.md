## Detailed Outline for SyRI structural variant detection
Batch process for SyRI analysis
0) Give paths to the command
1) Load reference and query genome (fasta) = genome.fa
2) Process reference and query genome
2.1) Remove fragments and retain chromosomes only

3) Run SyRI round1
3.1) Perform whole genome alignment using minimap2
3.2) Run SyRI
 #If this command gives error reverse those chromosomes that don't match, proceed to #4.
      #look at the error message it prints, and figure out which ones to reverse.
4) Reverse complement strands that got inverted from mapping result *manual ID required if SyRI fails to identify the correct ones. 

5) Run SyRI round2
5.1) Perform whole genome alignment using minimap2 on reverse-complemented genome
5.2) Run SyRI
5.3) Visualize the alignment in PDF format (built-in SyRI function) if you wish. 

6) Output divergence score (coords file) from SyRI analysis (Python)

7) Format SyRI output into tidy table

8) Merge outputs and make tables in R Studio
8.1) Import data frame
8.2) Make data frame into tibble
8.3) Rename variables
8.4) Merge tables
8.5) Plot histogram
