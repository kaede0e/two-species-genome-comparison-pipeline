#Segmental duplications at downstream/upstream of inversion breakpoints
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

#Data Import function
## Define directory with your data formatted as directory/genera/
genera <- c("Acer", "Actinidia", "Arabis", "Arachis",
            "Citrus", "Corylus", "Corymbia", "Cucumis", 
            "Eukalyptus", "Fragaria", "Glycine", "Gossypium", 
            "Ipomoea", "Juglans", "Malus", "Medicago", 
            "Phaseolus", "Populus", "Prunus", "Pyrus", 
            "Quercus", "Rhododendron", "Rosa", "Rubus", 
            "Salix", "Salvia", "Solanum", 
            "Vaccinium", "Vigna", "Vitis") #26 species pair had blast hit
empty_blast_output <- c("Actinidia", "Cucumis", "Eukalyptus", "Rubus")

#Data Import function 
datalist = list()

#### Blast output #####
for (genus in genera){
  blast_output_table <- read.table(paste0(directory, "/", genus, "/refgenome_inv_blast_score_aln.txt"), head = FALSE) %>%
    rename(qseqID = V1, 
           refseqID = V2, 
           percent_ID = V3, 
           length= V4,
           mismatch = V5,
           gap_open = V6,
           qstart = V7,
           qend = V8, 
           refstart = V9,
           refend = V10,
           E_value = V11, 
           bitscore = V12) %>% 
    add_column(Genus = genus) %>% 
    separate(qseqID, into = c("chr_1", "coordinates"), sep = ":") %>%
    separate(coordinates, into = c("region_start_1", "region_end_1"), sep = "-", convert = TRUE)
  bedtools_bpoverlap_breakpoint_regions_bed <- read.table(paste0(directory, "/", genus, "/start_breakpoints_regions_onesided.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus) 
  blast_output_table_joined <- left_join(bedtools_bpoverlap_breakpoint_regions_bed, blast_output_table)
  datalist[[genus]] <- blast_output_table_joined # add it to your list
}

blast_output = do.call(rbind, datalist)

write_csv(blast_output, "Master_inv_breakpoints_aln_blast_output_10kb_window.csv")

##### Cleaning up dataset #####
blast_10kb_window <- read.csv("Master_inv_breakpoints_aln_blast_output_10kb_window.csv")
blast_10kb_window_table <- blast_10kb_window %>%
  filter(Genus != "Raphanus") %>% 
  filter(Genus != "Luffa") %>%
  mutate(region_length = region_end_1 - region_start_1) %>% 
  filter(region_length >= 9000)

blast_10kb_window_table %>%
  ggplot()+
  geom_histogram(aes(x=percent_ID))+
  xlab("sequence similarity of 10kb window (%)")+ ylab("count")
  
duplicates <- blast_10kb_window_table %>% group_by(chr_1, region_start_1, region_end_1, Genus) %>% count() %>% filter(n != 1)
duplicate_hits <- left_join(blast_10kb_window_table, duplicates)
blast_longest_hits <- duplicate_hits %>% filter(n >= 1) %>% group_by(chr_1, region_start_1, region_end_1, Genus) %>% summarise(length = max(length))
unique_hits <- blast_10kb_window_table %>% group_by(chr_1, region_start_1, region_end_1, Genus) %>% count() %>% filter(n == 1) %>% select(-n)
blast_unique_hits <- left_join(unique_hits, blast_10kb_window_table) %>% select(chr_1, region_start_1, region_end_1, Genus, length)
unique_blast_10kb_window_table <- rbind(blast_unique_hits, blast_longest_hits)
Blast_10kb_window_table <- left_join(unique_blast_10kb_window_table, blast_10kb_window_table) %>% distinct()

write.csv(Blast_10kb_window_table, "Master_Blast_10kb_window_table.csv")

#### Plotting #### 
Blast_10kb_window_table <- read.csv("Master_Blast_10kb_window_table.csv") %>%   
  mutate(aln_direction = refend - refstart)
Blast_10kb_window_table$direction <- ifelse(Blast_10kb_window_table$aln_direction>0, 1, -1) #make direction of aln binary
Blast_10kb_window_table %>% 
  ggplot()+
  geom_histogram(aes(x=percent_ID))
Blast_10kb_window_table %>% 
  filter(length > 1000) %>%
  filter(percent_ID > 95) %>%
  ggplot()+
  geom_histogram(aes(x=length, fill = direction, colour = direction, group = direction))+
  xlab("Aligned breakpoint region length (bp)")+ ylab("number of BLAST hits")

Blast_10kb_window_table %>%
  ggplot()+
  geom_histogram(aes(x=aln_direction))

 
