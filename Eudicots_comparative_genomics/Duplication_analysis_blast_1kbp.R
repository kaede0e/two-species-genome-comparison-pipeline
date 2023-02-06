#Segmental duplications at inversion breakpoints within & between downstream/upstream
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

#Data Import function
## Define directory with your data formatted as directory/genera/
genera <- c("Acer", "Actinidia", "Arabis", "Arachis",
            "Citrus", "Corylus", "Corymbia", "Cucumis", 
            "Eukalyptus", "Fragaria", "Glycine", "Gossypium", 
            "Ipomoea", "Juglans", "Luffa", "Malus", "Medicago", 
            "Phaseolus", "Populus", "Prunus", "Pyrus", 
            "Quercus", "Rhododendron", "Raphanus", "Rosa", "Rubus", 
            "Salix", "Salvia", "Solanum", 
            "Vaccinium", "Vigna", "Vitis")

#Data Import
#### Blast output #####
blasted_1k_coordinates <- read_csv("blasted_1kbp_coordinates.csv") %>% 
    rename(region_start_1 = start, region_end_1 = end)
blast_output_1kbp <- read.table("Master_inv_breakpoints_within_aln_blast_output_1kb_window.txt") %>% 
  as_tibble(.)%>%
  distinct()%>%
  rename(Genus = V1, 
         qseqID = V2, 
         refseqID = V3, 
         percent_ID = V4, 
         length= V5,
         mismatch = V6,
         gap_open = V7,
         qstart = V8,
         qend = V9, 
         refstart = V10,
         refend = V11,
         E_value = V12, 
         bitscore = V13)%>%
  separate(qseqID, into = c("chr_1", "coordinates"), sep = ":") %>%
  separate(refseqID, into = c("chr_2", "coordinates2"), sep = ":") %>%
  separate(coordinates, into = c("region_start_1", "region_end_1"), sep = "-", convert = TRUE)%>%
  separate(coordinates2, into = c("region_start_2", "region_end_2"), sep = "-", convert = TRUE)


blast_output_table_joined <- left_join(blasted_1k_coordinates, blast_output_1kbp)

blast_output_table_joined %>% 
  filter(length <1000) %>% #matched length should not exceed 1kbp
  ggplot()+
  geom_histogram(aes(x=percent_ID))+
  xlab("sequence similarity of 10kb window (%)")+ ylab("count")

duplicates <- blast_output_table_joined %>% group_by(Genus, chr_1, region_start_1, region_end_1, chr_2, region_start_2, region_end_2) %>% count() %>% filter(n != 1)
duplicate_hits <- left_join(blast_output_table_joined, duplicates)
blast_longest_hits <- duplicate_hits %>% filter(n >= 1) %>% group_by(Genus, chr_1, region_start_1, region_end_1, chr_2, region_start_2, region_end_2) %>% summarise(length = max(length))
unique_hits <- blast_output_table_joined %>% group_by(Genus, chr_1, region_start_1, region_end_1, chr_2, region_start_2, region_end_2, length) %>% count() %>% filter(n == 1) %>% select(-n)
blast_unique_hits <- left_join(unique_hits, blast_output_table_joined) %>% select(Genus, chr_1, chr_1, region_start_1, region_end_1, chr_2, region_start_2, region_end_2, length)
unique_blast_1kb_aln_table <- rbind(blast_unique_hits, blast_longest_hits)
Blast_1kb_aln_table <- left_join(unique_blast_1kb_aln_table, blast_output_table_joined) %>% distinct()

write.csv(Blast_1kb_aln_table, "Master_inv_breakpoints_within_aln_blast_output_1kb_unique_hits_with_refseq.csv")

Blast_1kb_aln_table <- read.csv("Master_inv_breakpoints_within_aln_blast_output_1kb_unique_hits_with_refseq.csv") %>%   
  as_tibble(.)%>%
  mutate(aln_direction = refend - refstart)
Blast_1kb_aln_table$direction <- ifelse(Blast_1kb_aln_table$aln_direction>0, 1, -1) #make direction of aln binary
Blast_1kb_aln_table %>% 
  filter(length < 1000) %>%
  ggplot()+
  geom_histogram(aes(x=percent_ID))
Blast_1kb_aln_table %>% 
  filter(length < 1000) %>%
  ggplot()+
  geom_histogram(aes(x=length))

Blast_1kb_10kb_aln_table <- full_join(blasted_coordinates_summary, Blast_1kb_aln_table)%>% 
  select(Genus, chr_1, region_start_1, region_end_1, region_start_2, region_end_2, tenk_start, tenk_end, length, percent_ID, direction)


#### Plotting #### 
Blast_1kb_10kb_aln_table <- read.csv("Master_Blast_1kb_10kb_aln_table.csv") %>% as_tibble()
Blast_1kb_10kb_aln_table %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(length < 1000) %>% #matched length should not exceed 1kbp
  filter(region_start_1 >=0) %>% #any 1kbp should not start below zero
  filter(region_start_2 >=0) %>%
  filter(region_start_1 != region_start_2)%>% #remove self aln
  group_by(Genus, chr_1, tenk_start, tenk_end)%>% count() %>%
  mutate(BLAST_result = case_when(
    n <= 10 ~ "No BLAST hit",
    n >= 11 ~ "1 or more BLAST hits",
  ))%>%
  group_by(BLAST_result)%>% count()%>%
  ggplot(aes(x=BLAST_result, y=n, fill = BLAST_result))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = c("#509D85", "#999999"))+
  theme_classic()+
  theme(legend.position = "na")+
  xlab("")+ ylab("Number of BLAST hits within 10 kbp breakpoint regions")
