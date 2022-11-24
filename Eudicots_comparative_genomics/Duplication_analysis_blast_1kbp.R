#Segmental duplications at inversion breakpoints within & between downstream/upstream
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)
getwd()
setwd("/Volumes/Backup Plus/Comparative genomics - inversion genomics/Scripts/Shell scripts/Comparative genomics")

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
blasted_1k_coordinates <- read_csv("blasted_1kbp_coordinates.csv")
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
  separate(coordinates, into = c("start", "end"), sep = "-", convert = TRUE)


blast_output_table_joined <- left_join(blasted_1k_coordinates, blast_output_1kbp)

blast_output_table_joined %>% 
  filter(length <1000) %>% #matched length should not exceed 1kbp
  ggplot()+
  geom_histogram(aes(x=percent_ID))+
  xlab("sequence similarity of 10kb window (%)")+ ylab("count")

duplicates <- blast_output_table_joined %>% group_by(Genus, chr_1, start, end, refseqID) %>% count() %>% filter(n != 1)
duplicate_hits <- left_join(blast_output_table_joined, duplicates)
blast_longest_hits <- duplicate_hits %>% filter(n >= 1) %>% group_by(Genus, chr_1, start, end, refseqID) %>% summarise(length = max(length))
unique_hits <- blast_output_table_joined %>% group_by(Genus, chr_1, start, end, refseqID, length) %>% count() %>% filter(n == 1) %>% select(-n)
blast_unique_hits <- left_join(unique_hits, blast_output_table_joined) %>% select(Genus, chr_1, start, end, refseqID, length)
unique_blast_1kb_aln_table <- rbind(blast_unique_hits, blast_longest_hits)
Blast_1kb_aln_table <- left_join(unique_blast_1kb_aln_table, blast_output_table_joined) %>% distinct()

write.csv(Blast_1kb_aln_table, "Master_inv_breakpoints_within_aln_blast_output_1kb_unique_hits.csv")

#### Plotting #### 
Blast_1kb_aln_table <- read.csv("Master_inv_breakpoints_within_aln_blast_output_1kb_unique_hits.csv") %>%   
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

Blast_1kb_aln_table %>% 
  filter(length < 1000) %>% #matched length should not exceed 1kbp
  filter(length > 10) %>%
  filter(percent_ID > 95) %>%
  ggplot()+
  geom_histogram(aes(x=length, fill = direction,  group = direction, alpha = 0.4))+
  theme_classic()+
  xlab("Aligned breakpoint region length (bp)")+ ylab("Number of BLAST hits")+
  scale_fill_gradientn(colours = pal)

Blast_1kb_aln_table %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(BLAST_result = case_when(
    length == 0 ~ "No blast hit",
    direction == 1 ~ "Forward", 
    direction == -1 ~ "Reverse"
  ))%>%
  group_by(BLAST_result)%>% count()%>%
  ggplot(aes(x=BLAST_result, y=n, fill = BLAST_result))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = c("#509D85", "#999999", "#B2C15D"))+
  scale_x_discrete(limits=c("Forward", "Reverse", "No blast hit"))+
  theme_classic()+
  theme(legend.position = "na")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab("Alignment")+ ylab("Number of BLAST hits within 10 kbp breakpoint regions")+ labs(fill = "Alignment") 


