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

#Extract inversion breakpoint regions in 1kbp small pieces (I can't figure out loop)
  end_1k_1 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1,
              end = region_end_1 + 1000)
  end_1k_2 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 1000,
              end = region_end_1 + 2000)
  end_1k_3 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 2000,
              end = region_end_1 + 3000)
  end_1k_4 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 3000,
              end = region_end_1 + 4000)
  end_1k_5 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 4000,
              end = region_end_1 + 5000)
  end_1k_6 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 5000,
              end = region_end_1 + 6000)
  end_1k_7<- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 6000,
              end = region_end_1 + 7000)
  end_1k_8 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 7000,
              end = region_end_1 + 8000)
  end_1k_9 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 8000,
              end = region_end_1 + 9000)
  end_1k_10 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 9000,
              end = region_end_1 + 10000)
  end_1k_coordinates <- rbind(end_1k_1, end_1k_2, end_1k_3, end_1k_4, end_1k_5, end_1k_6, end_1k_7, end_1k_8, end_1k_9, end_1k_10)

  start_1k_1 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 1000,
              upper_2k = region_start_1)
  start_1k_2 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 2000,
              upper_2k = region_start_1 - 1000)
  start_1k_3 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 3000,
              upper_2k = region_start_1 - 2000)
  start_1k_4 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 4000,
              upper_2k = region_start_1 - 3000)
  start_1k_5 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 5000,
              upper_2k = region_start_1 - 4000)
  start_1k_6 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 6000,
              upper_2k = region_start_1 - 5000)
  start_1k_7 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 7000,
              upper_2k = region_start_1 - 6000)
  start_1k_8 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 8000,
              upper_2k = region_start_1 - 7000)
  start_1k_9 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 9000,
              upper_2k = region_start_1 - 8000)
  start_1k_10 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 10000,
              upper_2k = region_start_1 - 9000)
  start_1k_coordinates <- rbind(start_1k_1, start_1k_2, start_1k_3, start_1k_4, start_1k_5, start_1k_6, start_1k_7, start_1k_8, start_1k_9, start_1k_10)
  
for (genus in genera){
  end_1k_coordinates %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_end_1k_coordinates_', genus, '.csv'))
}
for (genus in genera){
  start_1k_coordinates %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_start_1k_coordinates_', genus, '.csv'))
}


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


