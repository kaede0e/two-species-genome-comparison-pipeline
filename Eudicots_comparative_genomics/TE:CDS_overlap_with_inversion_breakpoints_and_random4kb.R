#CDS/TE in inversion_breakpoint_4k/random_4k/inversion_breakpoints
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

datalist_bedtools_bpoverlap = list()
datalist_bedtools_count = list()#not used for our purposes

# we are looking for the proportion of CDS or TE covering the genome (relevant for bedtools intersect -wo and not -c count)
#### CDS #####
#Overlap with inv breakpoint 4kbp regions vs. genomic CDS proportion
for (genus in genera){
  ## for inversion breakpoints data
  bedtools_bpoverlap_breakpoint_regions_bed <- read.table(paste0(directory, "/", genus, "/breakpoints_regions_onesided.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  bedtools_bpoverlap_merged_cds_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_cds_inv_breakpoint_onesided.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           CDS_start = V5, 
           CDS_end = V6, 
           "nCDS(bp)" =V7) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  bedtools_bpoverlap_merged_cds_inv_breakpoints_including_zeros <- full_join(bedtools_bpoverlap_breakpoint_regions_bed, bedtools_bpoverlap_merged_cds_inv_breakpoints)
  datalist_bedtools_bpoverlap[[genus]] <- bedtools_bpoverlap_merged_cds_inv_breakpoints_including_zeros # add it to your list
  
}
Master_bedtools_bpoverlap_breakpoints_merged_CDS = do.call(rbind, datalist_bedtools_bpoverlap)

##### TE #####
#Overlap with inv breakpoints vs. genomic TE proportion
for (genus in genera){
  ## for inversion breakpoints data
  bedtools_bpoverlap_breakpoint_regions_bed <- read.table(paste0(directory, "/", genus, "/breakpoints_regions_onesided.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  bedtools_bpoverlap_merged_TE_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_TE_inv_breakpoints_onesided.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           TE_start = V5, 
           TE_end = V6, 
           "nTE(bp)" =V7) %>% 
    add_column(Genus = genus, type = "inv_breakpoints") %>%
    mutate(region_length = region_end_1 - region_start_1)
  bedtools_bpoverlap_merged_TE_inv_breakpoints_including_zeros <- full_join(bedtools_bpoverlap_breakpoint_regions_bed, bedtools_bpoverlap_merged_TE_inv_breakpoints)
  bedtools_count_TE_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_inv_breakpoints.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  datalist_bedtools_bpoverlap[[genus]] <- bedtools_bpoverlap_merged_TE_inv_breakpoints_including_zeros # add it to your list
  
}
Master_bedtools_bpoverlap_breakpoints_merged_TE = do.call(rbind, datalist_bedtools_bpoverlap)

##### Plotting ##### 
Master_bedtools_bpoverlap_breakpoints_merged_CDS <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_CDS_onesided.csv") 
Master_bedtools_bpoverlap_breakpoints_merged_TE <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_TE_onesided.csv") 
  
Master_bedtools_bpoverlap_breakpoints_merged_CDS <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_CDS_final.csv") #when 4kb of either ends are considered as breakpoint regions
Master_bedtools_bpoverlap_breakpoints_merged_TE <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_TE_final.csv") 
  
CDS_proportion2 <- Master_bedtools_bpoverlap_breakpoints_merged_CDS %>% #if using the genomic CDS%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>% #extra analyses we did but omitted.
  filter(type != "random_4k")%>% 
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nCDS(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_CDS_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  rbind(., genomic_CDS_proportion)%>% 
  ggplot()+
  geom_point(aes(x = type, y = proportion_of_CDS_per_region, colour = type))+
  geom_line(aes(x = type, y = proportion_of_CDS_per_region, group = Genus))+
  theme_classic()+
  scale_color_manual(values = pal3)+
  scale_x_discrete(labels = c("Inversion breakpoints", "Whole genome"))+
  xlab("")+ ylab("Proportion of CDS")

TE_proportion2 <- Master_bedtools_bpoverlap_breakpoints_merged_TE %>% #if using the whole genome TE%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>% #extra analyses we did but omitted. 
  filter(type != "random_4k") %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nTE(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_TE_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>% 
  rbind(., genomic_TE_proportion) %>%
  ggplot()+
  geom_point(aes(x = type, y = proportion_of_TE_per_region, colour = type))+
  geom_line(aes(x = type, y = proportion_of_TE_per_region, group = Genus))+
  theme_classic()+
  scale_color_manual(values = pal3)+
  scale_x_discrete(labels = c("Inversion breakpoints", "Whole genome"))+
  xlab("")+ ylab("Proportion of TE")

joined_proportion <- left_join(CDS_proportion, TE_proportion)
joined_proportion2 <- left_join(CDS_proportion2, TE_proportion2)
Master_bedtools_bpoverlap_summary <- rbind(joined_proportion, joined_proportion2) 


