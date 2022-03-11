#CDS/TE counts in inversion/syntenic region
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

directory <- "/PATH/TO/PARENT/FOLDER/comparative_genomics" #I have a folder called "comparative_genomics" which contains a subfolder named "Prunus" where I store all the relevant outputs from the pipeline. 
genera <- c("Prunus")

#Import function for CDS
for (genus in genera){
## for inversion data
  bedtools_bpoverlap_merged_cds_inv <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_cds_inv.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           CDS_start = V5, 
           CDS_end = V6, 
           "nCDS(bp)" =V7) %>% 
    add_column(Genus = genus, type = "inv")
  bedtools_count_cds_inv <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_inv.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "inv")

## for synteny data
  bedtools_count_cds_syn <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_syn.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "syn")
  bedtools_bpoverlap_merged_cds_syn <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_cds_syn.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           CDS_start = V5, 
           CDS_end = V6, 
           "nCDS(bp)" =V7) %>% 
    add_column(Genus = genus, type = "syn")
  
}

#Import function for TE
for (genus in genera){
## for inversion data
  bedtools_bpoverlap_merged_TE_inv <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_TE_inv.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           TE_start = V5, 
           TE_end = V6, 
           "nTE(bp)" =V7) %>% 
    add_column(Genus = genus, type = "inv") %>%
    mutate(region_length = region_end_1 - region_start_1)
  bedtools_count_TE_inv <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_inv.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "inv")
  
## for synteny data
  bedtools_bpoverlap_merged_TE_syn <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_TE_syn.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           TE_start = V5, 
           TE_end = V6, 
           "nTE(bp)" =V7) %>% 
    add_column(Genus =genus, type = "syn") %>%
    mutate(region_length = region_end_1 - region_start_1)
  bedtools_count_TE_syn <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_syn.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "syn")

}

#bind dataframes (for Prunus)
Prunus_bedtools_bpoverlap <- rbind(bedtools_bpoverlap_merged_cds_inv, bedtools_bpoverlap_merged_cds_syn)
Prunus_bedtools_count <- rbind(bedtools_count_cds_inv, bedtools_count_cds_syn)
Prunus_bedtools_bpoverlap_TE <- rbind(bedtools_bpoverlap_merged_TE_inv, bedtools_bpoverlap_merged_TE_syn)
Prunus_bedtools_bpoverlap_withTE <- full_join(Prunus_bedtools_bpoverlap, Prunus_bedtools_bpoverlap_TE)
Prunus_bedtools_count_TE <- rbind(bedtools_count_TE_inv, bedtools_count_TE_syn)
Prunus_bedtools_count_withTE <- full_join(Prunus_bedtools_count, Prunus_bedtools_count_TE)

# data analysis and visualization 
## 1) How many CDS/TE are in syn/inv? 
Prunus_bedtools_count_withTE %>%
+ group_by(type)%>%
+ summarise(total_number_of_CDS = sum(`nCDS(seq)`), total_number_of_TE = sum(`nTE(seq)`))


