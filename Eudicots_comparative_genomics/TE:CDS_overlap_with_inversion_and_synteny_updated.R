#CDS/TE in inversion/syntenic region
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
            "Ipomoea", "Juglans", "Malus", "Medicago", 
            "Phaseolus", "Populus", "Prunus", "Pyrus", 
            "Quercus", "Rhododendron", "Rosa", "Rubus", 
            "Salix", "Salvia", "Solanum", 
            "Vaccinium", "Vigna", "Vitis")
directory <- "/Volumes/Backup Plus/Comparative genomics - inversion genomics/Scripts/Shell scripts/Comparative genomics"

#### CDS #####

#Overlap with SyRI inv/syn region
datalist_bedtools_bpoverlap = list()
datalist_bedtools_count = list()
for (genus in genera){
## for inversion data
  bedtools_bpoverlap_inv_bed <- read.table(paste0(directory, "/", genus, "/inversion_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "inv")
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
  bedtools_bpoverlap_merged_cds_inv_including_zeros <- full_join(bedtools_bpoverlap_inv_bed, bedtools_bpoverlap_merged_cds_inv)
  bedtools_count_cds_inv <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_inv.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "inv")
## for synteny data
  bedtools_bpoverlap_syn_bed <- read.table(paste0(directory, "/", genus, "/syntenic_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    )%>%
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
  bedtools_bpoverlap_merged_cds_syn_including_zeros <- full_join(bedtools_bpoverlap_syn_bed, bedtools_bpoverlap_merged_cds_syn)
  bedtools_count_cds_syn <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_syn.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "syn")
  bedtools_bpoverlap_merged_both <- rbind(bedtools_bpoverlap_merged_cds_inv_including_zeros,bedtools_bpoverlap_merged_cds_syn_including_zeros)
  datalist_bedtools_bpoverlap[[genus]] <- bedtools_bpoverlap_merged_both # add it to your list
  bedtools_count_both <- rbind(bedtools_count_cds_inv,bedtools_count_cds_syn)
  datalist_bedtools_count[[genus]] <- bedtools_count_both # add it to your list
  
}
Master_bedtools_bpoverlap_merged_CDS = do.call(rbind, datalist_bedtools_bpoverlap)
Master_bedtools_count_CDS = do.call(rbind, datalist_bedtools_count)

##### TE #####

datalist <- list()
for (genus in genera){
  ## for raw TE merged if duplicates
  raw_TE <- read.table(paste0(directory, "/", genus, "/merged_TE.bed"), head = FALSE) %>%
    rename(chr_1 = V1,
           TE_start_1 = V2, 
           TE_end_1 = V3) %>%
    add_column(Genus = genus)%>%
    mutate(TE_length = TE_end_1 - TE_start_1)
  datalist[[genus]] <- raw_TE
}
raw_TE_data = do.call(rbind, datalist)
genome_length_subset <- genome_length %>%
  filter(ref_or_qry == "ref")
Master_data_table_TE <- left_join(raw_TE_summarised_data, genome_length_subset)
Master_data_table_summary <- left_join(Master_data_table_TE, divergence) #combined with divergence data

#Overlap with SyRI inv/syn region
datalist_bedtools_bpoverlap = list()
datalist_bedtools_count = list()
for (genus in genera){
## for inversion data
  bedtools_bpoverlap_inv_bed <- read.table(paste0(directory, "/", genus, "/inversion_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "inv")
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
  bedtools_bpoverlap_merged_TE_inv_including_zeros <- full_join(bedtools_bpoverlap_inv_bed, bedtools_bpoverlap_merged_TE_inv)
  bedtools_count_TE_inv <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_inv.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "inv")
## for synteny data
  bedtools_bpoverlap_syn_bed <- read.table(paste0(directory, "/", genus, "/syntenic_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    )%>%
    add_column(Genus = genus, type = "syn")
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
  bedtools_bpoverlap_merged_TE_syn_including_zeros <- full_join(bedtools_bpoverlap_syn_bed, bedtools_bpoverlap_merged_TE_syn)
  bedtools_count_TE_syn <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_syn.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "syn")
  bedtools_bpoverlap_merged_both <- rbind(bedtools_bpoverlap_merged_TE_inv_including_zeros,bedtools_bpoverlap_merged_TE_syn_including_zeros)
  datalist_bedtools_bpoverlap[[genus]] <- bedtools_bpoverlap_merged_both # add it to your list
  bedtools_count_both <- rbind(bedtools_count_TE_inv,bedtools_count_TE_syn)
  datalist_bedtools_count[[genus]] <- bedtools_count_both # add it to your list
  
}
Master_bedtools_bpoverlap_merged_TE = do.call(rbind, datalist_bedtools_bpoverlap)
Master_bedtools_count_TE = do.call(rbind, datalist_bedtools_count)


#####Plotting##### 
#Make sure it's `nCDS(bp)` not single quote!! ''
Master_bedtools_bpoverlap_merged_CDS <- read_csv("Master_bedtools_bpoverlap_CDS_final.csv")
Master_bedtools_bpoverlap_merged_TE <- read_csv("Master_bedtools_bpoverlap_TE_final.csv")

CDS_proportion <- Master_bedtools_bpoverlap_merged_CDS %>%
  distinct() %>% 
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nCDS(bp)`, na.rm = TRUE)
  )%>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_CDS_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  ggplot()+
  geom_point(aes(x = type, y = proportion_of_CDS_per_region, colour = type))+
  geom_line(aes(x = type, y = proportion_of_CDS_per_region, group = Genus))+
  scale_color_manual(values = pnw_palette("Bay",2))

TE_proportion <- Master_bedtools_bpoverlap_merged_TE %>%
  distinct() %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nTE(bp)`, na.rm = TRUE)
  )%>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_TE_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  ggplot()+
  geom_point(aes(x = type, y = proportion_of_TE_per_region, colour = type))+
  geom_line(aes(x = type, y = proportion_of_TE_per_region, group = Genus))+
  scale_color_manual(values = pnw_palette("Bay",2))



