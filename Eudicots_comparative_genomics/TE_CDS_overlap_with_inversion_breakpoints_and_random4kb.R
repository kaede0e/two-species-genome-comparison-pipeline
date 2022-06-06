#CDS/TE in inversion_breakpoint_4k/random_4k/inversion_breakpoints
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)
getwd()
setwd("/Volumes/Backup Plus/Comparative genomics - inversion genomics/Scripts/Shell scripts/Comparative genomics")

#Data Import function 
datalist_bedtools_bpoverlap = list()
datalist_bedtools_count = list()

#### CDS #####

#Overlap with inv breakpoint 4kbp regions vs. random 4kbp regions on the reference genome
for (genus in genera){
  ## for inversion breakpoints data
  bedtools_bpoverlap_breakpoint_regions_bed <- read.table(paste0(directory, "/", genus, "/breakpoints_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  bedtools_bpoverlap_merged_cds_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_cds_inv_breakpoints.txt"), head = FALSE) %>%
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
  bedtools_count_cds_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_inv_breakpoints.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  ## for random 4kbp regions data
  bedtools_bpoverlap_random_regions_bed <- read.table(paste0(directory, "/", genus, "/bedtools_random_4k_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "random_4k")
  bedtools_bpoverlap_merged_cds_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_cds_random_4k_wholechr.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           CDS_start = V5, 
           CDS_end = V6, 
           "nCDS(bp)" =V7) %>% 
    add_column(Genus = genus, type = "random_4k")
  bedtools_bpoverlap_merged_cds_random_4k_including_zeros <- full_join(bedtools_bpoverlap_random_regions_bed, bedtools_bpoverlap_merged_cds_random_4k)
  bedtools_count_cds_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_random_4k_wholechr.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "random_4k")
  ## for alignable random 4kbp regions data
  bedtools_bpoverlap_random_aligned_regions_bed <- read.table(paste0(directory, "/", genus, "/bedtools_random_4k_aligned_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "alignable_random_4k")
  bedtools_bpoverlap_merged_cds_alignable_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_cds_random_4k.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           CDS_start = V5, 
           CDS_end = V6, 
           "nCDS(bp)" =V7) %>% 
    add_column(Genus = genus, type = "alignable_random_4k")
  bedtools_bpoverlap_merged_cds_alignable_random_4k_including_zeros <- full_join(bedtools_bpoverlap_random_aligned_regions_bed, bedtools_bpoverlap_merged_cds_alignable_random_4k)
  bedtools_count_cds_alignable_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_random_4k.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nCDS(seq)" = V4) %>% 
    add_column(Genus = genus, type = "alignablerandom_4k")
  bedtools_bpoverlap_merged_both <- rbind(
    bedtools_bpoverlap_merged_cds_inv_breakpoints_including_zeros,
    bedtools_bpoverlap_merged_cds_random_4k_including_zeros,
    bedtools_bpoverlap_merged_cds_alignable_random_4k_including_zeros)
  datalist_bedtools_bpoverlap[[genus]] <- bedtools_bpoverlap_merged_both # add it to your list
  bedtools_count_both <- rbind(bedtools_count_cds_inv_breakpoints,bedtools_count_cds_random_4k,bedtools_count_cds_alignable_random_4k)
  datalist_bedtools_count[[genus]] <- bedtools_count_both # add it to your list
  
}
Master_bedtools_bpoverlap_breakpoints_merged_CDS = do.call(rbind, datalist_bedtools_bpoverlap)
Master_bedtools_count_breakpoints_CDS = do.call(rbind, datalist_bedtools_count)


##### TE #####
# we are looking for the proportion of TE covering the genome (relevant for bedtools intersect -wo) you can ignore count (-c)
#Overlap with inv breakpoints vs. random 4kbp regions on the reference genome
for (genus in genera){
  ## for inversion breakpoints data
  bedtools_bpoverlap_breakpoint_regions_bed <- read.table(paste0(directory, "/", genus, "/breakpoints_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "inv_breakpoints")
  bedtools_bpoverlap_merged_TE_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_TE_inv_breakpoints.txt"), head = FALSE) %>%
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
  ## for random 4kb data
  bedtools_bpoverlap_random_regions_bed <- read.table(paste0(directory, "/", genus, "/bedtools_random_4k_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "random_4k")
  bedtools_bpoverlap_merged_TE_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_TE_random_4k_wholechr.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           TE_start = V5, 
           TE_end = V6, 
           "nTE(bp)" =V7) %>% 
    add_column(Genus =genus, type = "random_4k") %>%
    mutate(region_length = region_end_1 - region_start_1)
  bedtools_bpoverlap_merged_TE_random_4k_including_zeros <- full_join(bedtools_bpoverlap_random_regions_bed, bedtools_bpoverlap_merged_TE_random_4k)
  bedtools_count_TE_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_random_4k_wholechr.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "random_4k")
  ## for alignable random 4kbp regions data
  bedtools_bpoverlap_random_aligned_regions_bed <- read.table(paste0(directory, "/", genus, "/bedtools_random_4k_aligned_regions.bed"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
    ) %>% 
    add_column(Genus = genus, type = "alignable_random_4k")
  bedtools_bpoverlap_merged_TE_alignable_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_bpoverlap_merged_TE_random_4k.txt"), head = FALSE) %>%
    select(V1, V2, V3, V4, V5, V6, V7) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           chr_1b = V4, 
           TE_start = V5, 
           TE_end = V6, 
           "nTE(bp)" =V7) %>% 
    add_column(Genus = genus, type = "alignable_random_4k")%>%
    mutate(region_length = region_end_1 - region_start_1) 
  bedtools_bpoverlap_merged_TE_alignable_random_4k_including_zeros <- full_join(bedtools_bpoverlap_random_aligned_regions_bed, bedtools_bpoverlap_merged_TE_alignable_random_4k)
  bedtools_count_TE_alignable_random_4k <- read.table(paste0(directory, "/", genus, "/bedtools_count_TE_random_4k.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "nTE(seq)" = V4) %>% 
    add_column(Genus = genus, type = "alignablerandom_4k")
  bedtools_bpoverlap_merged_both <- rbind(
    bedtools_bpoverlap_merged_TE_inv_breakpoints_including_zeros,
    bedtools_bpoverlap_merged_TE_random_4k_including_zeros,
    bedtools_bpoverlap_merged_TE_alignable_random_4k_including_zeros)
  datalist_bedtools_bpoverlap[[genus]] <- bedtools_bpoverlap_merged_both # add it to your list
  bedtools_count_both <- rbind(bedtools_count_TE_inv_breakpoints,bedtools_count_TE_random_4k,bedtools_count_TE_alignable_random_4k)
  datalist_bedtools_count[[genus]] <- bedtools_count_both # add it to your list
  
}
Master_bedtools_bpoverlap_breakpoints_merged_TE = do.call(rbind, datalist_bedtools_bpoverlap)
Master_bedtools_count_breakpoints_TE = do.call(rbind, datalist_bedtools_count)

#####Plotting##### 
Master_bedtools_bpoverlap_breakpoints_merged_CDS <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_CDS_final.csv")
Master_bedtools_bpoverlap_breakpoints_merged_TE <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_TE_final.csv")
  
CDS_proportion2 <- Master_bedtools_bpoverlap_breakpoints_merged_CDS %>%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "random_4k") %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nCDS(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_CDS_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  ggplot()+
  geom_point(aes(x = type, y = proportion_of_CDS_per_region, colour = type))+
  geom_line(aes(x = type, y = proportion_of_CDS_per_region, group = Genus))

TE_proportion2 <- Master_bedtools_bpoverlap_breakpoints_merged_TE %>%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nTE(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_TE_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  ggplot()+
  geom_point(aes(x = type, y = proportion_of_TE_per_region, colour = type))+
  geom_line(aes(x = type, y = proportion_of_TE_per_region, group = Genus))

joined_proportion <- left_join(CDS_proportion, TE_proportion)
joined_proportion2 <- left_join(CDS_proportion2, TE_proportion2)
Master_bedtools_bpoverlap_summary <- rbind(joined_proportion, joined_proportion2) 

