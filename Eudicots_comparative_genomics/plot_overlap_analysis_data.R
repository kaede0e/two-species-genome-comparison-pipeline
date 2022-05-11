#CDS/TE counts in inversion/syntenic/inversion_breakpoints/random_4k region
library(ggplot2)
library(tidyverse)
library(dplyr)

#####Plotting##### 
Master_bedtools_bpoverlap_merged_CDS <- read_csv("Master_bedtools_bpoverlap_CDS.csv")
Master_bedtools_bpoverlap_merged_TE <- read_csv("Master_bedtools_bpoverlap_TE.csv")
Master_bedtools_bpoverlap_breakpoints_merged_cds <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_cds.csv")
Master_bedtools_bpoverlap_breakpoints_merged_TE <- read_csv("Master_bedtools_bpoverlap_breakpoints_merged_TE.csv")

CDS_proportion <- Master_bedtools_bpoverlap_merged_CDS %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(region_start_1, region_end_1, region_length, type, Genus)%>%
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

TE_proportion <- Master_bedtools_bpoverlap_merged_TE %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(region_start_1, region_end_1, region_length, type, Genus)%>%
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

CDS_proportion2 <- Master_bedtools_bpoverlap_breakpoints_merged_cds %>%
  group_by(Genus, type)%>%
  summarise(proportion_of_CDS_per_region = sum(`nCDS(bp)`)/sum(region_length))%>%
  #filter(proportion_of_CDSs <= 0.3)%>% # Arachis was outlier
  ggplot()+
  #geom_point(aes(x = type, y = proportion_of_CDS_per_region, colour = type))+
  #geom_line(aes(x = type, y = proportion_of_CDS_per_region, group = Genus))
  geom_boxplot(aes(x=type, y=proportion_of_CDS_per_region, colour = type))

TE_proportion2 <- Master_bedtools_bpoverlap_breakpoints_merged_TE %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nTE(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_TE_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  ggplot()+
  #geom_point(aes(x = type, y = proportion_of_TE_per_region, colour = type))+
  #geom_line(aes(x = type, y = proportion_of_TE_per_region, group = Genus))
  geom_boxplot(aes(x=type, y=proportion_of_TE_per_region, colour = type))

# combined summary dataframe: 
joined_proportion <- left_join(CDS_proportion, TE_proportion)
joined_proportion2 <- left_join(CDS_proportion2, TE_proportion2)
Master_bedtools_bpoverlap_summary <- rbind(joined_proportion, joined_proportion2) #this summary file is found in the Github

