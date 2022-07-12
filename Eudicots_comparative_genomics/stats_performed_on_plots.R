library(tidyverse)

#####Statistical analysis#####

#Figures 
##For codes used in plotting the data, refer to the scripts deposited in /Eudicots_comparative_genomics under separate analysis title. 
#1) Density of sequence identity score inside inverted region (blue) vs. syntenic region (red) for 32 species pair. 
syn_inv_plot %>% 
  filter(Genus =="Fragaria")%>%
  var.test(region_percent_ID ~ type, data = ., alternative = "two.sided") #F-statistics to test variance
  t.test(region_percent_ID ~ type, data = ., paired = FALSE) #unpaired t-test 
                                                            #if the above F-stat shows null - same distribution for both types, which was only true for a handful of species pair

#2) What parameter(s) affects the number of inversions? Number of inversions plotted against
#2a) sequence divergence
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  aov(total_inv_regions_number ~ percent_divergence, 
      data = .) %>%
  summary() #Does the number of inversion depend on species divergence computed from sequence identity? 
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  aov(total_inv_regions_number ~ percent_divergence *
        Evidence_of_hybridization, 
      data = .) %>%
  summary()
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  aov(total_inv_regions_number ~ percent_divergence *
        `Self-compatibility`, 
      data = .) %>%
  summary() #Is there an interacting effect from evidence of hybridization or self-compatibility to the number of inversions?
#2b) genome length 
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  aov(total_inv_regions_number ~ Genome_length, 
      data = .) %>%
  summary()
#2c) TE content 
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  aov(total_inv_regions_number ~ TE_proportion, 
      data = .) %>%
  summary()

#3) paired t-test for proportion of CDS/TE in inversion vs. syntenic region
#3a) CDS proportion
Master_bedtools_bpoverlap_merged_CDS %>%
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
  t.test(proportion_of_CDS_per_region ~ type, data = ., paired = TRUE)
#3b) TE proportion
Master_bedtools_bpoverlap_merged_TE %>%
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
  t.test(proportion_of_TE_per_region ~ type, data = ., paired = TRUE)

#4) paired t-test for proportion of CDS/TE in 4kb breakpoint regions vs. random regions in the genome 
#4a) CDS proportion 
Master_bedtools_bpoverlap_breakpoints_merged_CDS %>%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nCDS(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_CDS_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  t.test(proportion_of_CDS_per_region ~ type, data = ., paired = TRUE)

Master_bedtools_bpoverlap_breakpoints_merged_TE %>%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>% #taking into account for the unaligned regions by SyRI did not make a difference to the analysis so excluded for simplicity.  
  filter(Genus != "Vaccinium") %>%
  filter(Genus != "Medicago") %>%
  filter(Genus != "Quercus") %>% #these three genera were excluded due to inconsistency of the random TE proportion compared to genomic TE content.
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(chr_1, region_start_1, region_end_1, region_length, type, Genus)%>%
  summarise(
    total_bpoverlap = sum(`nTE(bp)`, na.rm = TRUE)
  )%>%
  group_by(Genus, type)%>%
  summarise(
    proportion_of_TE_per_region = sum(total_bpoverlap)/sum(region_length)
  )%>%
  t.test(proportion_of_TE_per_region ~ type, data = ., paired = TRUE)

#5) Inversion breakpoints that overlapped with a gene 
#no stats performed.

#6) Duplicated segments and inverted duplicates identified by BLAST
#no stats performed. 
