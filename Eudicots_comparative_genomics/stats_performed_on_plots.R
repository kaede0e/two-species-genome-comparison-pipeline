library(tidyverse)

##### Statistical analysis #####
##For codes used in plotting the data, refer to the scripts deposited in /Eudicots_comparative_genomics under separate analysis title. 

#1) Linear regression analysis to determine statistical significance between number of inversions and: 
#1.1) sequence divergence
Master_data_table_summary <- read_csv("Master_data_table_summary_all_genomes_updated.csv")
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  aov(number_of_inversions ~ percent_genome_divergence,
      data = .)%>% 
  summary()
Master_data_table_summary %>% #proportion of genome inverted was tested - in supplemetary
  filter(ref_or_qry == "ref") %>%
  mutate(proportion_of_inverted_genome = inv_region_length/Genome_length) %>%
  aov(proportion_of_inverted_genome ~ percent_genome_divergence,
      data = .)%>% 
  summary()
#1.2) genome length
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  aov(Genome_length ~ percent_genome_divergence,
      data = .)%>% 
  summary()
#1.3) chromosome length 
inv_plot2 <- Master_score_table_updated %>%
  filter(type == "inv") %>%
  group_by(region_end_1, region_start_1, region_length, Genus, chr_1) %>%
  count()
inversion_stats2_by_chr <- inv_plot2 %>%
  group_by(Genus, chr_1)%>%
  count() %>%
  rename("total_inv_regions_number" = n)
combined_fai <- read_table("chr_length_for_32_genera.txt", col_names = FALSE)#.fai files used to obtain chromosome lengths from all genomes used. 
combined_fai_file <- combined_fai %>% rename("chr_1" = X1, "chr_length" = X2, "Genus" = X3)
plot_chr_length <- left_join(combined_fai_file, inversion_stats2_by_chr)
plot_chr_length %>% 
  lm(total_inv_regions_number ~ chr_length, 
     data=.)%>% 
  summary()
#1.4) geomic TE proportion
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length*100) %>%
  aov(number_of_inversions ~ TE_proportion, 
      data = .) %>%
  summary()

#1+) Rate of inversions fixation calculation 
number_of_inv_by_length%>% #Rate by inversion size categories
  left_join(., genome_divergence) %>%
  group_by(Genus, percent_genome_divergence, region_length_category) %>% 
  summarise(number_of_inversions_by_size = n()) %>%
  mutate(rate_of_inversion = number_of_inversions_by_size/percent_genome_divergence) %>%
  ggplot()+
  geom_boxplot(aes(x=region_length_category, y=rate_of_inversion, fill=region_length_category, alpha = 0.4, colour = region_length_category), outlier.shape = NA)+
  geom_jitter(aes(x=region_length_category, y=rate_of_inversion, colour=region_length_category))+
  theme_classic()+
  theme(legend.position = "na")+
  scale_colour_manual(values = pal5) +
  scale_fill_manual(values = pal5)+
  xlab("Inversion length category")+ ylab("#inversion per % sequence divergence")
number_of_inv_by_length%>%
  left_join(., genome_divergence) %>%
  group_by(Genus, percent_genome_divergence, region_length_category) %>% 
  summarise(number_of_inversions_by_size = n()) %>%
  mutate(rate_of_inversion = number_of_inversions_by_size/percent_genome_divergence)%>%
  group_by(region_length_category)%>%
  summarise(x = quantile(rate_of_inversion, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

#2) One-way ANOVA test to determine any significant relationship between the number of inversion/proportion of genome inverted and the following five parameters: 
#2.1) hybridization
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  left_join(., biological_info)%>%
  aov(number_of_inversions ~ Evidence_of_hybridization, 
      data = .) %>%
  summary()
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(proportion_of_inverted_genome = inv_region_length/Genome_length) %>%
  left_join(., biological_info)%>%
  aov(proportion_of_inverted_genome ~ Evidence_of_hybridization, 
      data = .) %>%
  summary() #Assume this proportion of inverted genome was repeated to every single analysis below. 
#2.2) reproductive strategy
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(proportion_of_inverted_genome = inv_region_length/Genome_length) %>%
  left_join(., biological_info)%>%
  aov(number_of_inversions ~ Reproductive_strategy, 
      data = .) %>%
  summary()
#2.3) assembly method - assembler
assembler <- read_csv("assembler.csv")
assembler %>% separate(Assembler, 
                       into = c("Assembler1", "Assembler2", "Assembler3"), 
                       sep = ", ") %>% 
  mutate(Assembler_category = case_when(
    Assembler1 == "Canu" | Assembler1 == "Celera" | Assembler1 == "Falcon" | Assembler1 == "HiCanu" | Assembler1 == "Hifiasm" | Assembler1 == "DeNovoMAGIC3" ~ "accurate", 
    Assembler1 == "ALLPATHS" | Assembler1 == "Arache" | Assembler1 == "Arachne" | Assembler1 == "SOAPdenovo" | Assembler1 == "Platanus" | Assembler1 == "Newbler" ~ "short-read-based", 
    Assembler1 == "MaSuRCA" | Assembler1 == "Shasta" | Assembler1 == "Flye" | Assembler1 == "MECAT" ~ "less accurate"
  ))
assembler %>%
  aov(number_of_inversions ~ Assembler_combined,
      data = .)%>% #manually Assembler_category from two species was combined to Assembler_combined
  summary()
#2.4) assembly method - sequencing platform
assembly_method %>%
  left_join(., Master_data_table_summary)%>%
  mutate(proportion_of_inverted_genome = inv_region_length/Genome_length) %>%
  aov(number_of_inversions ~ Sequencing_platform,
      data = .)%>% 
  summary()
#2.5) assembly method - physical mapping
assembly_method %>%
  left_join(., Master_data_table_summary)%>%
  mutate(proportion_of_inverted_genome = inv_region_length/Genome_length) %>%
  aov(number_of_inversions ~ Physical_mapping_combined,
      data = .)%>% 
  summary()

#3) paired t-test for proportion of CDS/TE in inversion vs. syntenic region
#3.1) CDS proportion
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
#3.2) TE proportion
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

#4) paired t-test for proportion of CDS/TE in 4kb breakpoint regions vs. genomic CDS/TE proportion
#4.1) CDS proportion 
Master_bedtools_bpoverlap_breakpoints_merged_CDS %>%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>%
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
  t.test(proportion_of_CDS_per_region ~ type, data = ., paired = TRUE)
#4.2) TE proportion
Master_bedtools_bpoverlap_breakpoints_merged_TE %>%
  distinct() %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(type != "alignable_random_4k") %>%
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
  t.test(proportion_of_TE_per_region ~ type, data = ., paired = TRUE)

#5) Inversion breakpoints that overlapped with a gene 
#no stats performed.

#6) Duplicated segments and inverted duplicates identified by BLAST
#no stats performed. 
