library(tidyverse)

#Figure 1: Distribution of nucleotide sequence divergence in inversion vs syntenic region
syn_inv_plot <- read_csv("syn_inv_plot.csv")
syn_inv_plot %>% 
  filter(Genus =="Fragaria")%>%
  var.test(region_percent_ID ~ type, data = ., alternative = "two.sided") #F-statistics to test variance
  t.test(region_percent_ID ~ type, data = ., paired = FALSE) #unpaired t-test 
                                                            #if the above F-stat shows null - same distribution for both types, which was only true for a handful of species pair

#Figure 2: Null model for number of inversions
Master_data_table_summary <- read_csv("Master_data_table_summary.csv")
Master_data_table_summary %>%
  left_join(., biological_info)%>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  aov(number_of_inversion_regions ~ TE_proportion *
        `Self-compatibility`,
      data = .) %>%
  summary()

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  filter(!is.na(`Self-compatibility`)) %>%
  ggplot()+
  geom_boxplot(aes(x=`Self-compatibility`, y=number_of_inversion_regions))

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  filter(!is.na(Evidence_of_hybridization)) %>%
  ggplot()+
  geom_boxplot(aes(x=Evidence_of_hybridization, y=number_of_inversion_regions))

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  #filter(Evidence_of_hybridization == "weak")%>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  ggplot(.,aes(x=percent_divergence, y=number_of_inversion_regions)) +
  geom_point(colour = '#00BFC4') + 
  geom_smooth(method="lm", colour = '#00BFC4')

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  #filter(`Self-compatibility`=="mixed")%>%
  aov(number_of_inversion_regions ~ percent_divergence +
        `Self-compatibility`, 
      data= .) %>%
  summary()

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  filter(`Self-compatibility`!="mixed")%>%
  aov(number_of_inversion_regions ~ `Self-compatibility`, 
      data= .) %>%
  summary()

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  aov(inv_region_length ~ syn_divergence + 
        TE_proportion + Genome_length, 
    data = .) %>%
  summary()

#Figure 3: TE proportion / TE length increases with genome size
Master_data_table_summary %>% 
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  lm(TE_proportion ~  Genome_length,
     data=.)

#Figure 4: paired t-test for CDS/TE proportion in inversion vs. syntenic region
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

Master_bedtools_bpoverlap_breakpoints_merged_TE %>%
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
  t.test(proportion_of_TE_per_region ~ type, data = ., paired = TRUE)

#Figure 5: proportion of CDS/TE in breakpoint regions vs. random regions 
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
  filter(type != "alignable_random_4k") %>%
  filter(Genus != "Vaccinium") %>%
  filter(Genus != "Medicago") %>%
  filter(Genus != "Quercus") %>%
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

#Figure 6: gene frequency at breakpoints vs. syntenic region
Master_bedtools_count_gene %>%
  group_by(Genus, type)%>%
  summarise(
    n = n(),
    freq_of_gene = sum(gene_hit)/n
  ) %>%
  t.test(freq_of_gene ~ type, data = ., paired = TRUE)
gene_freq_stats %>% 
  left_join(., biological_info)%>%
  left_join(., Master_data_table_summary)%>% 
  filter(ref_or_qry == "ref") %>%
  #filter(Evidence_of_hybridization == "strong") %>% 
  ggplot()+
  geom_point(aes(x = type, y = freq_of_gene, colour = Evidence_of_hybridization))+
  geom_line(aes(x = type, y = freq_of_gene, group = Genus))
  aov(freq_of_gene ~ type * Genus,
      data = .) %>%
  summary()


