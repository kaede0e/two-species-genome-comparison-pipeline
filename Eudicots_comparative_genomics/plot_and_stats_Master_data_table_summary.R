#what parameter(s) affects the number of inversions? 
library(tidyverse)

#####Plotting and statistical analysis#####
Master_data_table_summary <- read_csv("Master_data_table_summary_all_genomes_updated.csv")
Master_data_table_summary %>%
  left_join(., biological_info)%>%
  #filter(Evidence_of_hybridization == "weak")%>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  ggplot(.,aes(x=percent_divergence, y=number_of_inversion_regions)) +
  geom_point() + 
  geom_smooth(method="lm")

Master_data_table_summary %>%
  left_join(., biological_info)%>%
  #filter(`Self-compatibility`=="mixed")%>%
  aov(number_of_inversion_regions ~ percent_divergence +
        `Self-compatibility`, 
      data= .) %>%
  summary()

Master_data_table_summary %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  aov(inv_region_length ~ syn_divergence + 
        TE_proportion +
        Genome_length, 
    data = .) %>%
  summary()

#Figures
#1) Density of sequence identity score inside inverted region (blue) vs. syntenic region (red) for 32 species pair. 
syn_inv_plot %>%
  ggplot() +
  geom_density(aes(x=region_percent_ID, fill=type, alpha=0.4)) +
  scale_fill_manual(values = pnw_palette("Bay",2)) +
  facet_wrap("Genus", nrow = 5)

#2) Number of inversions plotted against
#2a) sequence divergence
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>% 
  ggplot()+ 
  geom_point(aes(x=percent_divergence, y=total_inv_regions_number))+
  geom_smooth(aes(x=percent_divergence, y=total_inv_regions_number), method = "lm")
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

#Supplementary Figures 
#S1a) Distribution of syntenic regions identified by SyRI with respect to its sequence identity 
syn_plot %>%
  ggplot() +
  geom_histogram(aes(x=region_percent_ID), binwidth = 0.002) +
  facet_wrap("Genus", nrow = 5) +
  geom_vline(aes(xintercept=syn_divergence))

#S1b) Distribution of inverted regions identified by SyRI with respect to its sequence identity 
syn_inv_plot %>%
  filter(type == "inv") %>%
  ggplot() +
  geom_histogram(aes(x=region_percent_ID), binwidth = 0.002) +
  facet_wrap("Genus", nrow = 5) 

#S2a) Number of inversions plotted against sequence divergence grouped by hybridization category
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  ggplot()+
  geom_point(aes(x=percent_divergence, y=total_inv_regions_number, colour = Evidence_of_hybridization))+ 
  geom_smooth(aes(x=percent_divergence, y=total_inv_regions_number, colour = Evidence_of_hybridization), method = "lm")+ 
  ylab("Number of inversions ≥1000bp") +xlab("Sequence divergence (%)")

#S2b) Total number of inversions grouped by hybridization category
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  filter(!is.na(Evidence_of_hybridization)) %>%
  ggplot()+
  geom_boxplot(aes(x=Evidence_of_hybridization, y=total_inv_regions_number))+ ylab("Number of inversions ≥1000bp")+ xlab("Evidence of hybridization")

#S2c) Number of inversions plotted against sequence divergence grouped by self-compatibility category
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  ggplot()+
  geom_point(aes(x=percent_divergence, y=total_inv_regions_number, colour = `Self-compatibility`))+ 
  geom_smooth(aes(x=percent_divergence, y=total_inv_regions_number, colour = `Self-compatibility`), method = "lm")+ 
  ylab("Number of inversions ≥1000bp") +xlab("Sequence divergence (%)")

#S2d) Total number of inversions grouped by self-compatibility category
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  filter(!is.na(`Self-compatibility`)) %>%
  ggplot()+
  geom_boxplot(aes(x=`Self-compatibility`, y=total_inv_regions_number))+ ylab("Number of inversions ≥1000bp")

#S3a) Genome length vs TE length 
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>% 
  ggplot()+ geom_point(aes(x=Genome_length, y=total_TE_length)) +geom_smooth(aes(x=Genome_length, y=total_TE_length), method = "lm")
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  aov(total_inv_regions_number ~ total_TE_length, 
      data = .) %>%
  summary()

#S3b) Genome length vs genomic TE proportion
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>% 
  ggplot()+ geom_point(aes(x=Genome_length, y=TE_proportion)) +geom_smooth(aes(x=Genome_length, y=TE_proportion), method = "lm")
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  aov(total_inv_regions_number ~ TE_proportion, 
      data = .) %>%
  summary()

#S3c) genomic TE content by genus 
genomic_TE_proportion <- Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(genomic_TE_proportion = total_TE_length / Genome_length)%>%
  ggplot(., aes(x=reorder(Genus,-genomic_TE_proportion), y=genomic_TE_proportion))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Genus")
