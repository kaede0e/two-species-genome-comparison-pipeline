library(dplyr)
library(tidyverse)
library(ggplot2)

#Supplementary Figures with codes for plots
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
  lm(total_TE_length ~  Genome_length,
     data=.)

#S3b) Genome length vs genomic TE proportion
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>% 
  ggplot()+ geom_point(aes(x=Genome_length, y=TE_proportion)) +geom_smooth(aes(x=Genome_length, y=TE_proportion), method = "lm")
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length) %>%
  lm(TE_proportion ~  Genome_length,
     data=.)

#S3c) genomic TE content by genus 
genomic_TE_proportion <- Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(genomic_TE_proportion = total_TE_length / Genome_length)%>%
  ggplot(., aes(x=reorder(Genus,-genomic_TE_proportion), y=genomic_TE_proportion))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Genus")


