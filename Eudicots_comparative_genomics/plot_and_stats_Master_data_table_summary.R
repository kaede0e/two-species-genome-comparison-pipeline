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

#Supplementary Figures 
#S1a) Distribution of syntenic regions identified by SyRI with respect to its sequence identity 

#S1b) Distribution of inverted regions identified by SyRI with respect to its sequence identity 

#S2a) Number of inversions plotted against sequence divergence grouped by self-compatibility category

#S2b) Total number of inversions grouped by self-compatibility category

#S2c) Number of inversions plotted against sequence divergence grouped by hybridization category

#S2d) Total number of inversions grouped by hybridization category
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  filter(!is.na(Evidence_of_hybridization)) %>%
  ggplot()+
  geom_boxplot(aes(x=Evidence_of_hybridization, y=total_inv_regions_number))+ ylab("Number of inversions ≥1000bp")+ xlab("Evidence of hybridization")

