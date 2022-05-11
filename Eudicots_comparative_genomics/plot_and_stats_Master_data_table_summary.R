library(tidyverse)

# 
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
