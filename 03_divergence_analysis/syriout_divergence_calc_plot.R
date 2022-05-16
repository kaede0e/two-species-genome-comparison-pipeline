library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

alignment_score <- read.table("table.txt", head=TRUE) %>%
  rename(block_start_1 = X0, 
         block_end_1 =X1, 
         block_start_2 =X2, 
         block_end_2 = X3, 
         length_1 = X4, 
         length_2 = X5, 
         percent_ID = X6,
         direction_1 = X7, 
         direction_8 = X8, 
         chr_1 = X9, 
         chr_2 = X10,
         alignment= X11) %>% 
  select(-alignment)
alignment_score$Genus <- "Genus1"
inversion_blocks = read.table("invOut_table.txt", head = TRUE)
inversion_blocks$type <- "inv"
syntenic_blocks = read.table("synOut_table.txt", head = TRUE)
syntenic_blocks$type <- "syn"
inversion_blocks <- inversion_blocks %>% select(-chr_2)
syntenic_blocks <- syntenic_blocks %>% select(-chr_2)
inversion_joined <- inner_join(alignment_score, inversion_blocks)
synteny_joined <- inner_join(alignment_score, syntenic_blocks)
inv_joined <- as_tibble(inversion_joined) #converting data frame into tibble
syn_joined <- as_tibble(synteny_joined)
"Genus1_table" <- full_join(syn_joined, inv_joined)
  
##If you have done SyRI analysis on multiple genera, repeat above with different genus names and concatenate all tables in one master table: 
Master_data_table <- rbind(Genus1_table, Genus2_table ... )

Master_score_table <- Master_data_table %>%
  mutate(normalized_percent_ID = percent_ID*0.01*length_1) %>% 
  filter(length_1 > 1000) #filter out any regions shorter than 1000 bp to eliminate misleading smaller alignments and focus on the larger alignment regions. 
  
syn_plot <- Master_score_table %>%
  filter(type == "syn") %>%
  group_by(region_start_1, region_end_1, Genus) %>% 
  summarise(region_percent_ID = sum(normalized_percent_ID, na.rm = TRUE)/sum(length_1, na.rm = TRUE))
syn_plot <- full_join(syn_plot, divergence, by = "Genus")
divergence <- syn_plot %>% 
  group_by(Genus) %>%
  summarise(syn_divergence = mean(region_percent_ID))

inv_plot2 <- Master_score_table %>%
  filter(type == "inv") %>%
  group_by(region_start_1, region_end_1, Genus) %>% 
  summarise(region_percent_ID = sum(normalized_percent_ID, na.rm = TRUE)/sum(length_1, na.rm = TRUE)) 
inv_plot2 <- Master_score_table %>%
  filter(type == "inv") %>%
  mutate(region_length = region_end_1 - region_start_1)%>%
  group_by(region_end_1, region_start_1, region_length, Genus) %>%
  count()
inversion_stats2 <- inv_plot2 %>%
  group_by(Genus)%>%
  count() %>%
  rename("total_inv_regions_number" = n)
inversion_stats2 <- left_join(inversion_stats2, genome_length, by = "Genus")
inversion_stats_plot2 <- full_join(inversion_stats2, divergence, by = "Genus")

##Visualize on plot
syn_plot %>%
  ggplot() +
  geom_histogram(aes(x=region_percent_ID), binwidth = 0.002) +
  facet_wrap("Genus", nrow = 5) +
  geom_vline(aes(xintercept=syn_divergence))
ggsave("synteny_score_distribution_species_v1.pdf", width = 20, height = 16)

inv_plot2 %>% # grouped by inversion region 
  ggplot() +
  geom_histogram(aes(x=region_percent_ID), binwidth = 0.002) +
  facet_wrap("Genus", nrow = 5)
ggsave("inversion_score_distribution_species_v1.pdf", width = 20, height = 16)

syn_inv_plot <- rbind(syn_plot, inv_plot2)
syn_inv_plot %>%
  ggplot() +
  geom_density(aes(x=region_percent_ID, fill=type, alpha=0.4)) +
  facet_wrap("Genus", nrow = 5)
ggsave("synteny_inversion_score_distribution_32species.pdf", width = 20, height = 16)

