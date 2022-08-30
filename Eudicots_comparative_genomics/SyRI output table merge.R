library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

#Data import: change the genus name as appropriate when you import tables. 
setwd("..../sub_folder/[genus_name]")
genera <- c("Acer", "Actinidia", "Arabis", "Arachis", 
                "Citrus", "Corylus", "Corymbia", "Cucumis", 
                "Eukalyptus", "Fragaria", 
                "Glycine", "Gossypium", "Ipomoea", 
                "Juglans", "Luffa", "Malus", "Medicago", 
                "Phaseolus", "Populus", "Prunus", "Pyrus", 
                "Quercus", 
                "Rhododendron", "Rosa", "Raphanus", "Rubus", 
                "Salix", "Salvia", "Solanum", 
                "Vaccinium", "Vigna", "Vitis") #32 genus pair

#for (genus in sub_folder){
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
  alignment_score$Genus <- "Vitis"
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
  "Vitis_table" <- full_join(syn_joined, inv_joined)
#}

#Concatenate all tables in one master table: 
Master_data_table <- rbind(Acer_table, Actinidia_table, Arabis_table, Arachis_table, 
                           Citrus_table, Corylus_table, Corymbia_table, Cucumis_table, 
                           Eukalyptus_table, Fragaria_table, 
                           Glycine_table, Gossypium_table, Ipomoea_table, 
                           Juglans_table, Luffa_table, Malus_table, Medicago_table, 
                           Phaseolus_table, Populus_table, Prunus_table, Pyrus_table, 
                           Quercus_table, 
                           Rhododendron_table, Rosa_table, Raphanus_table, Rubus_table, 
                           Salix_table, Salvia_table, Solanum_table, 
                           Vaccinium_table, Vigna_table, Vitis_table)
write.csv(Master_data_table, "Master_data_table_unfiltered_32_genera.csv")

#### Species divergence calculation ####
Master_score_table_updated <- Master_data_table %>% #Master_score_table_updated is UNfiltered by length or whatsoever.
  mutate(normalized_percent_ID = percent_ID*0.01*length_1) %>% 
  mutate(region_length = region_end_1 - region_start_1)

syn_plot <- Master_score_table_updated %>%
  filter(type == "syn") %>%
  filter(region_length >1000) %>%
  group_by(region_start_1, region_end_1, Genus, type) %>% 
  summarise(region_percent_ID = sum(normalized_percent_ID, na.rm = TRUE)/sum(length_1, na.rm = TRUE))
divergence <- syn_plot %>% 
  group_by(Genus) %>%
  summarise(syn_divergence = mean(region_percent_ID))
syn_plot <- full_join(syn_plot, divergence, by = "Genus")

inv_plot2 <- Master_score_table_updated %>%
  filter(type == "inv") %>%
  #filter(region_length >1000000) %>%
  group_by(region_start_1, region_end_1, Genus, type) %>% 
  summarise(region_percent_ID = sum(normalized_percent_ID, na.rm = TRUE)/sum(length_1, na.rm = TRUE)) 
inv_plot2 <- Master_score_table_updated %>%
  filter(type == "inv") %>%
  filter(region_length >1000)%>% 
  group_by(region_end_1, region_start_1, region_length, Genus, chr_1) %>%
  count()
inversion_stats2 <- inv_plot2 %>%
  group_by(Genus)%>%
  count() %>%
  rename("total_inv_regions_number" = n)

inversion_stats2 <- left_join(inversion_stats2, genome_length, by = "Genus")
inversion_stats2 <- left_join(inversion_stats2, biological_info)
inversion_stats_plot2 <- full_join(inversion_stats2, divergence, by = "Genus")

inversion_stats2_by_chr <- inv_plot2 %>%
  group_by(Genus, chr_1)%>%
  count() %>%
  rename("total_inv_regions_number" = n)
combined_fai <- read_table("chr_length_for_32_genera.txt", col_names = FALSE)
combined_fai_file <- combined_fai %>% rename("chr_1" = X1, "chr_length" = X2, "Genus" = X3)
plot_chr_length <- left_join(combined_fai_file, inversion_stats2_by_chr)

#### Visualize on plot ####
Master_score_table_updated %>% 
  filter(type == "inv") %>%
  group_by(Genus, chr_1, type, region_start_1, region_end_1, region_length) %>%
  distinct()%>%
  summarise(number_of_blocks_aligned = n()) %>%
  ggplot()+
  geom_density(aes(x=region_length, fill = type, alpha = 0.4))+
  scale_fill_manual(values = pnw_palette("Bay",1)) +
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_continuous(trans='log10')+
  xlab("Length of inversion")
  
Master_score_table_updated_inv <- Master_score_table_updated %>% 
  filter(type == "inv") %>%
  group_by(Genus, chr_1, type, region_start_1, region_end_1, region_length) %>%
  distinct()%>%
  summarise(number_of_blocks_aligned = n())
number_of_inv_by_length <- Master_score_table_updated_inv %>% 
  mutate(region_length_category = case_when(
    region_length < 1000 ~ "<1kbp",
    region_length >= 1000 & region_length <10000 ~ "1kbp-10kbp",
    region_length >= 10000 & region_length <100000 ~ "10kbp-100kbp", 
    region_length >= 100000 & region_length <1000000 ~ "100kbp-1Mbp", 
    region_length >= 1000000 ~ ">1Mbp"
  ))
number_of_inv_by_length$region_length_category <- factor(number_of_inv_by_length$region_length_category, levels = c("<1kbp", "1kbp-10kbp", "10kbp-100kbp", "100kbp-1Mbp", ">1Mbp"))
number_of_inv_by_length %>%
  group_by(region_length_category)%>%
  count()%>%
  mutate(proportion_of_inversions = n/6140)%>%
  ggplot(aes(x=region_length_category, y=proportion_of_inversions, fill = region_length_category))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = pal5) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab("")+ ylab("Proportion of identified inversions")+ labs(fill = "Inversion length") #THIS is all together.
number_of_inv_by_genus <- Master_score_table_updated_inv %>% 
  group_by(Genus) %>% 
  count() %>%
  rename("number_of_inversions" = n)

genome_length %>% 
  filter(ref_or_qry == "ref") %>% 
  mutate(proportion_of_genome_in_inversions = inv_region_length/Genome_length)%>% 
  arrange(., proportion_of_genome_in_inversions) %>%
  ggplot()+
  geom_col(aes(x=reorder(Genus, proportion_of_genome_in_inversions), y=proportion_of_genome_in_inversions, fill = ref_or_qry, alpha = 0.4))+
  xlab("Genus")+ ylab("Proportion of genome in inversions")+
  scale_fill_manual(values = pnw_palette("Bay",1)) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
Master_data_table_summary %>% filter(ref_or_qry == "ref") %>% 
  mutate(proportion_of_genome_in_inversions = inv_region_length/Genome_length)%>% 
  arrange(., number_of_inversions) %>%
  ggplot()+
  geom_col(aes(x=reorder(Genus, number_of_inversions), y=number_of_inversions, fill = ref_or_qry, alpha = 0.4))+
  xlab("Genus")+ ylab("Numbe of inversions")+
  scale_fill_manual(values = pnw_palette("Bay",1)) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
genome_divergence <- divergence %>% 
  transmute(Genus, percent_genome_divergence = (1-syn_divergence)*100)

syn_plot %>%
  ggplot() +
  geom_histogram(aes(x=region_percent_ID), binwidth = 0.002) +
  facet_wrap("Genus", nrow = 5) +
  geom_vline(aes(xintercept=syn_divergence))
ggsave("synteny_score_distribution_32species_v3.pdf", width = 20, height = 16)

inv_plot2 %>% # grouped by inversion region 
  ggplot() +
  geom_histogram(aes(x=region_percent_ID), binwidth = 0.002) +
  facet_wrap("Genus", nrow = 5)
ggsave("inversion_score_distribution_32species_v3.pdf", width = 20, height = 16)

syn_inv_plot <- read.csv("syn_inv_plot_updated.csv")
syn_inv_plot %>%
  ggplot() +
  geom_density(aes(x=region_percent_ID, fill=type, alpha=0.4)) +
  scale_fill_manual(values = pnw_palette("Bay",2)) +
  facet_wrap("Genus", nrow = 5)
ggsave("synteny_inversion_score_distribution_32species_v3.pdf", width = 20, height = 16)

genomic_TE_proportion <- Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(genomic_TE_proportion = total_TE_length / Genome_length)%>%
  ggplot(., aes(x=reorder(Genus,-genomic_TE_proportion), y=genomic_TE_proportion))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +xlab("Genus")

Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  ggplot()+
  geom_point(aes(x=percent_genome_divergence, y=number_of_inversions))+
  geom_smooth(aes(x=percent_genome_divergence, y=number_of_inversions, fill = TRUE), method="lm", colour = pal[25])+
  theme_classic()+
  theme(legend.position = "na")+
  scale_fill_manual(values = pal[25])+
  xlab("Genome sequence divergence (%)")+ ylab("Number of inversions")
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  ggplot()+
  geom_point(aes(x=Genome_length, y=number_of_inversions))+
  geom_smooth(aes(x=Genome_length, y=number_of_inversions), method="lm")+
  xlab("Genome length (bp)")+ ylab("Number of inversions")
Master_data_table_summary %>%
  filter(ref_or_qry == "ref") %>%
  mutate(TE_proportion = total_TE_length/Genome_length*100)%>% 
  ggplot()+
  geom_point(aes(x=TE_proportion, y=number_of_inversions))+
  geom_smooth(aes(x=TE_proportion, y=number_of_inversions), method="lm")+
  xlab("Genomic TE proportion (%)")+ ylab("Number of inversions")

plot_chr_length %>% 
  ggplot()+
  geom_point(aes(x=chr_length, y=total_inv_regions_number))+ 
  geom_smooth(aes(x=chr_length, y=total_inv_regions_number), method ="lm")+
  xlab("Chromosome length (bp)")+ ylab("Number of inversions")

number_of_inv_by_length%>% #rate by inversion size categories
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
  
#example supplementary plot of box plot with jitter: 
assembly_method %>%
  left_join(., Master_data_table_summary)%>%
  mutate(proportion_of_inverted_genome = inv_region_length/Genome_length) %>%
  ggplot()+
  geom_boxplot(aes(x=Sequencing_platform, y=number_of_inversions, fill = TRUE, alpha = 0.4), outlier.shape = NA)+
  geom_jitter(aes(x=Sequencing_platform, y=number_of_inversions))+
  theme_classic()+
  scale_fill_manual(values = pal[25])+
  scale_colour_manual(values = pal[25])+
  ylab("Number of inversions")

