#gene count at inversion breakpoints
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)
getwd()
setwd("/Volumes/Backup Plus/Comparative genomics - inversion genomics/Scripts/Shell scripts/Comparative genomics")

#Data Import function
## Define directory with your data formatted as directory/genera/
genera <- c("Acer", "Actinidia", "Arabis", "Arachis",
            "Citrus", "Corylus", "Corymbia", "Cucumis", 
            "Eukalyptus", "Fragaria", "Glycine", "Gossypium", 
            "Ipomoea", "Juglans", "Malus", "Medicago", 
            "Phaseolus", "Populus", "Prunus", "Pyrus", 
            "Quercus", "Rhododendron", "Rosa", "Rubus", 
            "Salix", "Salvia", "Solanum", 
            "Vaccinium", "Vigna", "Vitis")

#Data Import function 
datalist_bedtools_count = list()

#### Gene #####
for (genus in genera){
  ## for inversion breakpoints
  bedtools_count_gene_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_count_gene_at_breakpoints.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "gene_hit"= V4) %>% 
    add_column(Genus = genus, type = "inv_breakpoint")
  ## for random points in the syntenic region of the genome
  bedtools_count_gene_syn_aligned_random <- read.table(paste0(directory, "/", genus, "/bedtools_count_genes_at_random_points_in_syn_region.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "gene_hit"= V4) %>% 
    add_column(Genus = genus, type = "syn_random_point")
  bedtools_count_all <- rbind(bedtools_count_gene_inv_breakpoints, bedtools_count_gene_syn_aligned_random)
  datalist_bedtools_count[[genus]] <- bedtools_count_all # add it to your list
  
}

Master_bedtools_count_gene = do.call(rbind, datalist_bedtools_count)

#### CDS ##### for testing whether this matches with proportion of CDS in syntenic region.
for (genus in genera){
  ## for random points in the syntenic region of the genome
  bedtools_count_cds_syn_aligned_random <- read.table(paste0(directory, "/", genus, "/bedtools_count_cds_at_random_points_in_syn_region.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "cds_hit"= V4) %>% 
    add_column(Genus = genus, type = "syn_random_point")
  #bedtools_count_all <- rbind(bedtools_count_gene_inv_breakpoints, bedtools_count_gene_aligned_random, bedtools_count_gene_syn_aligned_random)
  datalist_bedtools_count[[genus]] <- bedtools_count_cds_syn_aligned_random # add it to your list
  
}

Master_bedtools_count_cds_test = do.call(rbind, datalist_bedtools_count)

## Testing bedtools -count plots ##
Master_bedtools_count_gene %>%
  group_by(Genus, type)%>%
  summarise(
     n = n(),
     freq_of_gene = sum(gene_hit)/n
   ) %>%
  #filter(type =="inv_breakpoint")%>% 
  #filter(freq_of_gene < 0.2) %>%
  #filter(freq_of_gene > 0.5) %>%
  #filter(Genus == "Arabis"|Genus == "Arachis" |Genus == "Corymbia" | Genus == "Glycine"| Genus == "Malus"| Genus == "Pyrus" | Genus == "Rosa" | Genus == "Salvia")%>% 
  filter(Genus == "Cucumis"|Genus == "Ipomoea" |Genus == "Medicago" | Genus == "Phaseolus"| Genus == "Salix"| Genus == "Solanum"| Genus == "Vigna")%>% 
  ggplot()+
  geom_point(aes(x = type, y = freq_of_gene, colour = type))+
  geom_line(aes(x = type, y = freq_of_gene, group = Genus))
  #geom_boxplot(aes(x = type, y = freq_of_gene, colour = type))+
  #geom_jitter(aes(x = type, y = freq_of_gene, colour = type))+
  scale_color_manual(values = pal2)
  
  Master_bedtools_count_cds_test%>% 
    group_by(Genus)%>%
    summarise(
      n = n(),
      freq_CDS_in_syntenic_region = sum(cds_hit)/n
    ) %>%
    ggplot()+
    geom_jitter(aes(x=Genus, y=freq_CDS_in_syntenic_region))

  
##### Plotting #####
  Master_bedtools_count_gene <- read_csv("Master_bedtools_count_gene.csv")
  Master_bedtools_count_gene %>%
    filter(type == "inv_breakpoint") %>%
    group_by(Genus) %>%
    summarise(
      n = n(),
      `proportion of inversions spanning a gene` = sum(gene_hit)/n, 
      `% of inversions spanning a gene` = `proportion of inversions spanning a gene` * 100
    ) %>%
    arrange(., `% of inversions spanning a gene`)%>%
    ggplot()+
    scale_x_discrete(limits = c("Arachis", "Gossypium", "Corymbia", "Malus", "Salvia", "Rosa", "Pyrus", "Arabis", "Glycine", "Corylus", "Eukalyptus", "Vitis", "Quercus", "Citrus", "Rhododendron", "Acer", "Prunus", "Juglans", "Actinidia", "Fragaria", "Vaccinium", "Rubus", "Populus", "Phaseolus", "Salix", "Cucumis", "Solanum", "Medicago", "Ipomoea", "Vigna"))+
    geom_col(aes(x=Genus, y=`% of inversions spanning a gene`), fill = pal[10], colour = pal[10])+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

