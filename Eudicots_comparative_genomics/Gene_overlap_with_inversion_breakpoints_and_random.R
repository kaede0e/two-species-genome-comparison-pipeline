#gene count at inversion breakpoints / inversion breakpoints count at any genes
library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

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
  ## how many inversion breakpoints intersect with any genes
  bedtools_count_gene_inv_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_count_gene_at_breakpoints.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "gene_hit"= V4) %>% 
    add_column(Genus = genus, type = "inv_breakpoint")
  ## how many random points in the syntenic region of the genome intersect with any genes
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

#### Inversion breakpoints #####
for (genus in genera){
  ## how many genes intersect with inversion breakpoints
  bedtools_count_inv_gene_breakpoints <- read.table(paste0(directory, "/", genus, "/bedtools_count_inversion_breakpoints_at_genes.txt"), head = FALSE) %>%
    rename(chr_1 = V1, 
           region_start_1 = V2, 
           region_end_1 = V3, 
           "inv_hit"= V4) %>% 
    add_column(Genus = genus, type = "gene_inv_overlap")
  datalist_bedtools_count[[genus]] <- bedtools_count_inv_gene_breakpoints # add it to your list
}

Master_bedtools_count_gene2 = do.call(rbind, datalist_bedtools_count)

  
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
    geom_col(aes(x=Genus, y=`% of inversions spanning a gene`), fill = pal[33], colour = NA)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  Master_bedtools_count_gene2 <- read_csv("Master_bedtools_count_gene2.csv")
  Master_bedtools_count_gene2 %>%
    rename(inv_hits = inv_hit) %>% 
    mutate(inv_hit = case_when(
      inv_hits == 0 ~ 0, 
      inv_hits > 0 ~ 1 #replace any double hits with '1' so anything with at least one inversion overlap is counted properly.
    )) %>%
    group_by(Genus) %>%
    summarise(
      n = n(), 
      `proportion of genes overlapping with breakpoints` = sum(inv_hit)/n,
      `% of genes overlapping with breakpoints` = `proportion of genes overlapping with breakpoints` * 100
    )%>%
    arrange(., `% of genes overlapping with breakpoints`)%>%
    ggplot()+
    scale_x_discrete(limits = c("Arachis", "Pyrus", "Glycine", "Gossypium", "Eukalyptus", "Malus", "Arabis", "Rosa", "Corymbia", "Salvia", "Citrus", "Vitis", "Rubus", "Fragaria", "Juglans", "Vigna", "Quercus", "Prunus", "Populus", "Solanum", "Phaseolus", "Ipomoea", "Acer", "Cucumis", "Actinidia", "Medicago", "Corylus", "Rhododendron", "Salix", "Vaccinium"))+
    geom_col(aes(x=Genus, y=`% of genes overlapping with breakpoints`), fill = pal[33], colour = NA)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
