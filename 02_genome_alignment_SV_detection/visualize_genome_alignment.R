library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

#The "table_original.txt" is the coords file you obtain from genome alignment file (sam), the output of minimap2 converted by samtocoords.py script. 
coords_file <- read.table("table_original.txt", head=TRUE) %>%
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
  select(-alignment) %>% 
  filter(length_1 > 3000)

pdf("coords_file_for_visualization_original.v0.pdf", height=24,width=24)
coords_file %>% 
  ggplot(., aes(x=block_start_1, y=block_start_2)) +
  facet_grid(chr_1 ~ chr_2) +
  geom_point(aes(colour = direction_8))
dev.off()
