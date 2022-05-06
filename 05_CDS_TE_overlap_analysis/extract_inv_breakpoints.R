library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

### Extract Inversion breakpoints from SyRI output ### 
Master_score_table_updated_inv <- Master_score_table_updated %>%
  filter(type =="inv")%>%
  filter(length_1 >= 5000) ## resulted in 0 inv from Vigna. 
end_2k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_2k = block_end_1 - 2000,
            upper_2k = block_end_1 + 2000)
start_2k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_2k = block_start_1 - 2000,
            upper_2k = block_start_1 + 2000)
breakpoints_regions <- rbind(end_2k, start_2k)
for (genus in genera){
  breakpoints_regions %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_regions_', genus, '.csv')) #Creates a separate brekapoints_regions.csv file for each genus you've analyzed
}
