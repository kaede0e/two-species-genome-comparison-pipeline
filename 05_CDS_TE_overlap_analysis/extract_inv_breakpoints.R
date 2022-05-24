library(tidyverse)

### Extract Inversion breakpoints from SyRI output ### 
Master_score_table_updated_inv <- Master_score_table_updated %>%
  filter(type =="inv")%>%
  filter(length_1 >= 5000)
end_4k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_2k = block_end_1 - 2000,
            upper_2k = block_end_1 + 2000)
start_4k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_2k = block_start_1 - 2000,
            upper_2k = block_start_1 + 2000)
breakpoints_regions <- rbind(end_4k, start_4k)
for (genus in genera){
  breakpoints_regions %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_regions_', genus, '.csv')) #Creates a separate brekapoints_regions.csv file for each genus you've analyzed
}
