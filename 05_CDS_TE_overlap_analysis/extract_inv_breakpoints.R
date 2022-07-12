library(tidyverse)
library(dplyr)

### Extract Inversion breakpoints from SyRI output ### 
Master_score_table_updated_inv <- Master_score_table_updated %>%
  filter(type =="inv")%>%
  filter(length_1 >= 5000)
#4kb window: 
end_4k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_end = region_end_1, 
            upper_end = region_end_1 + 4000)
start_4k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_end = region_start_1 -4000, 
            upper_end = region_start_1)
breakpoints_regions <- rbind(end_4k, start_4k)
for (genus in genera){
  breakpoints_regions %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_regions_', genus, '.csv')) #Creates a separate brekapoints_regions.csv file for each genus you've analyzed
}
#10kb window: 
end_10k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_end = region_end_1, 
            upper_end = region_end_1 + 10000)
start_10k <- Master_score_table_updated_inv%>%
  transmute(Genus, chr_1,
            lower_end = region_start_1 -10000, 
            upper_end = region_start_1)
breakpoints_regions_onesided_10k <- rbind(end_10k, start_10k)
for (genus in genera){
  breakpoints_regions_onesided_10k %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_regions_onesided_10k_', genus, '.csv'))
}
