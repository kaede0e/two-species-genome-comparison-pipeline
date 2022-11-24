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

#1kb x10 window: 
  end_1k_1 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1,
              end = region_end_1 + 1000)
  end_1k_2 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 1000,
              end = region_end_1 + 2000)
  end_1k_3 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 2000,
              end = region_end_1 + 3000)
  end_1k_4 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 3000,
              end = region_end_1 + 4000)
  end_1k_5 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 4000,
              end = region_end_1 + 5000)
  end_1k_6 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 5000,
              end = region_end_1 + 6000)
  end_1k_7<- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 6000,
              end = region_end_1 + 7000)
  end_1k_8 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 7000,
              end = region_end_1 + 8000)
  end_1k_9 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 8000,
              end = region_end_1 + 9000)
  end_1k_10 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              start = region_end_1 + 9000,
              end = region_end_1 + 10000)
  end_1k_coordinates <- rbind(end_1k_1, end_1k_2, end_1k_3, end_1k_4, end_1k_5, end_1k_6, end_1k_7, end_1k_8, end_1k_9, end_1k_10)

  start_1k_1 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 1000,
              upper_2k = region_start_1)
  start_1k_2 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 2000,
              upper_2k = region_start_1 - 1000)
  start_1k_3 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 3000,
              upper_2k = region_start_1 - 2000)
  start_1k_4 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 4000,
              upper_2k = region_start_1 - 3000)
  start_1k_5 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 5000,
              upper_2k = region_start_1 - 4000)
  start_1k_6 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 6000,
              upper_2k = region_start_1 - 5000)
  start_1k_7 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 7000,
              upper_2k = region_start_1 - 6000)
  start_1k_8 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 8000,
              upper_2k = region_start_1 - 7000)
  start_1k_9 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 9000,
              upper_2k = region_start_1 - 8000)
  start_1k_10 <- Master_score_table_updated_inv%>%
    transmute(Genus, chr_1,
              lower_2k = region_start_1 - 10000,
              upper_2k = region_start_1 - 9000)
  start_1k_coordinates <- rbind(start_1k_1, start_1k_2, start_1k_3, start_1k_4, start_1k_5, start_1k_6, start_1k_7, start_1k_8, start_1k_9, start_1k_10)
  
for (genus in genera){
  end_1k_coordinates %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_end_1k_coordinates_', genus, '.csv'))
}
for (genus in genera){
  start_1k_coordinates %>%
    filter(Genus == genus) %>%
    write.csv(., paste0('breakpoints_start_1k_coordinates_', genus, '.csv'))
}
