# Specific measurements
percentage_computation_for_counting <- function(dataset){
    dataset %>%
        group_by(BCAC_ID, marker) %>% summarise(percentage = mean(perAll)) %>% distinct() %>%
        mutate(percentage = percentage / 100,
               percentage = case_when(percentage == 0 ~ percentage + 1e-5, TRUE ~ percentage),
               percentage = case_when(percentage == 1 ~ percentage - 1e-5, TRUE ~ percentage)) %>% 
        left_join(., phenotype_BRIDGES, by = "BCAC_ID") %>% 
        filter(!is.na(BRIDGES_ID))
}


# Number of cases per marker (all tumours)
TILs_scores_dataset_BRIDGES %>% 
  percentage_computation_for_counting(.) %>%
  group_by(marker) %>%
  summarise(count = n_distinct(BCAC_ID))


# Number of TNBC cases per marker
Number of distinct BCAC_IDs 
TILs_scores_dataset_BRIDGES %>% 
  percentage_computation_for_counting(.) %>%  
  filter(TNBC == 1) %>%
  group_by(marker) %>%
  summarise(count = n_distinct(BCAC_ID))


## Number of cases per study
TILs_scores_dataset_BRIDGES %>% 
  percentage_computation_for_counting(.) %>%
  group_by(study, marker) %>%
  summarise(
    count_patients = n(),
    age = paste0(median(AgeDiagIndex), " ", "(", min(AgeDiagIndex), " - ", max(AgeDiagIndex), ")")
  ) %>%
  arrange(desc(marker), study)
