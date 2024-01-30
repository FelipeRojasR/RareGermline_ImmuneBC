##############################
# Dependencies for project
list.of.packages <- c(
  "dplyr", "betareg", "DBI", "ggplot2", "foreach", "readr", "stringr",
  "doParallel", "ranger", "palmerpenguins", "tidyr",
  "kableExtra", "rbgen", "ggpubr", "ggplot2", "forestploter",
  "ggpmisc", "ggrepel", "car", "boot", "patchwork", 
  "gt", "ggpmisc", "rlang", "forcats")

`%!in%` = Negate(`%in%`)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dependencies = TRUE)
}

# loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}


##############################
# Importing data 
# Phenotype data
# Phenotype BRIDGES includes the a new column (binary) for samples being triple negative
phenotype_BRIDGES <- read_tsv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/concept_734_schmidt_euro_cases_bridges_pheno_v15_02.txt") %>% 
  filter(ethnicityClass == 1) %>%
  select(BRIDGES_ID, BCAC_ID, ER_statusIndex, AgeDiagIndex, HER2_statusIndex, PR_statusIndex, study) %>%
  mutate(ER_statusIndex = case_when(ER_statusIndex == 888 ~ NA_real_, TRUE ~ ER_statusIndex), 
         HER2_statusIndex = case_when(HER2_statusIndex == 888 ~ NA_real_, TRUE ~ HER2_statusIndex), 
         PR_statusIndex = case_when(PR_statusIndex == 888 ~ NA_real_, TRUE ~ PR_statusIndex), 
         AgeDiagIndex = case_when(AgeDiagIndex == 888 ~ NA_real_, TRUE ~ AgeDiagIndex)) %>%
  mutate(TNBC = case_when(ER_statusIndex == 0 & HER2_statusIndex == 0 & PR_statusIndex == 0 ~ 1, TRUE ~ 0))


# Link between BCAC database
LINK_WITH_BCAC <- read_csv(file = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/Official_Data_Aaron/FelipeFinal_link_BCAC_ID_221128.csv")


# Database of scores obtained from Aaron.
TILs_scores_dataset_BRIDGES <- read_csv(file = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/Official_Data_Aaron/Felipe_ImmuneScores_221128.csv") %>% 
  rename(BCAC_Image_Core_ID = imgShort) %>% 
  full_join(., LINK_WITH_BCAC, by = "BCAC_Image_Core_ID") %>% 
  relocate(BCAC_ID, .before = `Image Tag`) %>%
  select(BCAC_ID, `Analysis Region`, marker, artifactPercent, 
         countAll, countStroma, countTumour, 
         countRazaAll, countRazaTumour, countRazaStroma, 
         perAll, perTumour, perStroma) %>% 
  filter(artifactPercent < 0.4 & `Analysis Region` == "entire image")


##############################
# This function computes the average percentage of immune cells by marker by BCAC-ID
# Adds 1e-5 to any score equal to 0.
# Substracts 1e-5 to any score  equal to 1.
# This was done to ensure that the beta regression can be used.

# This function does not group per marker (group by line) because it expects to be 
# used in an object that only contains a single marker.
percentage_computation <- function(dataset){
  # Adjustment required for the beta regression.
  dataset %>%
    group_by(BCAC_ID) %>% summarise(percentage = mean(perAll)) %>% distinct() %>%
    mutate(percentage = percentage / 100,
           percentage = case_when(percentage == 0 ~ percentage + 0.00001, TRUE ~ percentage),
           percentage = case_when(percentage == 1 ~ percentage - 0.00001, TRUE ~ percentage))
}


##############################
# Immune markers included
markers <- c("CD8", "CD20", "FOXP3", "CD163")


##############################
# All genes and respective functions related to DNA damage repair mechanims. This names MUST match the names created in the germline objects.
gene_and_functions <- tibble(gene = c("ABRAXAS1", "AKT1", "ATM", "BABAM2", "BARD1", "BRCA1", "BRCA2", "BRIP1", 
                                      "CDH1", "CHEK2", "CHEK2_1100delC", "CHEK2.other", 
                                      "EPCAM", "FANCC", "FANCM", "GEN1", "MEN1", "MLH1", "MRE11", 
                                      "MSH2", "MSH6", "MUTYH", "NBN", "NF1", "PALB2", "PIK3CA", "PMS2", "PTEN",
                                      "RAD50", "RAD51C", "RAD51D", "RECQL", "RINT1", "STK11", "TP53", "XRCC2"), 
                             mechanism =c("DSB", "nonDDR", "Associated", "DSB", "DSB", "DSB", "DSB", "DSB", 
                                          "nonDDR", "Associated", "Associated", "Associated",
                                          "nonDDR", "DSB", "DSB", "DSB", "nonDDR", "SSB", "DSB", "SSB",
                                          "SSB", "SSB", "DSB", "nonDDR", "DSB", "nonDDR", "SSB", "nonDDR", "DSB", 
                                          "DSB", "DSB", "DSB", "DSB", "nonDDR", "Associated", "DSB")) %>%
  mutate(mechanism = factor(mechanism, levels = c("DSB", "SSB", "Associated", "nonDDR")))


##############################
# This dataframe contains all scores obtained across the IHC datasets.
TilScores_BridgesPhenotype_markers <- bind_rows(
  lapply(markers, function(markers){
    TILs_scores_dataset_BRIDGES %>% filter(marker == markers) %>% 
      percentage_computation(.) %>% 
      left_join(., phenotype_BRIDGES, by = "BCAC_ID") %>%
      filter(!is.na(BRIDGES_ID)) %>%
      mutate(marker = markers)
  }
  )
)


##############################
# Function to combine MSVs and PTVs with the BRIDGES dataset
merge_with_BRIDGES <- function(dataset, variants_dataset){
  dataset %>% 
    select(percentage, BRIDGES_ID, ER_statusIndex, AgeDiagIndex, study) %>%
    inner_join(., variants_dataset, by = "BRIDGES_ID")
}

create_germline_with_phenotype <- function(phenotype_dataset, germline_variants_dataset){
  dataset_output <- bind_rows(
    lapply(markers, function(marker_input) {
      merge_with_BRIDGES(
        phenotype_dataset %>% filter(marker == as.character(marker_input)),
        germline_variants_dataset) %>%
        mutate(marker = marker_input) %>%
        relocate(marker, .before = percentage)
    }))
  
  return(dataset_output)
}


##############################
# This creates all 
# All cancer cases
missense_merged_results <- create_germline_with_phenotype(phenotype_dataset = TilScores_BridgesPhenotype_markers,
                                                          germline_variants_dataset = pathogenic_missense_bridges)

truncating_merged_results <- create_germline_with_phenotype(phenotype_dataset = TilScores_BridgesPhenotype_markers, 
                                                            germline_variants_dataset = truncating_bridges.including.chek2mods)

# Triple negative breast cancer cases
TNBC_missense_merged_results <- create_germline_with_phenotype(phenotype_dataset = TilScores_BridgesPhenotype_markers %>% filter(TNBC == 1),
                                                               germline_variants_dataset = pathogenic_missense_bridges)

TNBC_truncating_merged_results <- create_germline_with_phenotype(phenotype_dataset = TilScores_BridgesPhenotype_markers %>% filter(TNBC == 1), 
                                                                 germline_variants_dataset = truncating_bridges.including.chek2mods)

## REGRESSION MODEL
regression_across_datasets <- function(germline_dataset_combined, adjustment_input = NULL, carriers_for_sensitivity = NULL){
  if(!is.null(carriers_for_sensitivity)){
    bind_rows(
      lapply(markers, function(markers) {
        BetaRegression_RareVariants.sensitivity(dataset = germline_dataset_combined %>% filter(marker == markers), 
                                                carrier_dataset = carriers_for_sensitivity %>% filter(marker == markers)) %>%
          transform_output_function(.) %>%
          mutate(marker = markers)
      }
      )
    )
  }
  
  else {
    bind_rows(
      lapply(markers, function(markers) {
        BetaRegression_RareVariants(dataset = germline_dataset_combined %>% filter(marker == markers), 
                                    adjustment = as.character(adjustment_input)) %>%
          transform_output_function(.) %>%
          mutate(marker = markers)
      }
      )
    )
  }
}

# Baseline models
missense_baseline_results <- regression_across_datasets(missense_merged_results, adjustment = "baseline")
truncating_baseline_results <- regression_across_datasets(truncating_merged_results, adjustment = "baseline")

# Adjusted models
missense_adjusted_results <- regression_across_datasets(missense_merged_results, adjustment = "adjusted")
truncating_adjusted_results <- regression_across_datasets(truncating_merged_results, adjustment = "adjusted")

# Baseline TNBC models 
missense_TNBC_results <- regression_across_datasets(TNBC_missense_merged_results, adjustment = "baseline")
truncating_TNBC_results <- regression_across_datasets(TNBC_truncating_merged_results, adjustment = "baseline")




# FULL PRIMARY RESULTS
full_primary_results <- list(
  missense_baseline_results %>% mutate(variants = "missense", model = "baseline", cases = "overall"),
  truncating_baseline_results %>% mutate(variants = "truncating", model = "baseline", cases = "overall"),
  
  missense_adjusted_results %>% mutate(variants = "missense", model = "adjusted", cases = "overall"), 
  truncating_adjusted_results %>% mutate(variants = "truncating", model = "adjusted", cases = "overall"),
  
  missense_TNBC_results %>% mutate(variants = "missense", model = "baseline", cases = "TNBC"), 
  truncating_TNBC_results %>% mutate(variants = "truncating", model = "baseline", cases = "TNBC")
) %>%
  do.call(rbind.data.frame, .) %>%
  group_by(marker, variants, model, cases) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  left_join(., gene_and_functions, by = "gene")

full_primary_results[full_primary_results == 888] <- NA_real_


# Baseline Figure
figure_paper.baseline <- full_primary_results %>%
  filter(model == "baseline" & cases == "overall") %>%
  mutate(ci.lb = case_when(ci.lb <= 0.3 ~ 0.3, TRUE ~ ci.lb),
         ci.ub = case_when(ci.ub >= 5 ~ 5, TRUE ~ ci.ub),
         gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene)), 
         variants = factor(variants, levels = c("missense", "truncating")),
         sig.pval = TRUE,
         # sig.pval = case_when(p.value <= 0.05 ~ TRUE, TRUE ~ FALSE),
         marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163"))) %>% 
  ggplot(mapping = aes(x = estimate, y = fct_rev(gene), color = fct_rev(marker))) +
  theme_classic() + 
  geom_point(position = position_dodge2(width = 0.7)) + 
  facet_grid(cols = vars(variants), scales = "free") + 
  geom_linerange(aes(xmin = ci.lb, xmax = ci.ub), position = position_dodge2(width = 0.7)) + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  ####
  geom_vline(xintercept = 0.3, linetype = "solid") + 
  geom_vline(xintercept = 5, linetype = "solid") + 
  ####
  labs(x = "Odds Ratio (95% CI)", y = "") + 
  guides(color = guide_legend(reverse = TRUE)) + 
  coord_trans(x = 'log10') + 
  scale_color_manual(values=c("#E99999", "#C80000", "#80A2C4",  "#01458A")) +
  scale_x_continuous(breaks = c(0.5,1,2,5,10,20)) + 
  theme(legend.position = "top")

figure_paper.baseline



carriers.fix.data.function <- function(dataset, mechanisms = NULL){
  if(mechanisms == TRUE){
    dataset %>%
      select(-c(percentage, ER_statusIndex, AgeDiagIndex, study, marker)) %>% 
      pivot_longer(-BRIDGES_ID, names_to = "gene", values_to = "carrier") %>%
      group_by(gene) %>%
      summarise(carr = sum(carrier == 1),
                non.carr = sum(carrier == 0)) %>%
      arrange(gene)
    
  } else {
    dataset %>%
      select(-c(percentage, ER_statusIndex, AgeDiagIndex, study, marker)) %>% 
      pivot_longer(-BRIDGES_ID, names_to = "gene", values_to = "carrier") %>%
      mutate(gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene))) %>%
      group_by(gene) %>%
      summarise(carr = sum(carrier == 1),
                non.carr = sum(carrier == 0)) %>%
      arrange(gene)
  }
  
}


fix_data_estimate_table <- function(dataset){
  dataset %>% 
    mutate(gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene)), 
           variants = factor(variants, levels = c("missense", "truncating")),
           marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163")),
           mechanism = factor(mechanism, levels = c("DSB", "SSB", "Associated", "nonDDR"))
    ) %>% arrange(mechanism, gene) %>%
    ungroup() %>%
    transmute(marker, mechanism, gene, estimate,
              CI = paste0(ci.lb, " - ", ci.ub),
              p.value)
}


joined_all <- full_join(full_primary_results %>%
                          filter(model == "baseline" & cases == "overall") %>%
                          select(-c(ci.lb, ci.ub, p.value, p.adj, cases, model)),
                        
                        full_primary_results %>%
                          filter(model == "adjusted" & cases == "overall") %>%
                          select(-c(ci.lb, ci.ub, p.value, p.adj, cases, model)),
                        
                        by = c("gene", "marker", "variants", "mechanism"), 
                        suffix = c(".baseline", ".adjusted"))




# OVERALL COLUMNS 
# INPUT: variants, marker, cases
# raw_file_variants
compute_supplementary_tables <- function(variants_input, marker_input, cases_input, raw_file_variants_input){
  if(cases_input == "overall") {
    baseline_columns_tmp <- full_primary_results %>%
      filter(variants == as.character(variants_input) & marker == as.character(marker_input) & cases == as.character(cases_input) & model == "baseline") %>%
      fix_data_estimate_table(.)
    
    adjusted_columns_tmp <- full_primary_results %>%
      filter(variants == as.character(variants_input) & marker == as.character(marker_input) & cases == as.character(cases_input) & model == "adjusted") %>%
      fix_data_estimate_table(.)
    
    # Number of carriers
    carriers_tmp <- carriers.fix.data.function(raw_file_variants_input %>% filter(marker == as.character(marker_input))) # gene, carr, non.carr
    
    # Comparison of estimates
    comparison_of_estimates_tmp <- joined_all %>%
      filter(variants == as.character(variants_input) & marker == as.character(marker_input)) %>%
      # How large is the new one compared to the old one. 
      mutate(mediation.val = round(abs((estimate.baseline - estimate.adjusted) / estimate.baseline)* 100, digits = 1)) %>%
      mutate(gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene))) %>%
      select(mechanism, gene, mediation.val) %>% arrange(gene)
    
    inner_join(baseline_columns_tmp, adjusted_columns_tmp, by = c("mechanism", "gene", "marker"), suffix = c(".baseline", ".adjusted")) %>%
      inner_join(., comparison_of_estimates_tmp, by = c("mechanism", "gene")) %>%
      inner_join(., carriers_tmp, by = c("gene")) %>%
      dplyr::relocate(marker, .after = last_col()) %>%
      dplyr::mutate(variants = as.character(variants_input),
                    cases = as.character(cases_input))
    
  } else if(cases_input == "TNBC"){
    baseline_columns_tmp <- full_primary_results %>%
      filter(variants == as.character(variants_input) & marker == as.character(marker_input) & cases == as.character(cases_input) & model == "baseline") %>%
      fix_data_estimate_table(.)
    
    # Number of carriers
    carriers_tmp <- carriers.fix.data.function(raw_file_variants_input %>% filter(marker == as.character(marker_input)) %>%
                                                 filter(BRIDGES_ID %in% phenotype_BRIDGES %>% filter(TNBC == 1))) # gene, carr, non.carr
    
    inner_join(baseline_columns_tmp, carriers_tmp, by = c("gene")) %>%
      dplyr::relocate(marker, .after = last_col()) %>%
      dplyr::mutate(variants = as.character(variants_input),
                    cases = as.character(cases_input))
  }
}

# CREATE TABLE
compute_supplementary_tables(variants_input = "truncating", 
                             marker_input = "CD8", 
                             cases_input = "TNBC", 
                             raw_file_variants_input = TNBC_truncating_merged_results)





# FIGURE - BASELINE - Figure 2
baseline.tmp <- full_primary_results %>%
  filter(cases == "overall" & model == "baseline" & p.value < 0.05) %>%
  select(gene, marker, variants, mechanism,
         estimate, ci.lb, ci.ub, p.value) %>%
  mutate(model = "unadjusted")

adjusted.tmp <- full_primary_results %>%
  filter(cases == "overall" & model == "adjusted") %>%
  inner_join(baseline.tmp, ., by = c("gene", "marker", "variants", "mechanism"), suffix = c(".baseline", ".adjusted")) %>%
  select(gene, marker, variants, mechanism,
         estimate.adjusted, ci.lb.adjusted, ci.ub.adjusted, p.value.adjusted) %>%
  rename(estimate = estimate.adjusted, 
         ci.lb = ci.lb.adjusted, 
         ci.ub = ci.ub.adjusted, 
         p.value = p.value.adjusted) %>%
  mutate(model = "adjusted")

adjusted.joined.rbind <- rbind(baseline.tmp, adjusted.tmp)


































# MECHANISMS OF ACTION
mechanisms_dataset_creation <- function(dataset){
  dataset %>%
    pivot_longer(-BRIDGES_ID, names_to = "gene", values_to = "carrier") %>% 
    left_join(., gene_and_functions, by = "gene") %>%
    group_by(BRIDGES_ID, mechanism) %>%
    summarise(sum.mec = sum(carrier)) %>%
    mutate(mechanism.carrier = case_when(sum.mec >= 1 ~ 1, TRUE ~ 0)) %>%
    select(-sum.mec) %>%
    tidyr::spread(-mechanism, mechanism.carrier, fill = 0) %>% ungroup()
}


# MECHANISMS OF ACTION
missense_mechanisms <- mechanisms_dataset_creation(pathogenic_missense_bridges)
truncating_mechanisms <- mechanisms_dataset_creation(truncating_bridges.including.chek2mods)


missense_merged_results_mechanisms <- create_germline_with_phenotype(phenotype_dataset = TilScores_BridgesPhenotype_markers,
                                                                     germline_variants_dataset = missense_mechanisms)

truncating_merged_results_mechanisms <- create_germline_with_phenotype(phenotype_dataset = TilScores_BridgesPhenotype_markers, 
                                                                       germline_variants_dataset = truncating_mechanisms)


# Baseline models
missense_baseline_results_mechanisms <- regression_across_datasets(missense_merged_results_mechanisms, adjustment = "baseline")
truncating_baseline_results_mechanisms <- regression_across_datasets(truncating_merged_results_mechanisms, adjustment = "baseline")

full_mechanisms_results <- list(
  missense_baseline_results_mechanisms %>% mutate(variants = "missense", model = "baseline", cases = "overall"),
  truncating_baseline_results_mechanisms %>% mutate(variants = "truncating", model = "baseline", cases = "overall")
) %>%
  do.call(rbind.data.frame, .) %>%
  group_by(marker, variants, model, cases) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

full_mechanisms_results[full_mechanisms_results == 888] <- NA_real_


# Carriers are missing (Supplementary tables)
carriers_mechanisms <- function(germline_dataset_combined){
  bind_rows(
    lapply(markers, function(markers) {
      carriers.fix.data.function(dataset = germline_dataset_combined %>% filter(marker == as.character(markers)), mechanisms = TRUE) %>%
        mutate(marker = markers)
    }
    )
  )
}

carriers_missense_mechanisms <- carriers_mechanisms(missense_merged_results_mechanisms) 
carriers_truncating_mechanisms <- carriers_mechanisms(truncating_merged_results_mechanisms) 








any_carriers.function <- function(dataset.msv, dataset.ptv){
  any_carrier <- rbind(carriers.msv <- dataset.msv %>% 
                         mutate(sum.genes.tmp = rowSums(select(., -c(marker, percentage, ER_statusIndex, AgeDiagIndex, study, BRIDGES_ID)))) %>%
                         filter(sum.genes.tmp != 0) %>%
                         select(c(BRIDGES_ID)), 
                       
                       carriers_CD8.ptv <- dataset.ptv %>% 
                         mutate(sum.genes.tmp = rowSums(select(., -c(marker, percentage, ER_statusIndex, AgeDiagIndex, study, BRIDGES_ID)))) %>%
                         filter(sum.genes.tmp != 0) %>%
                         select(c(BRIDGES_ID))) %>%
    distinct(BRIDGES_ID)
  
  return(any_carrier)
}


# These are the IDs for which at least one germline variants was found
IDS_with_at_least_one_germlinevariant <- bind_rows(
  lapply(markers, function(markers){
    any_carriers.function(missense_merged_results %>% filter(marker == as.character(markers)),
                          truncating_merged_results %>% filter(marker == as.character(markers))) %>%
      mutate(marker = markers)
    }))

missense_sensitivity_results <- regression_across_datasets(germline_dataset_combined = missense_merged_results, carriers_for_sensitivity = IDS_with_at_least_one_germlinevariant)
truncating_sensitivity_results <- regression_across_datasets(germline_dataset_combined = truncating_merged_results, carriers_for_sensitivity = IDS_with_at_least_one_germlinevariant)


full_sensitivity_results <- list(
  missense_sensitivity_results %>% mutate(variants = "missense"),
  truncating_sensitivity_results %>% mutate(variants = "truncating")
) %>%
  do.call(rbind.data.frame, .) %>%
  group_by(marker, variants) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  left_join(., gene_and_functions, by = "gene")

full_sensitivity_results[full_sensitivity_results == 888] <- NA_real_



carriers_sensitivity <- function(dataset, carrier_dataset){
  columns_ignore.tmp <- c("marker", "percentage", "BRIDGES_ID", "ER_statusIndex", "AgeDiagIndex", "study", "to.be.included")
  
  matrix_to_fill.carrier.sens <- matrix(nrow = length(which(!names(dataset) %in% columns_ignore.tmp)), 
                                        ncol = 3)
  
  colnames(matrix_to_fill.carrier.sens) <- c("gene", "carrier", "non_carrier")
  
  dataset.t <- dataset %>%
    mutate(ER_statusIndex = case_when(ER_statusIndex == 888 ~ NA_real_, TRUE ~ ER_statusIndex))
  
  for(i in which(!names(dataset.t) %in% columns_ignore.tmp)){
    patterns <- exprs(!!sym(colnames(dataset.t[i])) == 1 & BRIDGES_ID %in% carrier_dataset$BRIDGES_ID  ~ "carrier",
                      !!sym(colnames(dataset.t[i])) == 0 & BRIDGES_ID %in% carrier_dataset$BRIDGES_ID  ~ "eliminate",
                      !!sym(colnames(dataset.t[i])) == 0 & BRIDGES_ID %!in% carrier_dataset$BRIDGES_ID  ~ "non_carrier")
    
    datas <- dataset.t  %>%
      mutate(to.be.included = case_when(!!!patterns)) %>% # "non_carrier"
      filter(to.be.included != "eliminate") %>% # Here we filter out all the non-carriers that have a variant in any other gene
      group_by(to.be.included) %>%
      summarise(n = n()) %>%
      mutate(gene = colnames(dataset.t[i])) %>%
      pivot_wider(names_from = "to.be.included", values_from = "n")
    
    if(length(colnames(datas)) < 3){
      datas <- datas %>%
        mutate(carrier = 0) %>%
        select("gene", "carrier", "non_carrier")
    } else{
      datas
    }
    
    matrix_to_fill.carrier.sens[seq_along(which(!names(dataset) %in% columns_ignore.tmp))[i-6],] <- as.character(datas[1,])
  }
  
  carrier.sens <- as_tibble(matrix_to_fill.carrier.sens) %>%
    mutate(carrier = as.numeric(carrier),
           non_carrier = as.numeric(non_carrier))
  
  return(carrier.sens)
}



supplementary_tables_sensitivity <- function(raw_file_variants_input, marker_input, variants_input){
  estimates_tmp <- full_sensitivity_results %>% 
    filter(marker == as.character(marker_input) & variants == as.character(variants_input)) %>%
    fix_data_estimate_table(.)
  
  carriers_tmp <- carriers_sensitivity(raw_file_variants_input %>% filter(marker == as.character(marker_input)), 
                                       IDS_with_at_least_one_germlinevariant %>% filter(marker == as.character(marker_input))) %>%
    left_join(., gene_and_functions, by = "gene") %>%
    mutate(gene = factor(gene, 
                         levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene)), 
           mechanism = factor(mechanism, levels = c("DSB", "SSB", "Associated", "nonDDR"))) %>% 
    arrange(mechanism, gene)
  
  left_join(estimates_tmp, carriers_tmp, by = c("gene", "mechanism"))
  
}

supplementary_tables_sensitivity(missense_merged_results, marker_input = "CD8", variants_input = "missense")


