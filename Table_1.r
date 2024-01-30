# TABLE 1
# Baseline characteristics among IHC markers datasets (Table 1)
dataset_ids_presence <- TilScores_BridgesPhenotype_markers %>% select(BCAC_ID, marker) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = "marker", values_from = "presence", values_fill = 0) %>%
  mutate(all = 1)
  
# Samples with all markers available
all_markers <- filter(.data = dataset_ids_presence, CD8 == 1 & FOXP3 == 1 & CD20 == 1 & CD163 == 1)$BCAC_ID

# All samples included with any of the markers
any_markers <- filter(.data = dataset_ids_presence, all == 1)$BCAC_ID


table(phenotype_BRIDGES_fortables$ER_statusIndex, useNA = "always")

phenotype_BRIDGES_fortables <- read_tsv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/concept_734_schmidt_euro_cases_bridges_pheno_v15_02.txt") %>%
  select(BRIDGES_ID, BCAC_ID, study, AgeDiagIndex, Index_corr, BehaviourIndex, GradeIndex, 
         MorphologygroupIndex_corr, NodeStatusIndex, NodesIndex, SizeIndex_Cat, 
         ER_statusIndex, PR_statusIndex, HER2_statusIndex)

phenotype_BRIDGES_fortables[phenotype_BRIDGES_fortables == 888] <- NA


phenotype_BRIDGES_fortables <- phenotype_BRIDGES_fortables %>%
  mutate(all_markers_ids = case_when(BCAC_ID %in% all_markers ~ 1, TRUE ~ 0),
         any_markers_ids = case_when(BCAC_ID %in% any_markers ~ 1, TRUE ~ 0)) %>%
  filter(any_markers_ids == 1) %>%
  mutate(BehaviourIndex = case_when(BehaviourIndex == 1 ~ "Invasive",
                                    BehaviourIndex == 2 ~ "In-situ",
                                    BehaviourIndex == 888 ~ NA_character_),
         MorphologygroupIndex_corr = case_when(MorphologygroupIndex_corr == "Mixed (ductal & lobular)" ~ "Mixed",
                                               MorphologygroupIndex_corr %in% c("Medullary", "Mucinous", 
                                                                                "Other", "Papillary", "Tubular") ~ "Other",
                                               TRUE ~ MorphologygroupIndex_corr),
         `All cases` = case_when(BCAC_ID %in% filter(dataset_ids_presence, all == 1)$BCAC_ID ~ 1, TRUE ~ 0),
         CD8 = case_when(BCAC_ID %in% filter(dataset_ids_presence, CD8 == 1)$BCAC_ID ~ 1, TRUE ~ 0),
         FOXP3 = case_when(BCAC_ID %in% filter(dataset_ids_presence, FOXP3 == 1)$BCAC_ID ~ 1, TRUE ~ 0),
         CD20 = case_when(BCAC_ID %in% filter(dataset_ids_presence, CD20 == 1)$BCAC_ID ~ 1, TRUE ~ 0),
         CD163 = case_when(BCAC_ID %in% filter(dataset_ids_presence, CD163 == 1)$BCAC_ID ~ 1, TRUE ~ 0)) %>% 
  mutate(BehaviourIndex = factor(BehaviourIndex, levels = c("In-situ", "Invasive"),
                                 labels = c("In-situ", "Invasive")),
         GradeIndex = factor(GradeIndex, levels = c(1, 2, 3),
                             labels = c("well differentiated", "moderately differentiated", "poorly/un-differentiated")),
         MorphologygroupIndex_corr = factor(MorphologygroupIndex_corr, levels = c("Ductal", "Lobular", "Mixed", "Other"), 
                                            labels = c("Ductal", "Lobular", "Mixed", "Other")),
         NodeStatusIndex = factor(NodeStatusIndex, levels = c(0,1), 
                                  labels = c("Negative", "Positive")),
         SizeIndex_Cat = factor(SizeIndex_Cat, levels = c(1,2,3), 
                                labels = c("≤2cm", ">2cm and ≤5cm", ">5cm")),
         ER_statusIndex = factor(ER_statusIndex, levels = c(0,1), 
                                 labels = c("Negative", "Positive")),
         PR_statusIndex = factor(PR_statusIndex, levels = c(0,1), 
                                 labels = c("Negative", "Positive")),
         HER2_statusIndex = factor(HER2_statusIndex, levels = c(0,1), 
                                   labels = c("Negative", "Positive"))
         )

#### ####
table1::label(phenotype_BRIDGES_fortables$AgeDiagIndex) <- "Age, years"
table1::label(phenotype_BRIDGES_fortables$BehaviourIndex) <- "Behaviour"
table1::label(phenotype_BRIDGES_fortables$GradeIndex) <- "Grade"
table1::label(phenotype_BRIDGES_fortables$MorphologygroupIndex_corr) <- "Morphology"
table1::label(phenotype_BRIDGES_fortables$NodeStatusIndex) <- "Node Status"
table1::label(phenotype_BRIDGES_fortables$SizeIndex_Cat) <- "Tumor size"
table1::label(phenotype_BRIDGES_fortables$ER_statusIndex) <- "ER-status"
table1::label(phenotype_BRIDGES_fortables$PR_statusIndex) <- "PR-status"
table1::label(phenotype_BRIDGES_fortables$HER2_statusIndex) <- "HER2-status"

overall_table1 <- phenotype_BRIDGES_fortables %>%
  dplyr::select(-c(all_markers_ids, any_markers_ids)) %>%
  pivot_longer(-c(BRIDGES_ID, BCAC_ID, study, AgeDiagIndex, Index_corr, 
                  BehaviourIndex, GradeIndex, MorphologygroupIndex_corr, 
                  NodeStatusIndex, NodesIndex, SizeIndex_Cat, 
                  ER_statusIndex, PR_statusIndex, HER2_statusIndex),
               names_to = "markers", values_to = "carrier") %>%
  mutate(markers = factor(markers, levels = c("All cases", "CD8", "CD20", "FOXP3", "CD163"))) %>%
  filter(carrier != 0)


table1::table1(~ AgeDiagIndex + ER_statusIndex + PR_statusIndex + HER2_statusIndex + GradeIndex + MorphologygroupIndex_corr + NodeStatusIndex + BehaviourIndex + SizeIndex_Cat | markers, 
               data = overall_table1, overall = FALSE,
               topclass = "Rtable1-zebra")
