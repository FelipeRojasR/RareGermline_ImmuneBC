# TABLE 1
# Baseline characteristics among IHC markers datasets (Table 1)
```{r COUNT_PATIENTS}
IDS_by_marker <- rbind(CD8_tils_BRIDGES %>% select(BCAC_ID), 
      FOXP3_tils_BRIDGES %>% select(BCAC_ID), 
      CD20_tils_BRIDGES %>% select(BCAC_ID),
      CD163_tils_BRIDGES %>% select(BCAC_ID)) %>%
  distinct(BCAC_ID) %>%
  mutate(CD8 = case_when(BCAC_ID %in% CD8_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0),
         FOXP3 = case_when(BCAC_ID %in% FOXP3_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0),
         CD20 = case_when(BCAC_ID %in% CD20_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0),
         CD163 = case_when(BCAC_ID %in% CD163_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0)
  )

all_markers <- IDS_by_marker %>% filter(CD8 == 1 & FOXP3 == 1 & CD20 == 1 & CD163 == 1) %>% .$BCAC_ID
any_markers <- IDS_by_marker %>% filter(CD8 == 1 | FOXP3 == 1 | CD20 == 1 | CD163 == 1) %>% .$BCAC_ID
length(all_markers)
length(any_markers)


phenotype_BRIDGES_fortables <- read_tsv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/concept_734_schmidt_euro_cases_bridges_pheno_v15_02.txt") %>%
  select(BRIDGES_ID, BCAC_ID, study, AgeDiagIndex, Index_corr, BehaviourIndex, GradeIndex, 
         MorphologygroupIndex_corr, NodeStatusIndex, NodesIndex, SizeIndex_Cat, 
         ER_statusIndex, PR_statusIndex, HER2_statusIndex) %>%
  # filter(BCAC_ID %in% any_markers) %>%
  mutate(all_markers_ids = case_when(BCAC_ID %in% all_markers ~ 1, TRUE ~ 0)) %>%
  mutate(AgeDiagIndex = case_when(AgeDiagIndex == 888 ~ NA_real_, TRUE ~ AgeDiagIndex),
         BehaviourIndex = case_when(BehaviourIndex == 1 ~ "Invasive",
                                    BehaviourIndex == 2 ~ "In-situ",
                                    BehaviourIndex == 888 ~ NA_character_),
         GradeIndex = case_when(GradeIndex == 888 ~ NA_real_, TRUE ~ GradeIndex),
         MorphologygroupIndex_corr = case_when(MorphologygroupIndex_corr == "Mixed (ductal & lobular)" ~ "Mixed",
                                               MorphologygroupIndex_corr %in% c("Medullary", "Mucinous", 
                                                                                "Other", "Papillary", "Tubular") ~ "Other",
                                               MorphologygroupIndex_corr == 888 ~ NA_character_,
                                               TRUE ~ MorphologygroupIndex_corr),
         NodeStatusIndex = case_when(NodeStatusIndex == 888 ~ NA_real_, TRUE ~ NodeStatusIndex),
         SizeIndex_Cat = case_when(SizeIndex_Cat == 888 ~ NA_real_, TRUE ~ SizeIndex_Cat),
         ER_statusIndex = case_when(ER_statusIndex == 888 ~ NA_real_, TRUE ~ ER_statusIndex),
         PR_statusIndex = case_when(PR_statusIndex == 888 ~ NA_real_, TRUE ~ PR_statusIndex),
         HER2_statusIndex = case_when(HER2_statusIndex == 888 ~ NA_real_, TRUE ~ HER2_statusIndex),
         CD8 = case_when(BCAC_ID %in% CD8_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0),
         FOXP3 = case_when(BCAC_ID %in% FOXP3_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0),
         CD20 = case_when(BCAC_ID %in% CD20_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0),
         CD163 = case_when(BCAC_ID %in% CD163_tils_BRIDGES$BCAC_ID ~ 1, TRUE ~ 0)) %>% 
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
                                   labels = c("Negative", "Positive")),
         all_markers_ids = factor(all_markers_ids, levels = c(0,1),
                                  labels = c("Not all Markers", "All Markers"))
         )

#### ####
table1::label(phenotype_BRIDGES_fortables$AgeDiagIndex) <- "Age"
table1::label(phenotype_BRIDGES_fortables$BehaviourIndex) <- "Behaviour"
table1::label(phenotype_BRIDGES_fortables$GradeIndex) <- "Grade"
table1::label(phenotype_BRIDGES_fortables$MorphologygroupIndex_corr) <- "Morphology"
table1::label(phenotype_BRIDGES_fortables$NodeStatusIndex) <- "Node Status"
table1::label(phenotype_BRIDGES_fortables$SizeIndex_Cat) <- "Size"
table1::label(phenotype_BRIDGES_fortables$ER_statusIndex) <- "ER"
table1::label(phenotype_BRIDGES_fortables$PR_statusIndex) <- "PR"
table1::label(phenotype_BRIDGES_fortables$HER2_statusIndex) <- "HER2"
table1::label(phenotype_BRIDGES_fortables$all_markers_ids) <- "Markers"


pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For numeric variables, perform a standard 2-sample t-test
        p <- t.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
    }


overall_table1 <- phenotype_BRIDGES_fortables %>%
  pivot_longer(-c(BRIDGES_ID, BCAC_ID, study, AgeDiagIndex, Index_corr, 
                  BehaviourIndex, GradeIndex, MorphologygroupIndex_corr, 
                  NodeStatusIndex, NodesIndex, SizeIndex_Cat, 
                  ER_statusIndex, PR_statusIndex, HER2_statusIndex, all_markers_ids), 
               names_to = "markers", values_to = "carrier") %>%
  filter(carrier == 1) %>%
  mutate(markers = factor(markers, levels = c("CD8", "CD20", "FOXP3", "CD163")))

table1::table1(~ AgeDiagIndex + ER_statusIndex + HER2_statusIndex + GradeIndex + MorphologygroupIndex_corr + NodeStatusIndex + BehaviourIndex + SizeIndex_Cat | markers, 
               data = overall_table1, overall = FALSE,
               topclass = "Rtable1-zebra")

```