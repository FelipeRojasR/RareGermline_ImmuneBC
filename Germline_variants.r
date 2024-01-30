# Germline variants

# This is a standalone script

## Missense variants
### Manually curated variants
# 
# Filter out all the carriers of BRCA1, BRCA2 and TP53.
# The function that is currently online can be used for the overall analysis and only create a table with the clinvar, enigma and TP53 filtered MSVs.

# ATM (ClinVar)
# c.875C>T (p.Pro292Leu), c.2849T>G (p.Leu950Arg), c.3848T>C (p.Leu1283Pro), c.6200C>A (p.Ala2067Asp), c.6679C>T (p.Arg2227Cys), c.7271T>G (p.Val2424Gly), c.7570G>C (p.Ala2524Pro), c.8122G>A (p.Asp2708Asn), c.8147T>C (p.Val2716Ala), c.8494C>T (p.Arg2832Cys), c.8546G>C (p.Arg2849Pro), c.9022C>T (p.Arg3008Cys), c.9023G>A (p.Arg3008His)
atm_annotated_msv <- tibble(gene = rep("ATM", 13),
                            gene.coord = c("c.875C>T", "c.2849T>G", "c.3848T>C", "c.6200C>A", "c.6679C>T", "c.7271T>G", "c.7570G>C", "c.8122G>A", "c.8147T>C", "c.8494C>T", "c.8546G>C", "c.9022C>T", "c.9023G>A"),
                            protein.coord = c("p.Pro292Leu", "p.Leu950Arg", "p.Leu1283Pro", "p.Ala2067Asp", "p.Arg2227Cys", "p.Val2424Gly", "p.Ala2524Pro", "p.Asp2708Asn", "p.Val2716Ala", "p.Arg2832Cys", "p.Arg2849Pro",
                                              "p.Arg3008Cys", "p.Arg3008His"))

# CHEK2 (ClinVar)
# c.470T>G (p.Ile157Ser), c.433C>T (p.Arg145Trp)
chek2_annotated_msv <- tibble(gene = rep("CHEK2", 2), 
                              gene.coord = c("c.470T>G", "c.433C>T"),
                              protein.coord = c("p.Ile157Ser", "p.Arg145Trp"))

# BRCA1 (ClinVar and ENIGMA)
# c.5339T>C (p.Leu1780Pro), c.5216A>G (p.Asp1739Gly), c.5213G>A (p.Gly1738Glu), c.5207T>C (p.Val1736Ala), c.5095C>T (p.Arg1699Trp), c.5089T>C (p.Cys1697Arg), c.5072C>T (p.Thr1691Ile), c.5057A>G (p.His1686Arg), c.4964C>T (p.Ser1655Phe), c.181T>G (p.Cys61Gly), c.181T>C (p.Cys61Arg), c.130T>A (p.Cys44Ser), c.53T>C (p.Met18Thr)
brca1_annotated_msv <- tibble(gene = rep("BRCA1", 13),
                              gene.coord = c("c.5339T>C", "c.5216A>G", "c.5213G>A", "c.5207T>C", "c.5095C>T", "c.5089T>C", "c.5072C>T", "c.5057A>G", "c.4964C>T", "c.181T>G", "c.181T>C", "c.130T>A", "c.53T>C"),
                              protein.coord = c("p.Leu1780Pro", "p.Asp1739Gly", "p.Gly1738Glu", "p.Val1736Ala", "p.Arg1699Trp", "p.Cys1697Arg", "p.Thr1691Ile", "p.His1686Arg", "p.Ser1655Phe", "p.Cys61Gly", 
                                                "p.Cys61Arg", "p.Cys44Ser", "p.Met18Thr"))

# BRCA2 (ClinVar and ENIGMA)
# c.7529T>C (p.Leu2510Pro), c.7879A>T (p.Ile2627Phe), c.7940T>C (p.Leu2647Pro), c.7958T>C (p.Leu2653Pro), c.8023A>G (p.Ile2675Val), c.8057T>C (p.Leu2686Pro), c.8167G>C (p.Asp2723His), c.8243G>A (p.Gly2748Asp), c.9004G>A (p.Glu3002Lys), c.9154C>T (p.Arg3052Trp), c.9226G>A (p.Gly3076Arg), c.9371A>T (p.Asn3124Ile)
brca2_annotated_msv <- tibble(gene = rep("BRCA2", 12),
                              gene.coord = c("c.7529T>C", "c.7879A>T", "c.7940T>C", "c.7958T>C", "c.8023A>G", "c.8057T>C", "c.8167G>C", "c.8243G>A", "c.9004G>A", "c.9154C>T", "c.9226G>A", "c.9371A>T"),
                              protein.coord = c("p.Leu2510Pro", "p.Ile2627Phe", "p.Leu2647Pro", "p.Leu2653Pro", "p.Ile2675Val", "p.Leu2686Pro", "p.Asp2723His", "p.Gly2748Asp", "p.Glu3002Lys", "p.Arg3052Trp", 
                                                "p.Gly3076Arg", "p.Asn3124Ile"))

# TP53
TP53_UMD_variants <- read_tsv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/UMD_variants_EU.tsv")
colnames(TP53_UMD_variants) <- colnames(TP53_UMD_variants) %>% gsub(" ", "_", .)
# TP53, Ensembl transcript: ENST00000269305.4, NCBI transcript: NM_000546
# Only the information from the `Transcript t1 MN_000546.5 column` is preserved because is the one that 
# matches with the NCBI transcript ID available in the NEJM manuscript.
TP53_UMD_variants.pathogenic <- TP53_UMD_variants %>%
  # From the 6,870 variants included in the p53.fr database, 264 are classified as pathogenic
  filter(Pathogenicity == "Pathogenic") %>%
  # From the 264 variants, 150 variants were found in published research as well as databases analysis 
  # which provides sufficient evidence for classification of these variants as pathogenic.
  filter(Final_comment %in% c("Published research as well as database analysis provide sufficient evidence for classification of this variant as pathogenic.", 
                              "Published research as well as database analysis provide sufficient evidence for classification of this variant as pathogenic. Splicing defects associated with this variant have been observed in multiple tumors.")) %>%
  # The NEJM paper states that the variants described were matched to the build37 (hg19) of the human genome.
  transmute(gene = "TP53",
            gene.coord = cDNA_variant, 
            protein.coord = Protein_P1_TP53_alpha__NP_000537.3)
# These columns may be of interest: HG19_Variant, Transcript_t1_MN_000546.5, Protein_P1_TP53_alpha__NP_000537.3

### NEJM variants
# TP53	c.524G>A
TP53_NEJM_variants.pathogenic <- readxl::read_excel("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/BRIDGES_TP53_missense_variants_NEJM_classifications.xlsx") %>%
  transmute(
    gene = "TP53",
    gene.coord = HGVSc,
    protein.coord = NA
  )


annotated_msv <- rbind(atm_annotated_msv, chek2_annotated_msv, brca1_annotated_msv, brca2_annotated_msv, TP53_NEJM_variants.pathogenic)


# PPM1D: Deleted 
# FAM175A: Renamed to ABRAXAS1
# MRE11A: Renamed to MRE11
# BRE: Renamed to BABAM2
fix_genes_BRIDGES <- function(dataset){
  dataset %>% 
    filter(genes != "PPM1D") %>% 
    mutate(genes = case_when(genes == "FAM175A" ~ "ABRAXAS1", 
                             genes == "MRE11A" ~ "MRE11",
                             genes == "BRE" ~ "BABAM2",
                             TRUE ~ genes)
    )
}


# From the annotation file I need the HGVSc separated by ':' and delete the 'p.'
missense_annotation <- read_csv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/felipe_missense_variants_annotation.csv") %>% 
  select(Chromosome, Position, HGVSc, Ref, Alt, Gene, Impact, CADD_phred, variant, aa_pos) %>%
  separate(HGVSc, c("additional", "gene.coord"), sep = ":", remove = FALSE) %>%
  select(-c(additional))

manually_curated_variants <- missense_annotation %>% filter(Gene == "ATM" & gene.coord %in% filter(annotated_msv, gene == "ATM")$gene.coord | 
                                                              Gene == "CHEK2" & gene.coord %in% filter(annotated_msv, gene == "CHEK2")$gene.coord |
                                                              Gene == "BRCA1" & gene.coord %in% filter(annotated_msv, gene == "BRCA1")$gene.coord |
                                                              Gene == "BRCA2" & gene.coord %in% filter(annotated_msv, gene == "BRCA2")$gene.coord |
                                                              Gene == "TP53" & gene.coord %in% filter(annotated_msv, gene == "TP53")$gene.coord)

###### MISSENSE BRIDGES CARRIERS ######
# CADD score 10 indicates that these are predicted to be the 10% most deleterious substitutions that you can do to the human genome, 
# a score of greater or equal 20 indicates the 1% most deleterious and so on. For CADD, we classified the variants based on the 
# phred-like score with a cutoff 20 [Niroula, et al. 2019. Plos Comput Biol],
missense_annotation_tofilter <- missense_annotation %>%
  mutate(pathogenic = case_when(CADD_phred >= 20 ~ "pathogenic", TRUE ~ "non-pathogenic")) %>%
  filter(pathogenic == "pathogenic")

# This combines the Phred-CADD score and the manual curation
pathogenic_missense_bridges <- read_csv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/schmidt_734_bridges_missense.csv",
                                        show_col_types = FALSE) %>%
  select(-c(...1)) %>% select("BRIDGES_ID", contains("_Variant")) %>% # select columns with _missense.variants label
  pivot_longer(-BRIDGES_ID, names_to = "genes", values_to = "variants") %>%
  mutate(pathogenic.phredCADD = case_when(variants %in% missense_annotation_tofilter$variant ~ 1, TRUE ~ 0),
         pathogenic.manual = case_when(variants %in% manually_curated_variants$variant ~ 1, TRUE ~ 0),
         genes = str_remove(genes, "_Variant")) %>%
  mutate(pathogenic.to.include = case_when(genes %in% unique(manually_curated_variants$Gene) & pathogenic.manual == 1 ~ 1,
                                           genes %!in% unique(manually_curated_variants$Gene) & pathogenic.phredCADD == 1 ~ 1,
                                           TRUE ~ 0)) %>%
  select(BRIDGES_ID, genes, pathogenic.to.include) %>% 
  fix_genes_BRIDGES(.) %>% group_by(BRIDGES_ID) %>% 
  tidyr::spread(-genes, pathogenic.to.include, fill = 0) %>% ungroup()



## Protein truncating variants

truncating_bridges <- read_csv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/schmidt_734_bridges_truncating.csv", 
                               show_col_types = FALSE) %>% select(-c(...1)) %>%
  pivot_longer(-BRIDGES_ID, names_to = "genes", values_to = "truncating_variants") %>%
  mutate(genes = str_remove(genes, "_truncating"), 
         carrier = case_when(truncating_variants >= 1 ~ 1, TRUE ~ truncating_variants)) %>%
  fix_genes_BRIDGES(.) %>% group_by(BRIDGES_ID) %>% select(-truncating_variants) %>% 
  tidyr::spread(-genes, carrier, fill = 0) %>% ungroup() # binary table, single row per Bridges ID

################

# Chr22: 28695869 (on Assembly GRCh38)
# Chr22: 29091857 (on Assembly GRCh37)
# https://gnomad.broadinstitute.org/variant/22-29091856-AG-A?dataset=gnomad_r2_1
# https://www.ensembl.org/Homo_sapiens/Variation/Explore?r=22:28695369-28696369;toggle_HGVS_names=open;v=rs555607708;vdb=variation;vf=194165085

# There are 12 cases with 2 chr22_29091856_AG_A mutations, those were converted to one to avoid inconsistency
truncating.chek2.1100delC <- read_csv("/DATA/users/f.rojas/Projects/immune_infiltration_germline/official_data_BCAC_06122022/schmidt_734_bridges_genotypes.csv") %>%
  mutate(chr22_29091856_AG_A = case_when(chr22_29091856_AG_A > 1 ~ 1, TRUE ~ chr22_29091856_AG_A)) %>%
  select(Bridges_ID, chr22_29091856_AG_A) %>% rename("BRIDGES_ID" = "Bridges_ID",
                                                     "CHEK2_1100delC" = "chr22_29091856_AG_A")
################
truncating_bridges.including.chek2mods <- truncating_bridges %>% 
  inner_join(., truncating.chek2.1100delC, by = "BRIDGES_ID") %>%
  mutate(CHEK2.other = case_when(CHEK2 == 1 & CHEK2_1100delC == 0 ~ 1, TRUE ~ 0))

rm(truncating.chek2.1100delC)
