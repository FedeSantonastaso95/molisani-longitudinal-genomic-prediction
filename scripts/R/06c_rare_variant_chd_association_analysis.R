################################################################################
# 06c_rare_variant_chd_association_analysis.R
#
# Description:
# This script performs cross-sectional and survival analyses for selected rare
# variants after extraction of genotype dosages from sequence-imputed genotypes
# using PLINK2.
#
# Genotype dosages were exported with PLINK2 using `--export A`. In this dataset,
# the exported dosage corresponds to the major allele. Therefore, genotype
# dosages are recoded so that:
#   0 = non-carrier of the rare allele
#   1 = heterozygous carrier
#   2 = homozygous carrier of the rare allele
#
# This script assumes that the input rare variant table has already been filtered
# upstream according to the predefined quality-control criteria, including
# imputation quality (INFO) and minor allele count (MAC).
#
# The script performs:
#   1. Cross-sectional logistic regression for CHD
#   2. Survival analysis for incident CHD using Cox regression
#
# Survival analysis is restricted to variants showing nominal evidence of
# association in the cross-sectional CHD analysis (P < 0.05).
#
# Inputs:
#   - PLINK2 dosage file (.raw) for selected rare variants
#   - Rare variant table already filtered upstream
#   - SNP list from the sequence-imputed genotype dataset
#   - Sample file from the sequence-imputed genotype dataset
#   - Phenotype/covariate table
#
# Outputs:
#   - Cross-sectional association results
#   - Survival association results
#
# Usage:
#   Rscript 06c_rare_variant_chd_association_analysis.R
#
################################################################################

# =========================
# Load libraries
# =========================
library(data.table)
library(dplyr)
library(survival)

# =========================
# Define input/output paths
# =========================
carriers_file  <- "/path/to/carriers_rare_variants.raw"
variant_file   <- "/path/to/rare_variants_final_filtered.tsv"
snplist_file   <- "/path/to/imputed_genotypes.snplist"
sample_file    <- "/path/to/imputed_genotypes.sample"
phenotype_file <- "/path/to/phenotypes_covariates_prs.txt"
output_dir     <- "/path/to/output_dir"

# Restrict to cardiometabolic traits considered in the analysis
selected_traits <- c(
  "LDL", "HDL", "TotChol", "TRI", "LPA", "ApoA", "ApoB",
  "glucose", "insulin", "BMI", "SBP", "Uric", "HeartRate", "ALT"
)

# =========================
# Helper functions
# =========================
clean_plink_variant_names <- function(x) {
  sub("_(?!.*_).*$", "", x, perl = TRUE)
}

recode_major_to_rare_dosage <- function(x) {
  # PLINK2 dosage corresponds to the major allele:
  # major dosage 2 -> rare dosage 0
  # major dosage 1 -> rare dosage 1
  # major dosage 0 -> rare dosage 2
  x <- round(x)
  out <- x
  out[x == 0] <- 2
  out[x == 2] <- 0
  out
}

collapse_homozygotes_if_rare <- function(x, min_hom = 4) {
  n_hom <- sum(x == 2, na.rm = TRUE)
  if (n_hom > 0 && n_hom < min_hom) {
    x[x == 2] <- 1
  }
  x
}

compute_mac <- function(x) {
  n_het <- sum(x == 1, na.rm = TRUE)
  n_hom <- sum(x == 2, na.rm = TRUE)
  n_het + 2 * n_hom
}

prepare_carrier_matrix <- function(carriers_file) {
  dt <- fread(carriers_file)
  
  if (ncol(dt) < 7) {
    stop("Unexpected PLINK .raw format: fewer than 7 columns found.")
  }
  
  dt[, 7:ncol(dt) := lapply(.SD, recode_major_to_rare_dosage), .SDcols = 7:ncol(dt)]
  
  old_names <- names(dt)[7:ncol(dt)]
  new_names <- clean_plink_variant_names(old_names)
  setnames(dt, old = old_names, new = new_names)
  
  dt
}

prepare_variant_map <- function(variant_file, snplist_file) {
  var <- fread(variant_file)
  snplist <- fread(snplist_file)
  
  var %>%
    left_join(
      snplist,
      by = c(
        "GENPOS"  = "V4",
        "CHROM"   = "V1",
        "ALLELE0" = "V5",
        "ALLELE1" = "V6"
      )
    ) %>%
    rename(bgen_ID = V2) %>%
    dplyr::select(-V3)
}

prepare_sample_table <- function(sample_file) {
  fread(sample_file) %>%
    .[-1, ]
}

prepare_phenotype_table <- function(phenotype_file, carriers_dt) {
  fread(phenotype_file, data.table = FALSE) %>%
    mutate(
      Dataexit_Bc_Fup2020  = as.character(Dataexit_Bc_Fup2020),
      Dataexit_Pc_Fup2020  = as.character(Dataexit_Pc_Fup2020),
      Dataexit_Crc_Fup2020 = as.character(Dataexit_Crc_Fup2020),
      Dataexit_Lc_Fup2020  = as.character(Dataexit_Lc_Fup2020),
      Dataexit_Rnc_Fup2020 = as.character(Dataexit_Rnc_Fup2020),
      Dataexit_Pnc_Fup2020 = as.character(Dataexit_Pnc_Fup2020)
    ) %>%
    mutate(
      Date_cancer_FUP2020 = case_when(
        Dataexit_Bc_Fup2020  != "2020-12-31" ~ Dataexit_Bc_Fup2020,
        Dataexit_Pc_Fup2020  != "2020-12-31" ~ Dataexit_Pc_Fup2020,
        Dataexit_Crc_Fup2020 != "2020-12-31" ~ Dataexit_Crc_Fup2020,
        Dataexit_Lc_Fup2020  != "2020-12-31" ~ Dataexit_Lc_Fup2020,
        Dataexit_Rnc_Fup2020 != "2020-12-31" ~ Dataexit_Rnc_Fup2020,
        Dataexit_Pnc_Fup2020 != "2020-12-31" ~ Dataexit_Pnc_Fup2020,
        TRUE ~ "2020-12-31"
      ),
      Fnf_cancer_FUP2020 = ifelse(
        Fnf_Bc_Fup2020 == 1 | Fnf_Pc_Fup2020 == 1 | Fnf_Crc_Fup2020 == 1 |
          Fnf_Lc_Fup2020 == 1 | Fnf_Rnc_Fup2020 == 1 | Fnf_Pnc_Fup2020 == 1,
        1, 0
      ),
      Fnf_cancer_FUP2020 = ifelse(is.na(Fnf_cancer_FUP2020), 0, Fnf_cancer_FUP2020)
    ) %>%
    left_join(carriers_dt, by = c("FID" = "FID", "IID" = "IID")) %>%
    mutate(
      Chd_Bs_New_final = ifelse(Chd_Bs_New == 0, 0, 1),
      Chd0_Fnof_Fup2020 = if_else(
        Chd_Bs_New_final == 1 & Chd0_Fnof_Fup2020 == 1,
        as.numeric(NA),
        as.numeric(Chd0_Fnof_Fup2020)
      )
    )
}

run_cross_sectional_logistic <- function(df, variant_map) {
  unique_combos <- unique(as.data.table(variant_map)[, .(gene_name, trait)])
  all_results <- list()
  
  df <- df %>%
    mutate(
      CHD_event_incident_prevalent = if_else(
        Chd0_Fnof_Fup2020 == 1 | Chd_Bs_New == 1,
        1, 0
      )
    )
  
  for (i in seq_len(nrow(unique_combos))) {
    this_gene  <- unique_combos$gene_name[i]
    this_trait <- unique_combos$trait[i]
    
    subset_variants <- variant_map %>%
      filter(gene_name == this_gene, trait == this_trait) %>%
      pull(bgen_ID) %>%
      unique()
    
    subset_variants <- subset_variants[!is.na(subset_variants) & subset_variants != ""]
    if (length(subset_variants) == 0) next
    
    for (variant in subset_variants) {
      if (!(variant %in% colnames(df))) next
      
      df_variant <- df
      df_variant[[variant]] <- collapse_homozygotes_if_rare(df_variant[[variant]], min_hom = 4)
      
      MAC <- compute_mac(df_variant[[variant]])
      n_homoz <- sum(df_variant[[variant]] == 2, na.rm = TRUE)
      
      if (MAC < 1) next
      
      formula <- as.formula(
        paste0(
          "CHD_event_incident_prevalent ~ factor(", variant,
          ") + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5"
        )
      )
      
      tryCatch({
        model <- glm(formula, data = df_variant, family = binomial)
        s <- summary(model)
        
        variant_row <- grep(paste0("^factor\\(", variant, "\\)"), rownames(s$coefficients))
        if (length(variant_row) == 0) next
        
        CI <- tryCatch(confint(model), error = function(e) NULL)
        OR <- exp(coef(model))
        CI_exp <- if (!is.null(CI)) exp(CI) else NULL
        p_vals <- coef(summary(model))[, "Pr(>|z|)"]
        
        result <- data.frame(
          gene     = this_gene,
          trait    = this_trait,
          variant  = variant,
          OR       = round(OR[variant_row], 2),
          CI_lower = if (!is.null(CI_exp)) round(CI_exp[variant_row, 1], 2) else NA,
          CI_upper = if (!is.null(CI_exp)) round(CI_exp[variant_row, 2], 2) else NA,
          P_value  = round(p_vals[variant_row], 4),
          MAC_analysis = MAC,
          N_homoz  = n_homoz
        )
        
        all_results[[paste(this_gene, this_trait, variant, sep = "_")]] <- result
      }, error = function(e) {
        message("Logistic model failed for ", variant, ": ", e$message)
      })
    }
  }
  
  rbindlist(all_results, fill = TRUE)
}

run_survival_cox <- function(df, variant_map) {
  unique_pairs <- unique(as.data.table(variant_map)[, .(gene_name, trait)])
  all_results <- list()
  
  for (i in seq_len(nrow(unique_pairs))) {
    this_gene  <- unique_pairs$gene_name[i]
    this_trait <- unique_pairs$trait[i]
    
    subset_variants <- variant_map %>%
      filter(gene_name == this_gene, trait == this_trait) %>%
      pull(bgen_ID) %>%
      unique()
    
    subset_variants <- subset_variants[!is.na(subset_variants) & subset_variants != ""]
    if (length(subset_variants) == 0) next
    
    for (variant in subset_variants) {
      if (!(variant %in% colnames(df))) next
      if (all(is.na(df[[variant]]))) next
      
      df_variant <- df
      df_variant[[variant]] <- collapse_homozygotes_if_rare(df_variant[[variant]], min_hom = 4)
      
      MAC <- compute_mac(df_variant[[variant]])
      n_homoz <- sum(df_variant[[variant]] == 2, na.rm = TRUE)
      
      if (MAC < 1) next
      
      formula <- as.formula(
        paste0(
          "Surv(time, status) ~ factor(", variant,
          ") + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5"
        )
      )
      
      tryCatch({
        model <- coxph(formula, data = df_variant)
        s <- summary(model)
        
        variant_row <- grep(paste0("^factor\\(", variant, "\\)"), rownames(s$coefficients))
        if (length(variant_row) == 0) next
        
        result <- data.frame(
          gene     = this_gene,
          trait    = this_trait,
          variant  = variant,
          HR       = round(s$coefficients[variant_row, "exp(coef)"], 2),
          CI_lower = round(s$conf.int[variant_row, "lower .95"], 2),
          CI_upper = round(s$conf.int[variant_row, "upper .95"], 2),
          P_value  = signif(s$coefficients[variant_row, "Pr(>|z|)"], 3),
          MAC_analysis = MAC,
          N_homoz  = n_homoz
        )
        
        all_results[[paste(this_gene, this_trait, variant, sep = "_")]] <- result
      }, error = function(e) {
        message("Cox model failed for ", variant, ": ", e$message)
      })
    }
  }
  
  rbindlist(all_results, fill = TRUE)
}

# =========================
# Load and prepare data
# =========================
df_carr <- prepare_carrier_matrix(carriers_file)
variant_map <- prepare_variant_map(variant_file, snplist_file)

sample_dt <- prepare_sample_table(sample_file)
pheno_dt  <- prepare_phenotype_table(phenotype_file, df_carr)

df_all <- sample_dt %>%
  left_join(pheno_dt, by = c("ID_1" = "FID"))

# =========================
# Cross-sectional analysis
# =========================
cross_results <- run_cross_sectional_logistic(df_all, variant_map)

cross_final <- variant_map %>%
  left_join(cross_results, by = c("gene_name" = "gene", "trait", "bgen_ID" = "variant")) %>%
  filter(trait %in% selected_traits) %>%
  filter(!is.na(OR)) %>%
  distinct()

fwrite(
  cross_final,
  file = file.path(output_dir, "rare_variants_cross_sectional_chd.txt"),
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# =========================
# Select variants for survival analysis
# =========================
# Survival analysis is restricted to variants showing nominal evidence
# of association in the cross-sectional CHD analysis (P < 0.05).
survival_variants <- cross_final %>%
  filter(!is.na(P_value), P_value <= 0.05) %>%
  dplyr::select(gene_name, trait, bgen_ID) %>%
  distinct()

# =========================
# Prepare survival dataset
# =========================
df_surv <- df_all %>%
  mutate(
    Recruitment_date = as.Date(gsub("/", "-", Recruitment_date), format = "%d-%m-%Y"),
    time = Annipers_Chd0_Fnof_Fup2020,
    status = Chd0_Fnof_Fup2020
  ) %>%
  filter(Chd_Bs_New_final != 1)

# =========================
# Survival analysis
# =========================
surv_results <- run_survival_cox(df_surv, survival_variants)

surv_final <- survival_variants %>%
  left_join(surv_results, by = c("gene_name" = "gene", "trait", "bgen_ID" = "variant"))

fwrite(
  surv_final,
  file = file.path(output_dir, "rare_variants_survival_chd.txt"),
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)