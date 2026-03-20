############################################################
# Supplementary analysis
# Cox interaction models for CHD PRS and rare variant carrier status
#
# Description:
# This script tests whether the association between CHD polygenic
# risk score (PRS) and incident coronary heart disease (CHD) differs
# according to carrier status of selected rare variants.
#
# Rare variants tested:
#   - rs3798220   (LPA / Lp(a))
#   - rs58757394  (SLC4A11 / SBP)
#   - rs730882080 (LDLR / LDL)
#
# Important genotype note:
#   - Genotype dosages were exported with PLINK2 using `--export A`
#   - In this dataset, exported dosages correspond to the major allele
#   - Dosages are therefore recoded so that:
#       0 = non-carrier of the minor allele
#       1 = heterozygous carrier
#       2 = homozygous carrier of the minor allele
#
# Main analyses:
#   1. Linear model:
#        PRS ~ rare variant carrier status + covariates
#   2. Logistic model:
#        rare variant carrier status ~ PRS + covariates
#   3. Cox models:
#        - rare variant only + covariates
#        - PRS only + covariates
#        - rare variant + PRS + covariates
#        - rare variant * PRS + covariates
#   4. Likelihood ratio test comparing joint vs interaction Cox model
#
# Covariate adjustment:
#   - Age
#   - Sex
#   - PC1, PC2, PC3, PC4, PC5
#
# Study sample:
#   - sequence-imputed samples only
#   - participants with non-missing CHD follow-up and covariates
#
# Outputs:
#   - Supplementary table for PRS × rare variant interaction models
#   - Comparison table for effect estimates before/after mutual adjustment
#   - Ancillary tables for PRS ~ carrier and carrier ~ PRS analyses
############################################################


## ----------------------------
## 0. Clean workspace
## ----------------------------
rm(list = ls())


## ----------------------------
## 1. Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(survival)


## ----------------------------
## 2. Input / output paths
## ----------------------------
carrier_file <- "/group/soranzo/f.santonastaso/prj04_quantitative_prs/burdenTest/carriers_rarVar.raw"

variant_file <- "/group/soranzo/f.santonastaso/prj04_quantitative_prs/burdenTest/2025.06.30_listVariants_Burden_ALL_final.txt"

snplist_file <- "/processing_data/shared_datasets/molisani/imputed_genotypes_TopMed/new_6ksequences_molisani_imputed_from_Molisanisequences/molisani_afterQC_imputed_FromMolisaniSequences_6811sequences_norefpanelsamples.snplist"

sample_file <- "/processing_data/shared_datasets/molisani/imputed_genotypes_TopMed/new_6ksequences_molisani_imputed_from_Molisanisequences/molisani_afterQC_imputed_FromMolisaniSequences_6811sequences_norefpanelsamples.sample"

pheno_file <- "/processing_data/shared_datasets/molisani/afterQC_data/phenotypes/db_ht_fup20250210_Phenotypes_Covariates_PRSeur_PRSmoli_ldpred2_mortality_hospitalizations_thirdDBUpdated.txt"

output_dir <- "/group/soranzo/f.santonastaso/prj03_molisani_longitudinalprs/survival_files"

output_interaction_table <- file.path(output_dir, "2026.03.20_Supplementary_Table_RV_PRS_interaction.txt")
output_comparison_table  <- file.path(output_dir, "2026.03.20_Cox_joint_model_comparison_RV_PRS.txt")
output_lm_table          <- file.path(output_dir, "2026.03.20_PRS_vs_carrier_linear.txt")
output_glm_table         <- file.path(output_dir, "2026.03.20_carrier_vs_PRS_logistic.txt")


## ----------------------------
## 3. User-defined settings
## ----------------------------
covars <- c("Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5")

variant_info <- tibble(
  rv_var = c("rv_3798", "rv_5875", "rv_7308"),
  genotype_var = c("rs3798220", "rs58757394", "rs730882080"),
  Variant = c("rs3798220", "rs58757394", "rs730882080"),
  Gene = c("LPA", "SLC4A11", "LDLR"),
  Trait = c("Lp(a)", "SBP", "LDL")
)


## ----------------------------
## 4. Helper: format p-values
## ----------------------------
format_pvalue <- function(p) {
  if (is.na(p)) return(NA_character_)
  out <- formatC(p, digits = 2, format = "e")
  out <- gsub("e\\+?(-?\\d+)", "×10^\\1", out)
  out
}


## ----------------------------
## 5. Helper: extract model terms
## ----------------------------
extract_term_cox <- function(fit, term) {
  s <- summary(fit)
  
  tibble(
    term = term,
    beta = s$coef[term, "coef"],
    HR = s$coef[term, "exp(coef)"],
    CI_low = s$conf.int[term, "lower .95"],
    CI_high = s$conf.int[term, "upper .95"],
    p = s$coef[term, "Pr(>|z|)"]
  )
}

extract_term_lm <- function(fit, term) {
  s <- summary(fit)
  ci <- confint(fit)
  
  tibble(
    term = term,
    beta = coef(fit)[term],
    CI_low = ci[term, 1],
    CI_high = ci[term, 2],
    p = coef(s)[term, "Pr(>|t|)"]
  )
}

extract_term_glm <- function(fit, term) {
  s <- summary(fit)
  ci <- suppressMessages(confint.default(fit))
  
  tibble(
    term = term,
    beta = coef(fit)[term],
    OR = exp(coef(fit)[term]),
    CI_low = exp(ci[term, 1]),
    CI_high = exp(ci[term, 2]),
    p = coef(s)[term, "Pr(>|z|)"]
  )
}


## ----------------------------
## 6. Helper: formatted summaries
## ----------------------------
fmt_hr_ci <- function(hr, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
}

fmt_or_ci <- function(or, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
}

fmt_beta_ci <- function(beta, lo, hi) {
  sprintf("%.3f (%.3f to %.3f)", beta, lo, hi)
}


## ----------------------------
## 7. Load carrier file
## ----------------------------
df_carr <- fread(carrier_file)

# Round dosage columns
df_carr[, 7:ncol(df_carr)] <- lapply(df_carr[, 7:ncol(df_carr)], round)

# Recode dosages:
# original export refers to major allele dosage
# final coding:
#   0 = non-carrier of minor allele
#   1 = heterozygous carrier
#   2 = homozygous carrier of minor allele
df_carr[, 7:ncol(df_carr)] <- lapply(df_carr[, 7:ncol(df_carr)], function(col) {
  if (is.numeric(col)) {
    col_new <- col
    col_new[col == 0] <- 2
    col_new[col == 2] <- 0
    return(col_new)
  } else {
    return(col)
  }
})

# Clean variant names by removing suffix after last underscore
old_names <- names(df_carr)[7:ncol(df_carr)]
new_names <- sub("_(?!.*_).*$", "", old_names, perl = TRUE)
names(df_carr)[7:ncol(df_carr)] <- new_names


## ----------------------------
## 8. Load variant annotation files
## ----------------------------
var_df <- fread(variant_file)
snplist_df <- fread(snplist_file)

var_annot <- var_df %>%
  left_join(
    snplist_df,
    by = c(
      "GENPOS" = "V4",
      "CHROM" = "V1",
      "ALLELE0" = "V5",
      "ALLELE1" = "V6"
    )
  ) %>%
  rename(bgen_ID = V2) %>%
  select(-V3)

setDT(var_annot)
unique_gene_trait <- unique(var_annot[, .(gene_name, trait)])


## ----------------------------
## 9. Load sample file
## ----------------------------
sample_df <- fread(sample_file)

# Remove first PLINK/BGEN dummy row
sample_df <- sample_df[-1, ]


## ----------------------------
## 10. Load phenotype data
## ----------------------------
df_pheno <- fread(
  pheno_file,
  data.table = FALSE
) %>%
  left_join(df_carr, by = c("FID", "IID"))


## ----------------------------
## 11. Restrict to sequence-imputed samples
## ----------------------------
df_seq <- sample_df %>%
  left_join(df_pheno, by = c("ID_1" = "FID"))


## ----------------------------
## 12. Keep analysis variables
## ----------------------------
df_analysis <- df_seq %>%
  select(
    IID,
    FID = ID_1,
    Chd0_Fnof_Fup2020,
    Annipers_Chd0_Fnof_Fup2020,
    rs3798220,
    rs58757394,
    rs730882080,
    PRED_eurW_CHD,
    Bio_Lp_A,
    SBP,
    LDL,
    Age,
    Sex,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5
  ) %>%
  mutate(
    Sex = as.factor(Sex),
    prs_std = as.numeric(scale(PRED_eurW_CHD)),
    rv_3798 = if_else(rs3798220 == 1, 1, 0, missing = NA_real_),
    rv_5875 = if_else(rs58757394 == 1, 1, 0, missing = NA_real_),
    rv_7308 = if_else(rs730882080 == 1, 1, 0, missing = NA_real_)
  ) %>%
  filter(
    !is.na(prs_std),
    !is.na(Annipers_Chd0_Fnof_Fup2020),
    !is.na(Chd0_Fnof_Fup2020)
  ) %>%
  filter(if_all(all_of(covars), ~ !is.na(.)))


cat("N included in rare variant analyses:", nrow(df_analysis), "\n")
cat("N CHD events:", sum(df_analysis$Chd0_Fnof_Fup2020 == 1, na.rm = TRUE), "\n")


## ----------------------------
## 13. Main function for variant analysis
## ----------------------------
run_variant_analysis <- function(data, rv_var, variant_label, gene_label, trait_label, covars) {
  
  dat <- data %>%
    filter(!is.na(.data[[rv_var]]))
  
  covar_str <- paste(covars, collapse = " + ")
  surv_str <- "Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020)"
  
  ## 13A. PRS ~ carrier status
  f_lm <- as.formula(
    paste0("prs_std ~ ", rv_var, " + ", covar_str)
  )
  fit_lm <- lm(f_lm, data = dat)
  
  lm_res <- extract_term_lm(fit_lm, rv_var) %>%
    mutate(
      Variant = variant_label,
      Gene = gene_label,
      Trait = trait_label,
      Model = "LM_PRS_on_carrier",
      N = nobs(fit_lm)
    )
  
  ## 13B. carrier status ~ PRS
  f_glm <- as.formula(
    paste0(rv_var, " ~ prs_std + ", covar_str)
  )
  fit_glm <- glm(f_glm, data = dat, family = binomial())
  
  glm_res <- extract_term_glm(fit_glm, "prs_std") %>%
    mutate(
      Variant = variant_label,
      Gene = gene_label,
      Trait = trait_label,
      Model = "GLM_carrier_on_PRS",
      N = nobs(fit_glm)
    )
  
  ## 13C. Cox models
  f_m0 <- as.formula(
    paste0(surv_str, " ~ ", rv_var, " + ", covar_str)
  )
  fit_m0 <- coxph(f_m0, data = dat)
  
  f_m1 <- as.formula(
    paste0(surv_str, " ~ prs_std + ", covar_str)
  )
  fit_m1 <- coxph(f_m1, data = dat)
  
  f_m2 <- as.formula(
    paste0(surv_str, " ~ ", rv_var, " + prs_std + ", covar_str)
  )
  fit_m2 <- coxph(f_m2, data = dat)
  
  f_m3 <- as.formula(
    paste0(surv_str, " ~ ", rv_var, " * prs_std + ", covar_str)
  )
  fit_m3 <- coxph(f_m3, data = dat)
  
  ## Extract Cox terms
  rv_m0 <- extract_term_cox(fit_m0, rv_var) %>%
    mutate(Model = "Cox_RV_only")
  
  prs_m1 <- extract_term_cox(fit_m1, "prs_std") %>%
    mutate(Model = "Cox_PRS_only")
  
  rv_m2 <- extract_term_cox(fit_m2, rv_var) %>%
    mutate(Model = "Cox_joint")
  
  prs_m2 <- extract_term_cox(fit_m2, "prs_std") %>%
    mutate(Model = "Cox_joint")
  
  int_term <- paste0(rv_var, ":prs_std")
  if (!(int_term %in% rownames(summary(fit_m3)$coef))) {
    int_term <- paste0("prs_std:", rv_var)
  }
  
  rv_m3 <- extract_term_cox(fit_m3, rv_var) %>%
    mutate(Model = "Cox_interaction")
  
  prs_m3 <- extract_term_cox(fit_m3, "prs_std") %>%
    mutate(Model = "Cox_interaction")
  
  int_m3 <- extract_term_cox(fit_m3, int_term) %>%
    mutate(Model = "Cox_interaction")
  
  cox_terms <- bind_rows(
    rv_m0,
    prs_m1,
    rv_m2,
    prs_m2,
    rv_m3,
    prs_m3,
    int_m3
  ) %>%
    mutate(
      Variant = variant_label,
      Gene = gene_label,
      Trait = trait_label,
      N = fit_m2$n,
      Events = fit_m2$nevent,
      `HR (95% CI)` = fmt_hr_ci(HR, CI_low, CI_high),
      `P-value` = vapply(p, format_pvalue, character(1))
    )
  
  ## Likelihood ratio test for interaction
  lrt_int_p <- anova(fit_m2, fit_m3, test = "LRT")[2, "Pr(>|Chi|)"]
  
  ## Compare effect estimates before and after mutual adjustment
  rv_unadj <- extract_term_cox(fit_m0, rv_var)
  rv_adj   <- extract_term_cox(fit_m2, rv_var)
  
  prs_unadj <- extract_term_cox(fit_m1, "prs_std")
  prs_adj   <- extract_term_cox(fit_m2, "prs_std")
  
  rv_compare <- tibble(
    Variant = variant_label,
    Gene = gene_label,
    Trait = trait_label,
    Comparison = "RV_effect_before_vs_after_PRS_adjustment",
    HR_unadjusted = rv_unadj$HR,
    CI_low_unadjusted = rv_unadj$CI_low,
    CI_high_unadjusted = rv_unadj$CI_high,
    p_unadjusted = rv_unadj$p,
    HR_adjusted = rv_adj$HR,
    CI_low_adjusted = rv_adj$CI_low,
    CI_high_adjusted = rv_adj$CI_high,
    p_adjusted = rv_adj$p,
    pct_change_beta = 100 * (rv_adj$beta - rv_unadj$beta) / abs(rv_unadj$beta),
    LRT_interaction_p = lrt_int_p,
    N = fit_m2$n,
    Events = fit_m2$nevent
  )
  
  prs_compare <- tibble(
    Variant = variant_label,
    Gene = gene_label,
    Trait = trait_label,
    Comparison = "PRS_effect_before_vs_after_RV_adjustment",
    HR_unadjusted = prs_unadj$HR,
    CI_low_unadjusted = prs_unadj$CI_low,
    CI_high_unadjusted = prs_unadj$CI_high,
    p_unadjusted = prs_unadj$p,
    HR_adjusted = prs_adj$HR,
    CI_low_adjusted = prs_adj$CI_low,
    CI_high_adjusted = prs_adj$CI_high,
    p_adjusted = prs_adj$p,
    pct_change_beta = 100 * (prs_adj$beta - prs_unadj$beta) / abs(prs_unadj$beta),
    LRT_interaction_p = lrt_int_p,
    N = fit_m2$n,
    Events = fit_m2$nevent
  )
  
  comparison_tab <- bind_rows(rv_compare, prs_compare)
  
  ## Console summary
  cat("\n==================================================\n")
  cat("Variant:", variant_label, "| Gene:", gene_label, "| Trait:", trait_label, "\n")
  cat("N =", fit_m2$n, "| Events =", fit_m2$nevent, "\n")
  
  cat("\nJoint Cox model:\n")
  print(
    cox_terms %>%
      filter(Model == "Cox_joint", term %in% c(rv_var, "prs_std")) %>%
      select(term, `HR (95% CI)`, `P-value`)
  )
  
  cat("\nInteraction LRT p-value:\n")
  print(lrt_int_p)
  
  list(
    lm_res = lm_res,
    glm_res = glm_res,
    cox_terms = cox_terms,
    comparison_tab = comparison_tab,
    fit_m0 = fit_m0,
    fit_m1 = fit_m1,
    fit_m2 = fit_m2,
    fit_m3 = fit_m3,
    lrt_int_p = lrt_int_p
  )
}


## ----------------------------
## 14. Run analyses for each variant
## ----------------------------
results_list <- pmap(
  list(
    rv_var = variant_info$rv_var,
    variant_label = variant_info$Variant,
    gene_label = variant_info$Gene,
    trait_label = variant_info$Trait
  ),
  ~ run_variant_analysis(
    data = df_analysis,
    rv_var = ..1,
    variant_label = ..2,
    gene_label = ..3,
    trait_label = ..4,
    covars = covars
  )
)

names(results_list) <- variant_info$Variant


## ----------------------------
## 15. Build final tables
## ----------------------------

### 15A. Linear model table
tab_prs_carrier <- bind_rows(lapply(results_list, `[[`, "lm_res")) %>%
  mutate(
    `Beta (95% CI)` = fmt_beta_ci(beta, CI_low, CI_high),
    `P-value` = vapply(p, format_pvalue, character(1))
  ) %>%
  select(
    Variant, Gene, Trait, Model, N,
    `Beta (95% CI)`, `P-value`,
    beta, CI_low, CI_high, p
  )

### 15B. Logistic model table
tab_carrier_prs_logistic <- bind_rows(lapply(results_list, `[[`, "glm_res")) %>%
  mutate(
    `OR (95% CI)` = fmt_or_ci(OR, CI_low, CI_high),
    `P-value` = vapply(p, format_pvalue, character(1))
  ) %>%
  select(
    Variant, Gene, Trait, Model, N,
    `OR (95% CI)`, `P-value`,
    OR, CI_low, CI_high, p
  )

### 15C. Cox comparison table
tab_cox_comparison <- bind_rows(lapply(results_list, `[[`, "comparison_tab")) %>%
  mutate(
    `HR_unadjusted (95% CI)` = fmt_hr_ci(HR_unadjusted, CI_low_unadjusted, CI_high_unadjusted),
    `HR_adjusted (95% CI)` = fmt_hr_ci(HR_adjusted, CI_low_adjusted, CI_high_adjusted),
    `P_unadjusted` = vapply(p_unadjusted, format_pvalue, character(1)),
    `P_adjusted` = vapply(p_adjusted, format_pvalue, character(1)),
    `LRT interaction P-value` = vapply(LRT_interaction_p, format_pvalue, character(1))
  ) %>%
  select(
    Variant, Gene, Trait, Comparison, N, Events,
    `HR_unadjusted (95% CI)`, `P_unadjusted`,
    `HR_adjusted (95% CI)`, `P_adjusted`,
    pct_change_beta,
    `LRT interaction P-value`,
    everything()
  )

### 15D. Supplementary interaction table
tab_interaction <- bind_rows(lapply(results_list, `[[`, "cox_terms")) %>%
  filter(Model == "Cox_interaction") %>%
  mutate(
    Predictor = case_when(
      term %in% c("rv_3798", "rv_5875", "rv_7308") ~ "Rare variant carrier status",
      term == "prs_std" ~ "PRS_CHD",
      grepl("prs_std", term) & grepl("rv_", term) ~ "PRS × rare variant interaction",
      TRUE ~ term
    )
  ) %>%
  select(
    Variant, Gene, Trait, N, Events,
    Predictor, `HR (95% CI)`, `P-value`,
    HR, CI_low, CI_high, p
  )


## ----------------------------
## 16. Print outputs
## ----------------------------
print(tab_prs_carrier)
print(tab_carrier_prs_logistic)
print(tab_cox_comparison)
print(tab_interaction)


## ----------------------------
## 17. Save outputs
## ----------------------------
fwrite(
  tab_prs_carrier,
  file = output_lm_table,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)

fwrite(
  tab_carrier_prs_logistic,
  file = output_glm_table,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)

fwrite(
  tab_cox_comparison,
  file = output_comparison_table,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)

fwrite(
  tab_interaction,
  file = output_interaction_table,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)