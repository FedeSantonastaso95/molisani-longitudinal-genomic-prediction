# =========================================================
# In this script, you will find the Cox interaction models
# testing whether the effect of polygenic risk score (PRS)
# on incident CHD differs according to carrier status of
# selected rare variants.
#
# IMPORTANT:
# - Genotype dosages were exported with PLINK2 using `--export A`.
# - In this dataset, the exported dosage corresponds to the major allele.
# - Therefore, dosages are recoded so that:
#     0 = non-carrier of the minor allele
#     1 = heterozygous carrier
#     2 = homozygous carrier of the minor allele
#
# This script performs:
# - Cox models with additive effects of PRS and rare variant carrier status
# - Cox models including the PRS × rare variant interaction term
# - Likelihood ratio tests comparing additive vs interaction models
#
# Variants tested:
# - rs3798220   (LPA / Lp(a))
# - rs58757394  (SLC4A11 / SBP)
# - rs730882080 (LDLR / LDL)
#
# Outputs:
# - Console summaries of interaction models
# - A supplementary table with HRs, confidence intervals, and p-values
# =========================================================

# =========================================================
# Libraries
# =========================================================
library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggpubr)
library(ggtext)
library(scales)
library(grid)
library(tibble)

# =========================================================
# Load carrier file
# =========================================================
df_carr <- fread("/path/to/carriers_rarVar.raw")

# Round all dosage columns
df_carr[, 7:ncol(df_carr)] <- lapply(df_carr[, 7:ncol(df_carr)], round)

# PLINK dosages correspond to the major allele in this dataset.
# Recode so that:
#   0 = non-carrier of the rare allele
#   1 = heterozygous carrier
#   2 = homozygous carrier of the rare allele
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

# Clean variant column names by removing the suffix after the last underscore
old_names <- names(df_carr)[7:ncol(df_carr)]
new_names <- sub("_(?!.*_).*$", "", old_names, perl = TRUE)
names(df_carr)[7:ncol(df_carr)] <- new_names

# =========================================================
# Load rare variants and SNP list
# =========================================================
var <- fread("/path/to/listVariants_Burden_ALL_final.txt")
snplist <- fread("/path/to/imputed_genotypes.snplist")

var1 <- var %>%
  left_join(
    snplist,
    by = c(
      "GENPOS"  = "V4",
      "CHROM"   = "V1",
      "ALLELE0" = "V5",
      "ALLELE1" = "V6"
    )
  ) %>%
  dplyr::rename(bgen_ID = V2) %>%
  dplyr::select(-V3)

setDT(var1)
unique_combos <- unique(var1[, .(gene_name, trait)])

# =========================================================
# Load sample file
# =========================================================
sample <- fread("/path/to/imputed_genotypes.sample")

# Remove the first line of the PLINK/BGEN sample file if it contains
# the standard dummy row
sample <- sample[-1, ]

# =========================================================
# Load phenotype data
# =========================================================
df <- fread(
  "/path/to/phenotypes_covariates_prs.txt",
  data.table = FALSE
) %>%
  left_join(df_carr, by = c("FID" = "FID", "IID" = "IID"))

# Restrict to sequence-imputed samples
df1 <- sample %>%
  left_join(df, by = c("ID_1" = "FID"))

# Keep only variables needed for interaction analyses
df1 <- df1 %>%
  select(
    IID,
    Chd0_Fnof_Fup2020,
    Annipers_Chd0_Fnof_Fup2020,
    rs3798220,
    rs58757394,
    rs730882080,
    PRED_eurW_CHD,
    Bio_Lp_A,
    SBP,
    LDL
  )

# =========================================================
# Prepare data for Cox models
# =========================================================
df_cox <- df1 %>%
  mutate(
    prs_std = as.numeric(scale(PRED_eurW_CHD)),
    rv_3798 = ifelse(rs3798220   == 1, 1, 0),
    rv_5875 = ifelse(rs58757394  == 1, 1, 0),
    rv_7308 = ifelse(rs730882080 == 1, 1, 0)
  ) %>%
  filter(
    !is.na(prs_std),
    !is.na(Annipers_Chd0_Fnof_Fup2020),
    !is.na(Chd0_Fnof_Fup2020)
  )

# =========================================================
# Fit additive and interaction Cox models
# =========================================================

# rs3798220 (LPA)
df_cox_3798 <- df_cox %>%
  filter(!is.na(rv_3798))

cox_3798_add <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    rv_3798 + prs_std,
  data = df_cox_3798
)

cox_3798_int <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    rv_3798 * prs_std,
  data = df_cox_3798
)

cat("\n===== rs3798220 (LPA) =====\n")
print(summary(cox_3798_int))
cat("\nLRT interaction (rs3798220 × PRS):\n")
print(anova(cox_3798_add, cox_3798_int, test = "LRT"))

# rs58757394 (SLC4A11)
df_cox_5875 <- df_cox %>%
  filter(!is.na(rv_5875))

cox_5875_add <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    rv_5875 + prs_std,
  data = df_cox_5875
)

cox_5875_int <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    rv_5875 * prs_std,
  data = df_cox_5875
)

cat("\n===== rs58757394 (SLC4A11) =====\n")
print(summary(cox_5875_int))
cat("\nLRT interaction (rs58757394 × PRS):\n")
print(anova(cox_5875_add, cox_5875_int, test = "LRT"))

# rs730882080 (LDLR)
df_cox_7308 <- df_cox %>%
  filter(!is.na(rv_7308))

cox_7308_add <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    rv_7308 + prs_std,
  data = df_cox_7308
)

cox_7308_int <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    rv_7308 * prs_std,
  data = df_cox_7308
)

cat("\n===== rs730882080 (LDLR) =====\n")
print(summary(cox_7308_int))
cat("\nLRT interaction (rs730882080 × PRS):\n")
print(anova(cox_7308_add, cox_7308_int, test = "LRT"))

# =========================================================
# Helper functions for supplementary table
# =========================================================
extract_term_full <- function(fit, term) {
  s <- summary(fit)
  tibble(
    HR      = s$coef[term, "exp(coef)"],
    CI_low  = s$conf.int[term, "lower .95"],
    CI_high = s$conf.int[term, "upper .95"],
    p       = s$coef[term, "Pr(>|z|)"]
  )
}

fmt_hrci_p <- function(HR, lo, hi, p) {
  paste0(
    sprintf("%.2f", HR),
    " (", sprintf("%.2f", lo), "–", sprintf("%.2f", hi), ")",
    "; p=", format.pval(p, digits = 3, eps = 1e-3)
  )
}

make_sup_row_full <- function(variant, gene, trait, fit_add, fit_int, rv_term) {
  
  n      <- fit_add$n
  events <- fit_add$nevent
  
  rv_add  <- extract_term_full(fit_add, rv_term)
  prs_add <- extract_term_full(fit_add, "prs_std")
  
  int_term <- paste0(rv_term, ":prs_std")
  rv_int   <- extract_term_full(fit_int, rv_term)
  prs_int  <- extract_term_full(fit_int, "prs_std")
  inter    <- extract_term_full(fit_int, int_term)
  
  lrt_p <- anova(fit_add, fit_int, test = "LRT")[2, "Pr(>|Chi|)"]
  
  tibble(
    Variant = variant,
    Gene    = gene,
    Trait   = trait,
    N       = n,
    Events  = events,
    
    RV_add_HR      = rv_add$HR,
    RV_add_CI_low  = rv_add$CI_low,
    RV_add_CI_high = rv_add$CI_high,
    RV_add_p       = rv_add$p,
    
    PRS_add_HR      = prs_add$HR,
    PRS_add_CI_low  = prs_add$CI_low,
    PRS_add_CI_high = prs_add$CI_high,
    PRS_add_p       = prs_add$p,
    
    RV_int_HR      = rv_int$HR,
    RV_int_CI_low  = rv_int$CI_low,
    RV_int_CI_high = rv_int$CI_high,
    RV_int_p       = rv_int$p,
    
    PRS_int_HR      = prs_int$HR,
    PRS_int_CI_low  = prs_int$CI_low,
    PRS_int_CI_high = prs_int$CI_high,
    PRS_int_p       = prs_int$p,
    
    INT_HR      = inter$HR,
    INT_CI_low  = inter$CI_low,
    INT_CI_high = inter$CI_high,
    INT_p       = inter$p,
    
    LRT_INT_p = lrt_p
  )
}

# =========================================================
# Build supplementary table
# =========================================================
sup_tab <- bind_rows(
  make_sup_row_full("rs3798220",   "LPA",     "Lp(a)", cox_3798_add, cox_3798_int, "rv_3798"),
  make_sup_row_full("rs58757394",  "SLC4A11", "SBP",   cox_5875_add, cox_5875_int, "rv_5875"),
  make_sup_row_full("rs730882080", "LDLR",    "LDL",   cox_7308_add, cox_7308_int, "rv_7308")
)

sup_int <- sup_tab %>%
  select(
    Variant, Gene, Trait,
    starts_with("RV_int_"),
    starts_with("PRS_int_"),
    starts_with("INT_"),
    LRT_INT_p
  )

write.csv(
  sup_int,
  "/path/to/survival_files/Supplementary_Table_RV_PRS_Cox.csv",
  row.names = FALSE,
  quote = TRUE
)