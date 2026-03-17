############################################################
# CHD combinatorial models with selected PRS and baseline traits
#
# Description:
# This script fits two Cox proportional hazards models for incident CHD:
#
#   1. Traits-only model
#      Includes only the polygenic predictors retained after LASSO selection
#
#   2. Full model
#      Includes the same polygenic predictors plus the corresponding
#      baseline phenotypic measurements, to evaluate whether measured
#      baseline traits capture part of the genetic risk
#
# Predictors included:
#   - PRS_CHD
#   - PGI_ApoB
#   - PGI_SBP
#
# Corresponding baseline traits:
#   - ApoB  -> Bio_Apo_B
#   - SBP   -> SBP
#
# Analysis design:
#   - incident-only CHD risk set
#   - baseline-prevalent CHD cases excluded
#   - secondary CHD events handled upstream by setting the follow-up
#     CHD event flag to NA when prevalent CHD is already present at baseline
#
# Model adjustment:
#   - Age
#   - Sex
#   - PC1-PC5
#
# Output:
#   - one table for the traits-only model
#   - one table for the full model
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(survival)
library(dplyr)
library(rlang)
library(data.table)
library(tibble)


## ----------------------------
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

outfile_full <- "path/to/output/2026.02.12_CombinatorialModel_Full_secondVersion.txt"
outfile_traits_only <- "path/to/output/2026.02.12_CombinatorialModel_TraitsOnly_secondVersion.txt"


## ----------------------------
## 1. Read and preprocess data
## ----------------------------
df <- fread(
  infile,
  data.table = FALSE
) %>%
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
    )
  ) %>%
  mutate(
    Fnf_cancer_FUP2020 = ifelse(is.na(Fnf_cancer_FUP2020), 0, Fnf_cancer_FUP2020)
  ) %>%
  mutate(
    # Baseline CHD coding:
    #   0 = No
    #   1 = Yes, documented
    #   2 = Yes, not documented
    #
    # In the analysis, both 1 and 2 are treated as prevalent CHD.
    Chd_Bs_New_final = case_when(
      Chd_Bs_New == 0 ~ 0,
      Chd_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Secondary CHD events are not counted as incident events:
    # if CHD is already prevalent at baseline and the follow-up event flag is 1,
    # the follow-up flag is set to NA.
    Chd0_Fnof_Fup2020 = if_else(
      Chd_Bs_New_final == 1 & Chd0_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Chd0_Fnof_Fup2020)
    )
  ) %>%
  # Remove baseline-prevalent CHD from the non-incident comparison group
  filter(!(Chd0_Fnof_Fup2020 == 0 & Chd_Bs_New_final == 1))


## ----------------------------
## 2. Define predictors and corresponding baseline phenotypes
## ----------------------------
df_legend <- data.frame(
  traits = c("CHD", "ApoB", "SBP"),
  pheno  = c(NA, "Bio_Apo_B", "SBP"),
  stringsAsFactors = FALSE
)


## ----------------------------
## 3. Fit CHD combinatorial models
## ----------------------------
fit_chd_combinatorial_models <- function(dat, endpoint_pattern = "Chd0") {
  
  # ------------------------------------------
  # Identify endpoint-specific date/event columns
  # ------------------------------------------
  columns_disease <- grep(endpoint_pattern, colnames(dat), value = TRUE)
  
  column_dataexit <- columns_disease[grepl("Dat", columns_disease)]
  column_fnf      <- columns_disease[grepl("Fn", columns_disease)]
  
  if (length(column_fnf) > 1) {
    column_fnf <- column_fnf[!grepl("Dat|Annipers", column_fnf)]
  }
  
  if (length(column_dataexit) != 1) {
    stop("Ambiguous CHD date column: ", paste(column_dataexit, collapse = ", "))
  }
  
  if (length(column_fnf) != 1) {
    stop("Ambiguous CHD event column: ", paste(column_fnf, collapse = ", "))
  }
  
  # ------------------------------------------
  # Build binary event indicator
  # ------------------------------------------
  df1 <- dat %>%
    mutate(
      status = case_when(
        .data[[column_fnf]] == 1 ~ 1,
        .data[[column_fnf]] == 0 ~ 0,
        TRUE ~ NA_real_
      )
    )
  
  # ------------------------------------------
  # Process recruitment date and follow-up time
  # ------------------------------------------
  df1$Recruitment_date <- gsub("/", "-", df1$Recruitment_date)
  df1$Recruitment_date <- as.Date(df1$Recruitment_date, format = "%d-%m-%Y")
  
  df2 <- df1 %>%
    mutate(
      !!sym(column_dataexit) := as.Date(.data[[column_dataexit]]),
      !!sym(column_dataexit) := if_else(
        is.na(.data[[column_dataexit]]),
        as.Date("2020-12-31"),
        .data[[column_dataexit]]
      ),
      time = as.numeric(difftime(.data[[column_dataexit]], Recruitment_date, units = "days")) / 365.25,
      Sex = as.factor(Sex),
      Recruiting_centre = as.factor(Recruiting_centre)
    ) %>%
    filter(!is.na(time), time >= 0)
  
  # ------------------------------------------
  # Add standardized PRS columns
  # ------------------------------------------
  df3 <- df2
  
  for (tr in df_legend$traits) {
    prs_col <- paste0("PRED_eurW_", tr)
    prs_sc  <- paste0(prs_col, "_scaled")
    
    if (!(prs_col %in% names(df3))) {
      stop("Missing PRS column: ", prs_col)
    }
    
    df3[[prs_sc]] <- as.numeric(scale(df3[[prs_col]]))
  }
  
  # ------------------------------------------
  # Add standardized baseline phenotypes
  # ------------------------------------------
  pheno_names <- df_legend$pheno[!is.na(df_legend$pheno)]
  
  df4 <- df3 %>%
    select(
      time,
      status,
      Age,
      Sex,
      PC1,
      PC2,
      PC3,
      PC4,
      PC5,
      all_of(paste0("PRED_eurW_", df_legend$traits, "_scaled")),
      all_of(pheno_names)
    )
  
  # Standardize observed baseline phenotypes and rename them
  for (i in seq_len(nrow(df_legend))) {
    if (!is.na(df_legend$pheno[i])) {
      pheno_col <- df_legend$pheno[i]
      new_name  <- paste0(df_legend$traits[i], "_scaled")
      
      df4[[new_name]] <- as.numeric(scale(df4[[pheno_col]]))
    }
  }
  
  # Complete-case dataset for modelling
  df4 <- df4 %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  # ------------------------------------------
  # Model formulas
  # ------------------------------------------
  
  # PRS-only model
  prs_scaled_cols <- paste0("PRED_eurW_", df_legend$traits, "_scaled")
  
  formula_traits_only <- as.formula(
    paste(
      "Surv(time, status) ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 +",
      paste(prs_scaled_cols, collapse = " + ")
    )
  )
  
  # Full model = PRS + baseline phenotypes
  baseline_scaled_cols <- c("ApoB_scaled", "SBP_scaled")
  full_predictors <- c(prs_scaled_cols, baseline_scaled_cols)
  
  formula_full <- as.formula(
    paste(
      "Surv(time, status) ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 +",
      paste(full_predictors, collapse = " + ")
    )
  )
  
  # ------------------------------------------
  # Fit models
  # ------------------------------------------
  fit_traits_only <- coxph(formula_traits_only, data = df4)
  fit_full        <- coxph(formula_full, data = df4)
  
  # ------------------------------------------
  # Helper: convert Cox summary to table
  # ------------------------------------------
  summary_to_table <- function(cox_model) {
    
    coef_table <- as.data.frame(summary(cox_model)$coefficient)
    colnames(coef_table) <- c("beta", "HR", "SE", "Z", "pval")
    coef_table <- tibble::rownames_to_column(coef_table, var = "variable")
    
    conf_table <- as.data.frame(summary(cox_model)$conf.int) %>%
      dplyr::select(-`exp(-coef)`)
    colnames(conf_table) <- c("HR", "CI95_lower", "CI95_upper")
    conf_table <- tibble::rownames_to_column(conf_table, var = "variable") %>%
      dplyr::select(-HR)
    
    out <- coef_table %>%
      dplyr::select(-beta, -SE, -Z) %>%
      full_join(conf_table, by = "variable") %>%
      mutate(
        HR = round(HR, 2),
        CI95_lower = round(CI95_lower, 2),
        CI95_upper = round(CI95_upper, 2),
        pval_num = pval,
        pval = format(pval, scientific = TRUE, digits = 3),
        stars = case_when(
          pval_num < 0.001 ~ "***",
          pval_num < 0.01  ~ "**",
          pval_num < 0.05  ~ "*",
          pval_num < 0.1   ~ ".",
          TRUE ~ "n.s."
        )
      ) %>%
      dplyr::select(variable, HR, CI95_lower, CI95_upper, pval, stars)
    
    out
  }
  
  table_traits_only <- summary_to_table(fit_traits_only)
  table_full        <- summary_to_table(fit_full)
  
  list(
    table_traits_only = table_traits_only,
    table_full = table_full,
    fit_traits_only = fit_traits_only,
    fit_full = fit_full
  )
}


## ----------------------------
## 4. Run analysis
## ----------------------------
res <- fit_chd_combinatorial_models(dat = df, endpoint_pattern = "Chd0")


## ----------------------------
## 5. Save output tables
## ----------------------------
fwrite(
  res$table_full,
  outfile_full,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

fwrite(
  res$table_traits_only,
  outfile_traits_only,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)