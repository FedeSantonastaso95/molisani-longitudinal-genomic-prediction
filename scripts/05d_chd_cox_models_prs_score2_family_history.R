############################################################
# Supplementary Table S12
# Cox models for PRS_CHD, SCORE2, and family history
#
# Description:
# This script fits Cox proportional-hazards models to evaluate the
# independent and joint associations of:
#   - PRS for CHD
#   - SCORE2
#   - family history of CHD
# with incident coronary heart disease (CHD).
#
#
# Models fitted:
#   - Model 0: PRS_CHD only
#   - Model 1: PRS_CHD + SCORE2
#   - Model 2: PRS_CHD + family history
#   - Model 3: PRS_CHD + SCORE2 + family history
#
# Analysis restrictions:
#   - QC-passed participants only
#   - incident CHD only
#   - exclude prevalent CHD at baseline
#   - exclude prevalent stroke at baseline
#   - exclude diabetes
#   - keep age 40-69 years
#
# Continuous predictors:
#   - PRS_CHD is standardized to 1 SD
#   - SCORE2 is standardized to 1 SD
#
# Input:
#   - main longitudinal dataset
#   - SCORE2 file computed previously
#   - auxiliary phenotype file with family history and diabetes
#
# Output:
#   A tab-delimited table formatted for Supplementary Table S12.
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(survival)
library(tibble)


## ----------------------------
## Input / output paths
## ----------------------------
main_file   <- "path/to/main_longitudinal_dataset.txt"
score2_file <- "path/to/output/2026.02.18_SCORE2_computed.txt"
aux_file    <- "path/to/auxiliary_phenotype_dataset.txt"

output_file <- "path/to/output/2026.02.23_Supplementary_Table_S12_CHD_PRS_SCORE2_FH.txt"


## ----------------------------
## Helper: format p-values
## ----------------------------
format_pvalue <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) {
    out <- formatC(p, format = "e", digits = 2)
    out <- gsub("e", "×10^", out, fixed = TRUE)
    return(out)
  }
  formatC(p, format = "f", digits = 3)
}


## ----------------------------
## Helper: extract model results
## ----------------------------
extract_cox_terms <- function(fit, model_label, predictor_labels) {
  sm <- summary(fit)
  
  coef_tab <- as.data.frame(sm$coefficients) %>%
    tibble::rownames_to_column("term")
  
  ci_tab <- as.data.frame(sm$conf.int) %>%
    tibble::rownames_to_column("term")
  
  out <- coef_tab %>%
    left_join(ci_tab, by = "term") %>%
    filter(term %in% names(predictor_labels)) %>%
    transmute(
      Model = model_label,
      Predictor = unname(predictor_labels[term]),
      HR = `exp(coef)`,
      CI_low = `lower .95`,
      CI_high = `upper .95`,
      p_value_num = `Pr(>|z|)`
    ) %>%
    mutate(
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      `P-value` = vapply(p_value_num, format_pvalue, character(1))
    ) %>%
    select(Model, Predictor, `HR (95% CI)`, `P-value`, HR, CI_low, CI_high, p_value_num)
  
  out
}


## ----------------------------
## 1. Load main dataset
## ----------------------------
df_main <- fread(
  main_file,
  data.table = FALSE
) %>%
  mutate(
    # Baseline CHD coding:
    #   0 = no
    #   1 = yes, documented
    #   2 = yes, not documented
    #
    # For analysis, both 1 and 2 are treated as prevalent CHD.
    Chd_Bs_New_final = case_when(
      Chd_Bs_New == 0 ~ 0,
      Chd_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Secondary CHD events are not counted as incident events:
    # if CHD is already prevalent at baseline and follow-up CHD flag is 1,
    # the follow-up event flag is set to NA.
    Chd0_Fnof_Fup2020 = if_else(
      Chd_Bs_New_final == 1 & Chd0_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Chd0_Fnof_Fup2020)
    )
  ) %>%
  # Remove baseline-prevalent CHD from the non-incident comparison group
  filter(!(Chd0_Fnof_Fup2020 == 0 & Chd_Bs_New_final == 1))


## ----------------------------
## 2. Load SCORE2
## ----------------------------
score2_df <- fread(
  score2_file,
  data.table = FALSE
)


## ----------------------------
## 3. Load auxiliary phenotype data
## ----------------------------
aux_df <- fread(
  aux_file,
  data.table = FALSE
) %>%
  select(
    Idth_Ms2022,
    Fam_Chd,
    T2d_Arw
  )


## ----------------------------
## 4. Merge data sources
## ----------------------------
df <- df_main %>%
  left_join(score2_df, by = c("FID", "Idth_Ms2022")) %>%
  left_join(aux_df, by = "Idth_Ms2022")


## ----------------------------
## 5. Define analysis sample
## ----------------------------
df_analysis <- df %>%
  # Exclude prevalent CHD
  filter(!(Chd_Bs_New %in% c(1, 2))) %>%
  # Exclude prevalent stroke
  filter(!(Cerebro_Bs_New %in% c(1, 2))) %>%
  # Exclude diabetes
  filter(T2d_Arw != 1) %>%
  # Keep SCORE2 age range
  filter(Age >= 40 & Age <= 69) %>%
  transmute(
    Idth_Ms2022,
    FID,
    time = as.numeric(Annipers_Chd0_Fnof_Fup2020),
    event = as.integer(Chd0_Fnof_Fup2020),
    PRS_CHD = as.numeric(scale(PRED_eurW_CHD)),
    SCORE2 = as.numeric(scale(SCORE2_score)),
    Fam_Chd = case_when(
      Fam_Chd == 1 ~ 1,
      Fam_Chd == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(time),
    !is.na(event)
  )


## ----------------------------
## 6. Use common complete-case dataset
## ----------------------------
#
# Using the same complete-case dataset across all four models ensures
# directly comparable estimates.
df_cc <- df_analysis %>%
  select(time, event, PRS_CHD, SCORE2, Fam_Chd) %>%
  tidyr::drop_na()


## ----------------------------
## 7. Fit Cox models
## ----------------------------
fit_model_0 <- coxph(
  Surv(time, event) ~ PRS_CHD,
  data = df_cc
)

fit_model_1 <- coxph(
  Surv(time, event) ~ PRS_CHD + SCORE2,
  data = df_cc
)

fit_model_2 <- coxph(
  Surv(time, event) ~ PRS_CHD + Fam_Chd,
  data = df_cc
)

fit_model_3 <- coxph(
  Surv(time, event) ~ PRS_CHD + SCORE2 + Fam_Chd,
  data = df_cc
)


## ----------------------------
## 8. Extract model results
## ----------------------------
predictor_map <- c(
  PRS_CHD = "PRS_CHD",
  SCORE2 = "SCORE2",
  Fam_Chd = "Family history"
)

tab_model_0 <- extract_cox_terms(
  fit = fit_model_0,
  model_label = "Model 0",
  predictor_labels = predictor_map
)

tab_model_1 <- extract_cox_terms(
  fit = fit_model_1,
  model_label = "Model 1",
  predictor_labels = predictor_map
)

tab_model_2 <- extract_cox_terms(
  fit = fit_model_2,
  model_label = "Model 2",
  predictor_labels = predictor_map
)

tab_model_3 <- extract_cox_terms(
  fit = fit_model_3,
  model_label = "Model 3",
  predictor_labels = predictor_map
)


## ----------------------------
## 9. Combine final table
## ----------------------------
final_table <- bind_rows(
  tab_model_0,
  tab_model_1,
  tab_model_2,
  tab_model_3
)

print(final_table)


## ----------------------------
## 10. Save output
## ----------------------------
fwrite(
  final_table,
  file = output_file,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)