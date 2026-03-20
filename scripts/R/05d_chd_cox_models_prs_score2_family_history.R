############################################################
# Supplementary Table S12
# Cox models for PRS_CHD, SCORE2, and family history
#
# Updated version:
#   - Cox models adjusted for Age, Sex, and PC1:PC5
#   - extract only target predictors for the final table
#
# Description:
# This script fits Cox proportional-hazards models to evaluate the
# independent and joint associations of:
#   - PRS for CHD
#   - SCORE2
#   - family history of CHD
# with incident coronary heart disease (CHD).
#
# Models fitted:
#   - Model 0: PRS_CHD + covariates
#   - Model 1: PRS_CHD + SCORE2 + covariates
#   - Model 2: PRS_CHD + family history + covariates
#   - Model 3: PRS_CHD + SCORE2 + family history + covariates
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
# Covariate adjustment:
#   - Age
#   - Sex
#   - PC1, PC2, PC3, PC4, PC5
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
## 0. Clean workspace
## ----------------------------
rm(list = ls())


## ----------------------------
## 1. Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(tibble)


## ----------------------------
## 2. Input / output paths
## ----------------------------
main_file   <- "path/to/main_longitudinal_dataset.txt"
score2_file <- "path/to/output/2026.02.18_SCORE2_computed.txt"
aux_file    <- "path/to/auxiliary_phenotype_dataset.txt"

output_file <- "path/to/output/2026.03.20_Supplementary_Table_S12_CHD_PRS_SCORE2_FH_adjusted.txt"


## ----------------------------
## 3. Helper: format p-values
## ----------------------------
format_pvalue <- function(p) {
  if (is.na(p)) return(NA_character_)
  out <- formatC(p, digits = 2, format = "e")
  out <- gsub("e\\+?(-?\\d+)", "×10^\\1", out)
  out
}


## ----------------------------
## 4. Helper: extract target predictors
## ----------------------------
extract_target_terms <- function(model, model_name, n_total, n_events) {
  
  coef_df <- as.data.frame(summary(model)$coefficients) %>%
    tibble::rownames_to_column(var = "variable")
  
  conf_df <- as.data.frame(summary(model)$conf.int) %>%
    tibble::rownames_to_column(var = "variable")
  
  out <- coef_df %>%
    select(variable, coef, `Pr(>|z|)`) %>%
    rename(
      beta = coef,
      p_value_num = `Pr(>|z|)`
    ) %>%
    left_join(
      conf_df %>%
        select(var, `exp(coef)`, `lower .95`, `upper .95`) %>%
        rename(
          HR = `exp(coef)`,
          CI_low = `lower .95`,
          CI_high = `upper .95`
        ),
      by = c("variable" = "var")
    ) %>%
    filter(variable %in% c(
      "PRS_CHD",
      "SCORE2",
      "factor(Fam_Chd)1"
    )) %>%
    mutate(
      Predictor = case_when(
        variable == "PRS_CHD" ~ "PRS_CHD",
        variable == "SCORE2" ~ "SCORE2",
        variable == "factor(Fam_Chd)1" ~ "Family history",
        TRUE ~ variable
      ),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      `P-value` = vapply(p_value_num, format_pvalue, character(1)),
      Model = model_name,
      N_total = n_total,
      N_events = n_events
    ) %>%
    select(
      Model,
      Predictor,
      `HR (95% CI)`,
      `P-value`,
      N_total,
      N_events,
      HR,
      CI_low,
      CI_high,
      p_value_num
    )
  
  out
}


## ----------------------------
## 5. Load main dataset
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
## 6. Load SCORE2
## ----------------------------
score2_df <- fread(
  score2_file,
  data.table = FALSE
)


## ----------------------------
## 7. Load auxiliary phenotype data
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
## 8. Merge data sources
## ----------------------------
df <- df_main %>%
  left_join(score2_df, by = c("FID", "Idth_Ms2022")) %>%
  left_join(aux_df, by = "Idth_Ms2022")


## ----------------------------
## 9. Define analysis sample
## ----------------------------
df_analysis <- df %>%
  # Exclude prevalent CHD
  filter(!(Chd_Bs_New %in% c(1, 2))) %>%
  # Exclude prevalent stroke
  filter(!(Cerebro_Bs_New %in% c(1, 2))) %>%
  # Exclude diabetes
  filter(T2d_Arw != 1 | is.na(T2d_Arw)) %>%
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
    ),
    Age = as.numeric(Age),
    Sex = as.factor(Sex),
    PC1 = as.numeric(PC1),
    PC2 = as.numeric(PC2),
    PC3 = as.numeric(PC3),
    PC4 = as.numeric(PC4),
    PC5 = as.numeric(PC5)
  ) %>%
  filter(
    !is.na(time),
    !is.na(event)
  )


## ----------------------------
## 10. Use common complete-case dataset
## ----------------------------
#
# Using the same complete-case dataset across all four models ensures
# directly comparable estimates.
df_cc <- df_analysis %>%
  select(
    time, event,
    PRS_CHD, SCORE2, Fam_Chd,
    Age, Sex, PC1, PC2, PC3, PC4, PC5
  ) %>%
  tidyr::drop_na()

n_total <- nrow(df_cc)
n_events <- sum(df_cc$event == 1, na.rm = TRUE)

cat("N included in models:", n_total, "\n")
cat("N events:", n_events, "\n")


## ----------------------------
## 11. Fit Cox models
## ----------------------------
fit_model_0 <- coxph(
  Surv(time, event) ~ PRS_CHD + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
  data = df_cc
)

fit_model_1 <- coxph(
  Surv(time, event) ~ PRS_CHD + SCORE2 + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
  data = df_cc
)

fit_model_2 <- coxph(
  Surv(time, event) ~ PRS_CHD + factor(Fam_Chd) + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
  data = df_cc
)

fit_model_3 <- coxph(
  Surv(time, event) ~ PRS_CHD + SCORE2 + factor(Fam_Chd) + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
  data = df_cc
)


## ----------------------------
## 12. Extract model results
## ----------------------------
models <- list(
  "Model 0: PRS_CHD + covariates" = fit_model_0,
  "Model 1: PRS_CHD + SCORE2 + covariates" = fit_model_1,
  "Model 2: PRS_CHD + family history + covariates" = fit_model_2,
  "Model 3: PRS_CHD + SCORE2 + family history + covariates" = fit_model_3
)

final_table <- bind_rows(
  lapply(names(models), function(model_name) {
    extract_target_terms(
      model = models[[model_name]],
      model_name = model_name,
      n_total = n_total,
      n_events = n_events
    )
  })
)

print(final_table)


## ----------------------------
## 13. Save output
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