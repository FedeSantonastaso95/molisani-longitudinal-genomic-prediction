############################################################
# CHD composite score model
#
# Description:
# This script builds a composite polygenic score for incident CHD
# using the beta coefficients from a multivariable Cox model including:
#   - PRED_eurW_ApoB
#   - PRED_eurW_SBP
#   - PRED_eurW_CHD
#   - Age
#   - Sex
#
# Analysis design:
#   - incident-only CHD risk set
#   - baseline-prevalent CHD cases are excluded
#   - secondary CHD events are handled upstream by setting the follow-up
#     CHD event flag to NA when prevalent CHD is already present at baseline
#
# Output:
#   1. a table with coefficients from the multivariable Cox model
#   2. a composite score derived from the model beta coefficients
#   3. a Cox model testing the standardized composite score
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(survival)
library(pROC)
library(caret)
library(nricens)
library(Hmisc)
library(survivalROC)
library(data.table)
library(dplyr)
library(tidyr)
library(textshape)
library(car)
library(tibble)


## ----------------------------
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

outfile_model <- "path/to/output/2026.02.12_CompositeScore_CoxModel_secondVersion.txt"
outfile_composite <- "path/to/output/2026.02.12_CompositeScore_FinalResults.txt"


## ----------------------------
## 1. Read and preprocess data
## ----------------------------
df <- fread(
  infile,
  data.table = FALSE
) %>%
  mutate(
    PRED_eurW_ApoB_scaled = as.numeric(scale(PRED_eurW_ApoB)),
    PRED_eurW_SBP_scaled  = as.numeric(scale(PRED_eurW_SBP))
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
## 2. Identify CHD follow-up date and event columns
## ----------------------------
column_name <- "Chd0"

columns_disease <- grep(column_name, colnames(df), value = TRUE)

column_dataexit <- columns_disease[grepl("Dat", columns_disease)]
column_fnf <- columns_disease[grepl("Fn", columns_disease)]

if (length(column_fnf) > 1) {
  column_fnf <- column_fnf[!grepl("Dat", column_fnf)]
}
if (length(column_fnf) > 1) {
  column_fnf <- column_fnf[!grepl("Annipers", column_fnf)]
}

if (length(column_dataexit) != 1) {
  stop("Ambiguous CHD date column: ", paste(column_dataexit, collapse = ", "))
}
if (length(column_fnf) != 1) {
  stop("Ambiguous CHD event column: ", paste(column_fnf, collapse = ", "))
}


## ----------------------------
## 3. Build CHD survival dataset
## ----------------------------
df1 <- df %>%
  mutate(
    # Standard 0/1 coding for Surv():
    #   1 = incident CHD
    #   0 = no incident CHD
    status = case_when(
      .data[[column_fnf]] == 1 ~ 1,
      .data[[column_fnf]] == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

df1$Recruitment_date <- as.Date(
  gsub("/", "-", df1$Recruitment_date),
  format = "%d-%m-%Y"
)

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


## ----------------------------
## 4. Standardize CHD PRS
## ----------------------------
df3 <- df2 %>%
  mutate(
    PRED_eurW_CHD_scaled = as.numeric(scale(PRED_eurW_CHD))
  )


## ----------------------------
## 5. Fit multivariable Cox model for composite score derivation
## ----------------------------
#
# Model includes:
#   - ApoB PRS
#   - SBP PRS
#   - CHD PRS
#   - Age
#   - Sex
df_model <- df3 %>%
  filter(
    !is.na(time),
    !is.na(status),
    !is.na(PRED_eurW_ApoB_scaled),
    !is.na(PRED_eurW_SBP_scaled),
    !is.na(PRED_eurW_CHD_scaled),
    !is.na(Age),
    !is.na(Sex)
  )

res_cox_multivariable <- coxph(
  Surv(time, status) ~
    PRED_eurW_ApoB_scaled +
    PRED_eurW_SBP_scaled +
    PRED_eurW_CHD_scaled +
    Age +
    factor(Sex),
  data = df_model
)


## ----------------------------
## 6. Extract model coefficients and confidence intervals
## ----------------------------
table_cox <- as.data.frame(summary(res_cox_multivariable)$coefficients)
colnames(table_cox) <- c("beta", "HR", "SE", "Z", "p")
table_cox <- tibble::rownames_to_column(table_cox, var = "variable")

table_cox <- table_cox %>%
  mutate(
    HR_CI_lower = exp(beta - 1.96 * SE),
    HR_CI_upper = exp(beta + 1.96 * SE)
  )

table_cox_final <- table_cox %>%
  select(variable, beta, HR, HR_CI_lower, HR_CI_upper, SE, p)

fwrite(
  table_cox_final,
  file = outfile_model,
  row.names = FALSE,
  na = "NA",
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)


## ----------------------------
## 7. Extract beta values for composite score construction
## ----------------------------
beta_values <- table_cox_final %>%
  filter(variable %in% c(
    "PRED_eurW_ApoB_scaled",
    "PRED_eurW_SBP_scaled",
    "PRED_eurW_CHD_scaled",
    "Age",
    "factor(Sex)1"
  )) %>%
  select(variable, beta) %>%
  pivot_wider(names_from = variable, values_from = beta)


## ----------------------------
## 8. Construct and standardize composite score
## ----------------------------
df4 <- df3 %>%
  mutate(
    combined_score =
      (PRED_eurW_ApoB_scaled * beta_values$PRED_eurW_ApoB_scaled) +
      (PRED_eurW_SBP_scaled  * beta_values$PRED_eurW_SBP_scaled) +
      (PRED_eurW_CHD_scaled  * beta_values$PRED_eurW_CHD_scaled) +
      (as.numeric(Sex)       * beta_values$`factor(Sex)1`)
  ) %>%
  mutate(
    combined_score_scaled = as.numeric(scale(combined_score))
  )


## ----------------------------
## 9. Test composite score in a Cox model
## ----------------------------
df_composite <- df4 %>%
  filter(
    !is.na(time),
    !is.na(status),
    !is.na(combined_score_scaled)
  )

res_cox_composite <- coxph(
  Surv(time, status) ~ combined_score_scaled,
  data = df_composite
)

summary_cox <- summary(res_cox_composite)

results <- data.frame(
  variable = rownames(summary_cox$coefficients),
  HR       = summary_cox$conf.int[, "exp(coef)"],
  CI_lower = summary_cox$conf.int[, "lower .95"],
  CI_upper = summary_cox$conf.int[, "upper .95"],
  p_value  = summary_cox$coefficients[, "Pr(>|z|)"]
)

print(results)

fwrite(
  results,
  file = outfile_composite,
  row.names = FALSE,
  na = "NA",
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)