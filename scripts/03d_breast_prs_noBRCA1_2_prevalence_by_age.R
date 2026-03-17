############################################################
# Breast cancer cumulative risk summary up to ages 40 and 50
#
# Important note:
# The total number of breast cancer cases included in this analysis does not
# match the overall number of 474 prevalent + incident BC cases available in
# the full dataset.
#
# This is expected for two reasons:
#   1. The present analysis is restricted to women with evaluable event or
#      censoring age up to 80 years, so women whose event age exceeds 80 are
#      excluded.
#   2. Some women with prevalent breast cancer at baseline do not have a
#      recorded age at onset and therefore cannot be included in analyses
#      based on age as the time scale.
#
# As a consequence, the analysis dataset used here includes a reduced number
# of breast cancer cases compared with the overall case count.
#
# This script:
#   1. defines the breast cancer survival dataset in women only
#   2. creates PRS quintiles for the breast PRS excluding BRCA1/BRCA2
#   3. estimates the cumulative proportion of breast cancer cases observed
#      by age 40 and by age 50 within each PRS quintile
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(survival)


## ----------------------------
## 1. Read input data
## ----------------------------
df <- fread(
  "path/to/phenotype_dataset.txt",
  data.table = FALSE
) 


## ----------------------------
## 2. Recode breast cancer outcome
## ----------------------------
#
# Outcome coding:
#   - 2 = breast cancer case
#   - 1 = censored
#
# Age_breast is used as the age-at-event or age-at-censoring:
#   - prevalent cases: reported age at onset
#   - incident cases: age at diagnosis during follow-up
#   - censored women: age at end of follow-up
df2 <- df %>%
  mutate(
    breast_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 & Mt1_Icd9 == 174,
      2,
      ifelse(Fnf_Bc_Fup2020 == 1, 2, 1)
    ),
    Age_breast = case_when(
      # Prevalent breast cancer at baseline
      breast_Bs_final == 2 &
        Mt_Diagnosis_New == 1 &
        Mt1_Icd9 == 174 ~ Age_Mt1_New,
      
      # Incident breast cancer during follow-up
      breast_Bs_final == 2 &
        Fnf_Bc_Fup2020 == 1 ~
        (as.numeric(
          as.Date(Dataexit_Bc_Fup2020, "%Y-%m-%d") -
            as.Date(Recruitment_date, "%d/%m/%Y")
        ) / 365.25) + Age,
      
      # Censored at end of follow-up
      breast_Bs_final == 1 ~
        (as.numeric(
          as.Date(Dataexit_Bc_Fup2020, "%Y-%m-%d") -
            as.Date(Recruitment_date, "%d/%m/%Y")
        ) / 365.25) + Age
    )
  )


## ----------------------------
## 3. Restrict to women and to age <= 80 years
## ----------------------------
df_bc <- df2 %>%
  filter(
    Sex == 0,
    !is.na(Age_breast),
    !is.na(breast_Bs_final),
    Age_breast <= 80
  )


## ----------------------------
## 4. Create PRS quintiles
## ----------------------------
#
# Quintiles are based on the breast PRS excluding BRCA1 and BRCA2.
df_bc <- df_bc %>%
  mutate(
    PRED_eurW_breast_noBRCA1_noBRCA2_category =
      ntile(PRED_eurW_breast_noBRCA1_noBRCA2, 5)
  )


## ----------------------------
## 5. Dataset censored at age 40
## ----------------------------
#
# event_40:
#   - 1 if breast cancer occurred by age 40
#   - 0 otherwise
#
# time_40:
#   - observed age, truncated at 40 years
df40 <- df_bc %>%
  mutate(
    event_40 = ifelse(breast_Bs_final == 2 & Age_breast <= 40, 1, 0),
    time_40  = pmin(Age_breast, 40)
  )

prev40 <- df40 %>%
  group_by(PRED_eurW_breast_noBRCA1_noBRCA2_category) %>%
  summarise(
    quintile = paste0("Q", first(PRED_eurW_breast_noBRCA1_noBRCA2_category)),
    N_total = n(),
    N_cases_40 = sum(event_40, na.rm = TRUE),
    prevalence_40 = N_cases_40 / N_total,
    .groups = "drop"
  )


## ----------------------------
## 6. Dataset censored at age 50
## ----------------------------
#
# event_50:
#   - 1 if breast cancer occurred by age 50
#   - 0 otherwise
#
# time_50:
#   - observed age, truncated at 50 years
df50 <- df_bc %>%
  mutate(
    event_50 = ifelse(breast_Bs_final == 2 & Age_breast <= 50, 1, 0),
    time_50  = pmin(Age_breast, 50)
  )

prev50 <- df50 %>%
  group_by(PRED_eurW_breast_noBRCA1_noBRCA2_category) %>%
  summarise(
    quintile = paste0("Q", first(PRED_eurW_breast_noBRCA1_noBRCA2_category)),
    N_total = n(),
    N_cases_50 = sum(event_50, na.rm = TRUE),
    prevalence_50 = N_cases_50 / N_total,
    .groups = "drop"
  )


## ----------------------------
## 7. Merge age-40 and age-50 summary tables
## ----------------------------
prevalence_table <- prev40 %>%
  full_join(prev50, by = c("quintile", "N_total"))

prevalence_table