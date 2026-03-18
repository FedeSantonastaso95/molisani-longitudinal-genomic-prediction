############################################################
# Breast cancer lifetime risk analysis up to age 80
#
# Important note:
# The total number of breast cancer cases included in this analysis does not
# match the overall number of 474 prevalent + incident BC cases available in
# the full dataset.
#
# This is expected for two reasons:
#   1. The present analysis is restricted to events observed up to age 80,
#      so women whose event age exceeds 80 are excluded.
#   2. Some women with prevalent breast cancer at baseline do not have a
#      recorded age at onset and therefore cannot be included in a survival
#      analysis based on age as the time scale.
#
# As a result, the analysis dataset used here includes:
#   - 435 breast cancer cases
#   - 9,191 censored women
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
# Outcome coding used for Surv():
#   - 1 = censored
#   - 2 = breast cancer case
#
# Age_breast represents the age at event or censoring:
#   - for prevalent BC cases: age at onset
#   - for incident BC cases: age at follow-up diagnosis
#   - for censored women: age at end of follow-up
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
## 3. Restrict to women and to the age window <= 80 years
## ----------------------------
df_bc <- df2 %>%
  filter(
    Sex == 0,
    !is.na(Age_breast),
    !is.na(breast_Bs_final),
    Age_breast <= 80
  )


## ----------------------------
## 4. Function: PRS quintiles + HRs + counts
## ----------------------------
make_quintile_table <- function(data, pred) {
  
  # Create PRS quintiles
  qs <- quantile(data[[pred]], probs = c(.2, .4, .6, .8), na.rm = TRUE)
  brks <- c(min(data[[pred]], na.rm = TRUE), qs, max(data[[pred]], na.rm = TRUE))
  cat_col <- paste0(pred, "_Q")
  
  d <- data %>%
    filter(!is.na(.data[[pred]])) %>%
    mutate(
      !!cat_col := cut(
        .data[[pred]],
        breaks = brks,
        labels = paste0("Q", 1:5),
        include.lowest = TRUE
      )
    )
  
  # Counts by quintile
  counts <- d %>%
    group_by(.data[[cat_col]]) %>%
    summarise(
      N_total = n(),
      N_cases = sum(breast_Bs_final == 2),
      N_censored = sum(breast_Bs_final == 1),
      .groups = "drop"
    ) %>%
    rename(quintile = !!cat_col)
  
  # Cox model with Q1 as reference
  fit <- coxph(
    as.formula(paste0("Surv(Age_breast, breast_Bs_final) ~ Age + ", cat_col)),
    data = d
  )
  
  s <- summary(fit)
  rows_q <- grepl(cat_col, rownames(s$coefficients))
  
  hr <- s$conf.int[rows_q, "exp(coef)"]
  lo <- s$conf.int[rows_q, "lower .95"]
  hi <- s$conf.int[rows_q, "upper .95"]
  p  <- s$coefficients[rows_q, "Pr(>|z|)"]
  qn <- sub(".*Q", "Q", rownames(s$coefficients)[rows_q])
  
  hr_tab <- data.frame(
    predictor = pred,
    model = "Quintiles",
    quintile = c("Q1", qn),
    HR = c(1, hr),
    low = c(1, lo),
    high = c(1, hi),
    p = c(NA, p)
  )
  
  out <- hr_tab %>%
    left_join(counts, by = "quintile") %>%
    mutate(
      HR_CI = ifelse(
        quintile == "Q1",
        "1.00 [ref]",
        sprintf("%.2f (%.2f–%.2f)", HR, low, high)
      ),
      p = ifelse(
        is.na(p),
        NA,
        gsub("e", "×10^", formatC(p, format = "e", digits = 2))
      )
    ) %>%
    select(predictor, model, quintile, HR_CI, p, N_total, N_cases, N_censored)
  
  out
}


## ----------------------------
## 5. Function: continuous PRS (+1 SD) + counts
## ----------------------------
make_continuous_table <- function(data, pred) {
  
  d <- data %>%
    filter(!is.na(.data[[pred]]))
  
  d[[paste0(pred, "_z")]] <- scale(d[[pred]])
  
  fit <- coxph(
    as.formula(
      paste0("Surv(Age_breast, breast_Bs_final) ~ Age + ", pred, "_z")
    ),
    data = d
  )
  
  s <- summary(fit)
  r <- which(rownames(s$coefficients) == paste0(pred, "_z"))
  
  data.frame(
    predictor = pred,
    model = "Continuous (+1 SD)",
    quintile = NA,
    HR_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      s$conf.int[r, "exp(coef)"],
      s$conf.int[r, "lower .95"],
      s$conf.int[r, "upper .95"]
    ),
    p = gsub(
      "e", "×10^",
      formatC(s$coefficients[r, "Pr(>|z|)"], format = "e", digits = 2)
    ),
    N_total = nrow(d),
    N_cases = sum(d$breast_Bs_final == 2),
    N_censored = sum(d$breast_Bs_final == 1)
  )
}


## ----------------------------
## 6. Build final summary table
## ----------------------------
preds <- c(
  "PRED_eurW_breast_noBRCA1_noBRCA2",
  "PRED_eurW_breast"
)

tab_Q <- bind_rows(lapply(preds, make_quintile_table, data = df_bc))
tab_C <- bind_rows(lapply(preds, make_continuous_table, data = df_bc))

final_table <- bind_rows(tab_C, tab_Q)

rownames(final_table) <- NULL


## ----------------------------
## 7. Inspect final table
## ----------------------------
final_table


## ----------------------------
## 8. Optional export
## ----------------------------
write.table(
  final_table,
  "path/to/output/HRtable_quintiles_and_continuous_breast.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)