############################################################
# Incident CHD model evaluation
#
# Description:
# This script builds the analysis dataset for incident CHD prediction,
# computes the composite score, and evaluates multiple Cox models using:
#   - Harrell's C-index with bootstrap confidence intervals
#   - delta C-index versus the SCORE2-only model
#   - time-dependent AUC at 10 years with bootstrap confidence intervals
#
# Additional feature:
#   SCORE2 is also residualized on Age and Sex, and the residualized
#   version is tested in alternative prediction models.
#
# Analysis design:
#   - incident-only CHD risk set
#   - baseline-prevalent CHD cases are excluded
#   - prevalent stroke is excluded
#   - participants with diabetes are excluded
#   - age restricted to 40-69 years
#
# Output:
#   A final summary table including:
#   - model label
#   - C-index (95% CI)
#   - delta C-index vs SCORE2 (95% CI)
#   - AUC at 10 years (95% CI)
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(survival)
library(nricens)
library(Hmisc)
library(survivalROC)
library(data.table)
library(dplyr)
library(tidyr)
library(textshape)
library(car)
library(boot)
library(riskRegression)
library(tibble)


## ----------------------------
## Input / output
## ----------------------------
infile_main  <- "path/to/main_longitudinal_dataset.txt"
infile_score <- "path/to/2026.02.18_SCORE2_computed.txt"
infile_aux   <- "path/to/auxiliary_phenotype_dataset.txt"

outfile_eval <- "path/to/output/2026.02.21_IncidentCHD_Prediction_EvaluationMetrics_SCORE2resid.txt"


## ----------------------------
## 0. Build analysis dataset (df7)
## ----------------------------

# Main longitudinal dataset
df <- fread(
  infile_main,
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
    #   0 = no
    #   1 = yes, documented
    #   2 = yes, not documented
    #
    # In the analysis, both 1 and 2 are treated as prevalent CHD.
    Chd_Bs_New_final = ifelse(Chd_Bs_New == 0, 0, 1),
    
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

# Add SCORE2
score <- fread(infile_score, data.table = FALSE)

df <- df %>%
  left_join(score, by = c("FID", "Idth_Ms2022"))

# Add family history and diabetes/smoking
df_aux <- fread(infile_aux, data.table = FALSE) %>%
  select(
    Idth_Ms2022,
    Framingham98,
    Fam_Chd,
    Fam_Chd_Youth,
    Smoking_Tabacco,
    T2d_Arw
  )

df2 <- df %>%
  left_join(df_aux, by = "Idth_Ms2022")

# Primary prevention only + no prevalent stroke + no diabetes
df3 <- df2 %>%
  filter(!(Chd_Bs_New %in% c(1, 2))) %>%
  filter(!(Cerebro_Bs_New %in% c(1, 2))) %>%
  filter(T2d_Arw != 1)

# SCORE2 age range
df4 <- df3 %>%
  filter(Age >= 40 & Age <= 69)

# Select variables used in model evaluation
df5 <- df4 %>%
  select(
    Idth_Ms2022,
    FID,
    Age,
    Sex,
    PC1:PC5,
    Bio_Apo_B,
    SBP,
    PRED_eurW_CHD,
    PRED_eurW_SBP,
    PRED_eurW_ApoB,
    Fam_Chd,
    SCORE2_score,
    Chd0_Fnof_Fup2020,
    Annipers_Chd0_Fnof_Fup2020
  ) %>%
  mutate(
    PRED_eurW_SBP_scaled  = as.numeric(scale(PRED_eurW_SBP)),
    PRED_eurW_ApoB_scaled = as.numeric(scale(PRED_eurW_ApoB)),
    PRED_eurW_CHD_scaled  = as.numeric(scale(PRED_eurW_CHD))
  )

# Build composite score from multivariable Cox model
res_cox_composite_base <- coxph(
  Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
    PRED_eurW_ApoB_scaled +
    PRED_eurW_SBP_scaled +
    PRED_eurW_CHD_scaled +
    Age +
    factor(Sex),
  data = df5
)

table_cox_composite_base <- as.data.frame(summary(res_cox_composite_base)$coefficients)
colnames(table_cox_composite_base) <- c("beta", "HR", "SE", "Z", "p")
table_cox_composite_base <- tibble::rownames_to_column(table_cox_composite_base, var = "variable")

table_cox_composite_base <- table_cox_composite_base %>%
  mutate(
    HR_CI_lower = exp(beta - 1.96 * SE),
    HR_CI_upper = exp(beta + 1.96 * SE)
  ) %>%
  select(variable, beta, HR, HR_CI_lower, HR_CI_upper, SE, p)

beta_values <- table_cox_composite_base %>%
  filter(variable %in% c(
    "PRED_eurW_ApoB_scaled",
    "PRED_eurW_SBP_scaled",
    "PRED_eurW_CHD_scaled",
    "Age",
    "factor(Sex)1"
  )) %>%
  select(variable, beta) %>%
  pivot_wider(names_from = variable, values_from = beta)

df6 <- df5 %>%
  mutate(
    combined_score =
      (PRED_eurW_ApoB_scaled * beta_values$PRED_eurW_ApoB_scaled) +
      (PRED_eurW_SBP_scaled  * beta_values$PRED_eurW_SBP_scaled) +
      (PRED_eurW_CHD_scaled  * beta_values$PRED_eurW_CHD_scaled) +
      as.numeric(Sex)        * beta_values$`factor(Sex)1`
  ) %>%
  mutate(
    combined_score_scaled = as.numeric(scale(combined_score))
  )

# Create combined-score quintiles
breaks <- quantile(
  df6$combined_score,
  probs = seq(0, 1, 0.2),
  na.rm = TRUE
)
breaks <- unique(breaks)

df7 <- df6 %>%
  mutate(
    combined_score_quintile = cut(
      combined_score,
      breaks = breaks,
      labels = FALSE,
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    combined_score_quintile = factor(combined_score_quintile),
    combined_score_scaled = as.numeric(scale(combined_score))
  )


## ----------------------------
## 1. Define time and event
## ----------------------------
df_eval <- df7 %>%
  mutate(
    time  = Annipers_Chd0_Fnof_Fup2020,
    event = Chd0_Fnof_Fup2020
  ) %>%
  as.data.frame()

base_cov <- c("PC1", "PC2", "PC3", "PC4", "PC5")


## ----------------------------
## 2. Helper: robust numeric conversion
## ----------------------------
force_numeric <- function(x) {
  if (is.matrix(x) && ncol(x) == 1) x <- x[, 1]
  if (is.list(x)) x <- unlist(x)
  as.numeric(x)
}


## ----------------------------
## 3. Coerce key variables to numeric
## ----------------------------
key_vars <- intersect(
  c(
    "time", "event", "Age", "Sex", "Fam_Chd",
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PRED_eurW_CHD_scaled", "SCORE2_score", "combined_score_scaled"
  ),
  names(df_eval)
)

for (v in key_vars) {
  df_eval[[v]] <- force_numeric(df_eval[[v]])
}

df_eval$event <- as.integer(df_eval$event)


## ----------------------------
## 4. Residualize SCORE2 on Age and Sex
## ----------------------------
lm_score2_as <- lm(
  SCORE2_score ~ Age + Sex,
  data = df_eval,
  na.action = na.exclude
)

df_eval <- df_eval %>%
  mutate(
    SCORE2_resid = as.numeric(resid(lm_score2_as)),
    SCORE2_resid_scaled = as.numeric(scale(SCORE2_resid))
  )


## ----------------------------
## 5. Define model sets
## ----------------------------
model_vars_surv <- list(
  cidx_1  = c("time", "event", base_cov, "Fam_Chd"),
  cidx_2  = c("time", "event", base_cov, "PRED_eurW_CHD_scaled"),
  cidx_3  = c("time", "event", base_cov, "Fam_Chd", "PRED_eurW_CHD_scaled"),
  cidx_4  = c("time", "event", base_cov, "PRED_eurW_CHD_scaled", "SCORE2_score"),
  cidx_5  = c("time", "event", base_cov, "SCORE2_score"),
  cidx_6  = c("time", "event", base_cov, "combined_score_scaled"),
  cidx_7  = c("time", "event", base_cov, "Fam_Chd", "PRED_eurW_CHD_scaled", "SCORE2_score"),
  cidx_8  = c("time", "event", base_cov, "PRED_eurW_CHD_scaled", "Sex", "Age"),
  cidx_9  = c("time", "event", base_cov, "SCORE2_resid_scaled"),
  cidx_10 = c("time", "event", base_cov, "PRED_eurW_CHD_scaled", "SCORE2_resid_scaled"),
  cidx_11 = c("time", "event", base_cov, "Fam_Chd", "PRED_eurW_CHD_scaled", "SCORE2_resid_scaled")
)


## ----------------------------
## 6. Build common complete-case dataset
## ----------------------------
needed <- unique(unlist(model_vars_surv))

df_common <- df_eval %>%
  dplyr::select(all_of(needed)) %>%
  tidyr::drop_na() %>%
  as.data.frame()


## ----------------------------
## 7. Flatten matrix columns if needed
## ----------------------------
mat_cols <- names(df_common)[sapply(df_common, is.matrix)]

if (length(mat_cols) > 0) {
  for (cn in mat_cols) {
    m <- df_common[[cn]]
    df_common[[cn]] <- NULL
    m <- as.data.frame(m)
    names(m) <- paste0(cn, "_", seq_len(ncol(m)))
    df_common <- cbind(df_common, m)
  }
}

nm <- names(df_common)
to_fix <- grep("_scaled_1$", nm, value = TRUE)

if (length(to_fix) > 0) {
  for (old in to_fix) {
    new <- sub("_1$", "", old)
    names(df_common)[names(df_common) == old] <- new
  }
}


## ----------------------------
## 8. Force numeric again after flattening
## ----------------------------
num_candidates <- intersect(
  c(
    "time", "event", "Age", "Sex", "Fam_Chd",
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PRED_eurW_CHD_scaled", "SCORE2_score", "combined_score_scaled",
    "SCORE2_resid", "SCORE2_resid_scaled"
  ),
  names(df_common)
)

for (v in num_candidates) {
  df_common[[v]] <- force_numeric(df_common[[v]])
}

df_common$event <- as.integer(df_common$event)


## ----------------------------
## 9. Helpers: Cox fit and C-index
## ----------------------------
fit_cox_from_vars <- function(data, vars) {
  predictors <- setdiff(vars, c("time", "event"))
  f <- as.formula(
    paste0("Surv(time, event) ~ ", paste(predictors, collapse = " + "))
  )
  coxph(f, data = data, x = TRUE, y = TRUE)
}

cindex_from_cox <- function(fit) {
  as.numeric(summary(fit)$concordance[1])
}


## ----------------------------
## 10. Fit all models
## ----------------------------
models_cox <- lapply(names(model_vars_surv), function(m) {
  fit_cox_from_vars(df_common, model_vars_surv[[m]])
})
names(models_cox) <- names(model_vars_surv)


## ----------------------------
## 11. C-index with bootstrap confidence intervals
## ----------------------------
R_c <- 200
set.seed(1)

boot_cindex_refit <- function(data, idx, vars) {
  d <- data[idx, , drop = FALSE]
  out <- tryCatch({
    fit <- fit_cox_from_vars(d, vars)
    cindex_from_cox(fit)
  }, error = function(e) NA_real_)
  
  if (length(out) != 1) out <- NA_real_
  out
}

C_boot <- bind_rows(lapply(names(model_vars_surv), function(m) {
  b <- boot(
    df_common,
    statistic = function(dat, idx) boot_cindex_refit(dat, idx, model_vars_surv[[m]]),
    R = R_c
  )
  
  bt <- as.numeric(b$t)
  ci <- quantile(bt, c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    Model  = m,
    Cindex = as.numeric(b$t0),
    C_low  = as.numeric(ci[1]),
    C_high = as.numeric(ci[2]),
    ok_C   = sum(!is.na(bt))
  )
}))


## ----------------------------
## 12. Delta C-index versus SCORE2-only
## ----------------------------
R_d <- 200
set.seed(1)

deltaC_refit <- function(data, idx, vars_m, vars_ref) {
  d <- data[idx, , drop = FALSE]
  out <- tryCatch({
    fit_m   <- fit_cox_from_vars(d, vars_m)
    fit_ref <- fit_cox_from_vars(d, vars_ref)
    cindex_from_cox(fit_m) - cindex_from_cox(fit_ref)
  }, error = function(e) NA_real_)
  
  if (length(out) != 1) out <- NA_real_
  out
}

ref <- "cidx_5"

DeltaC_boot <- bind_rows(lapply(names(model_vars_surv), function(m) {
  if (m == ref) {
    return(data.frame(
      Model = m,
      DeltaC = 0,
      dC_low = 0,
      dC_high = 0,
      ok_dC = R_d
    ))
  }
  
  b <- boot(
    df_common,
    statistic = function(dat, idx) deltaC_refit(dat, idx, model_vars_surv[[m]], model_vars_surv[[ref]]),
    R = R_d
  )
  
  bt <- as.numeric(b$t)
  ci <- quantile(bt, c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    Model  = m,
    DeltaC = as.numeric(b$t0),
    dC_low = as.numeric(ci[1]),
    dC_high = as.numeric(ci[2]),
    ok_dC  = sum(!is.na(bt))
  )
}))


## ----------------------------
## 13. AUC at 10 years
## ----------------------------
get_auc_ci <- function(models, data, tau, B = 500, seed = 1) {
  set.seed(seed)
  
  sc <- Score(
    object   = models,
    formula  = Surv(time, event) ~ 1,
    data     = data,
    metrics  = "AUC",
    times    = tau,
    conf.int = TRUE,
    B        = B,
    summary  = "risk"
  )
  
  as.data.frame(sc$AUC$score) %>%
    dplyr::select(model, AUC, lower, upper) %>%
    rename(Model = model) %>%
    mutate(
      tau = tau,
      AUC_CI = sprintf("%.3f (%.3f–%.3f)", AUC, lower, upper)
    )
}

B_auc <- 500

auc10 <- get_auc_ci(models_cox, df_common, tau = 10, B = B_auc, seed = 1) %>%
  transmute(
    Model,
    AUC10 = AUC,
    AUC10_low = lower,
    AUC10_high = upper,
    AUC10_CI = AUC_CI
  )


## ----------------------------
## 14. Final summary table
## ----------------------------
labels <- tibble(
  Model = names(model_vars_surv),
  Label = c(
    "Family history (FH)",
    "PRS CHD",
    "FH + PRS",
    "PRS + SCORE2",
    "SCORE2",
    "Composite score",
    "Full: FH + PRS + SCORE2",
    "PRS + Age + Sex",
    "SCORE2_resid (Age/Sex removed)",
    "PRS + SCORE2_resid",
    "FH + PRS + SCORE2_resid"
  )
)

final_table <- labels %>%
  left_join(C_boot, by = "Model") %>%
  left_join(DeltaC_boot, by = "Model") %>%
  left_join(auc10, by = "Model") %>%
  mutate(
    Cindex_CI = sprintf("%.3f (%.3f–%.3f)", Cindex, C_low, C_high),
    DeltaC_CI = sprintf("%.3f (%.3f–%.3f)", DeltaC, dC_low, dC_high)
  ) %>%
  arrange(factor(Model, levels = names(model_vars_surv))) %>%
  dplyr::select(
    Model,
    Label,
    Cindex_CI,
    DeltaC_CI,
    AUC10_CI
  )

print(final_table)

final_table_paper <- as.data.frame(final_table)

fwrite(
  final_table_paper,
  file = outfile_eval,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)