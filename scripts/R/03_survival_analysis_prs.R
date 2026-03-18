############################################################
# Survival models and table export
#
# This section:
#   1. prepares endpoint-specific survival datasets
#   2. fits Cox proportional hazards models using:
#        - continuous standardized PRS
#        - PRS quintiles (Q1 as reference)
#   3. stores hazard ratios and confidence intervals
#   4. exports results tables for:
#        - continuous PRS models
#        - quintile-based PRS models
#        - case counts after baseline exclusions
#
# Notes:
# - Baseline-prevalent cases are removed from each endpoint-specific risk set.
# - For CHD and stroke, secondary events are already handled upstream by
#   recoding follow-up event flags to NA when prevalent disease is present
#   at baseline.
# - Models are adjusted for:
#     * Age
#     * Sex, when appropriate
#     * PC1-PC5, when available
############################################################


## ----------------------------
## Prepare survival dataset for a single endpoint
## ----------------------------
prep_surv_df <- function(df,
                         pattern,
                         prs_outcome,
                         sex_filter = "all",
                         bs_col = NULL,
                         censor_date = as.Date("2020-12-31"),
                         cap_years = 16) {
  
  cd <- censor_date
  cy <- cap_years
  
  # Identify endpoint-specific columns from the column name pattern
  cols <- grep(pattern, names(df), value = TRUE)
  
  # Follow-up date column
  col_date <- cols[grepl("(Date|Dat)", cols, ignore.case = TRUE)][1]
  
  # Follow-up event flag column
  col_fnf <- cols[
    grepl("Fno?f", cols, ignore.case = TRUE) &
      !grepl("(Date|Dat)", cols, ignore.case = TRUE)
  ][1]
  
  # Explicit fallback for CHD naming if automatic matching fails
  if ((is.na(col_date) || is.na(col_fnf)) && pattern == "Chd0") {
    if ("Chd0_Fnof_Date_Fup2020" %in% names(df)) col_date <- "Chd0_Fnof_Date_Fup2020"
    if ("Chd0_Fnof_Fup2020" %in% names(df))      col_fnf  <- "Chd0_Fnof_Fup2020"
  }
  
  # Safety checks
  if (is.na(col_date) || is.na(col_fnf)) {
    stop(
      paste(
        "Cannot identify date/event columns for pattern:", pattern,
        "\nCandidate columns found:", paste(cols, collapse = ", ")
      )
    )
  }
  
  if (!is.null(bs_col) && !(bs_col %in% names(df))) {
    stop(paste("Baseline prevalence column not found:", bs_col))
  }
  
  # Build time-to-event variables
  #
  # status = 1 only if:
  #   - a valid exit date is available
  #   - exit date is different from the administrative censoring date
  #   - the endpoint-specific event flag is equal to 1
  #
  # If an event flag has already been recoded to NA upstream
  # (e.g. secondary CHD/stroke events), it will not contribute as an event.
  df1 <- df %>%
    mutate(
      Recruitment_date = as.Date(parse_date_safe(Recruitment_date)),
      exit_date_raw    = as.character(.data[[col_date]]),
      exit_date        = as.Date(parse_date_safe(exit_date_raw))
    ) %>%
    mutate(
      status   = ifelse(!is.na(exit_date) & exit_date != cd & .data[[col_fnf]] == 1, 1L, 0L),
      exit_eff = dplyr::coalesce(exit_date, cd),
      time     = as.numeric(difftime(exit_eff, Recruitment_date, units = "days")) / 365.25,
      time     = pmin(pmax(time, 0), cy)
    )
  
  # Apply sex restriction for sex-specific endpoints
  if ("Sex" %in% names(df1)) df1 <- df1 %>% mutate(Sex = factor(Sex))
  if (sex_filter == "female") df1 <- df1 %>% filter(Sex == 0)
  if (sex_filter == "male")   df1 <- df1 %>% filter(Sex == 1)
  
  # Remove baseline-prevalent cases to define an incident-only risk set
  if (!is.null(bs_col)) {
    df1 <- df1 %>% filter(is.na(.data[[bs_col]]) | .data[[bs_col]] == 0)
  }
  
  list(
    df = df1,
    prs_var = paste0("PRED_eurW_", prs_outcome)
  )
}


## ----------------------------
## Fit Cox models and extract results
## - continuous PRS model
## - quintile-based PRS model
## ----------------------------
fit_survival_models <- function(df,
                                prs_var,
                                outcome_label,
                                sex_filter = "all") {
  
  # Assign PRS quintiles
  df_q <- make_quintiles(df, prs_var)
  
  # Build covariate set:
  # - Age always
  # - PC1-PC5 if available
  # - Sex only when appropriate
  covs <- build_covariates(df_q, sex_filter = sex_filter)
  
  # Make sure PCs are numeric
  pc_terms <- covs[grepl("^PC[0-9]+$", covs)]
  if (length(pc_terms) > 0) {
    df_q <- df_q %>%
      mutate(across(all_of(pc_terms), ~ as.numeric(.)))
  }
  
  # Make sure Sex is treated as factor if included in the model
  if ("Sex" %in% covs) {
    df_q <- df_q %>% mutate(Sex = factor(Sex))
  }
  
  # Keep complete cases only for the variables used in the model
  needed_vars <- unique(c("time", "status", prs_var, "PRS_q", covs))
  df_model <- df_q %>%
    select(all_of(needed_vars)) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  # ==========================================================
  # 1. Continuous PRS model
  # ==========================================================
  #
  # PRS is standardized inside the model, so the HR corresponds to
  # a 1-SD increase in the PRS.
  rhs_cont <- paste(c(paste0("scale(", prs_var, ")"), covs), collapse = " + ")
  f_cont   <- as.formula(paste("Surv(time, status) ~", rhs_cont))
  
  fit_cont <- coxph(f_cont, data = df_model, ties = "efron")
  sm_cont  <- summary(fit_cont)
  
  ci_cont   <- sm_cont$conf.int
  coef_cont <- sm_cont$coefficients
  row_prs   <- grep("scale", rownames(ci_cont))
  
  cont_stats <- data.frame(
    outcome = outcome_label,
    prs_type = "continuous",
    term = rownames(ci_cont)[row_prs],
    HR = ci_cont[row_prs, "exp(coef)"],
    CI_low = ci_cont[row_prs, "lower .95"],
    CI_high = ci_cont[row_prs, "upper .95"],
    p_value = coef_cont[row_prs, "Pr(>|z|)"],
    N_total = nrow(df_model),
    N_cases = sum(df_model$status == 1, na.rm = TRUE),
    N_censored = sum(df_model$status == 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # ==========================================================
  # 2. Quintile-based PRS model
  # ==========================================================
  #
  # Q1 is used as the reference category.
  df_model$PRS_q <- relevel(df_model$PRS_q, ref = "Q1")
  
  rhs_cat <- paste(c("PRS_q", covs), collapse = " + ")
  f_cat   <- as.formula(paste("Surv(time, status) ~", rhs_cat))
  
  fit_cat <- coxph(f_cat, data = df_model, ties = "efron")
  sm_cat  <- summary(fit_cat)
  
  ci_cat   <- sm_cat$conf.int
  coef_cat <- sm_cat$coefficients
  rows_q   <- grep("^PRS_q", rownames(ci_cat))
  
  cat_stats <- data.frame(
    outcome = outcome_label,
    prs_type = "quintile",
    term = c("PRS_qQ1", rownames(ci_cat)[rows_q]),
    HR = c(1, ci_cat[rows_q, "exp(coef)"]),
    CI_low = c(1, ci_cat[rows_q, "lower .95"]),
    CI_high = c(1, ci_cat[rows_q, "upper .95"]),
    p_value = c(NA, coef_cat[rows_q, "Pr(>|z|)"]),
    N_total = nrow(df_model),
    N_cases = sum(df_model$status == 1, na.rm = TRUE),
    N_censored = sum(df_model$status == 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Count participants and events within each quintile
  quintile_counts <- df_model %>%
    group_by(PRS_q) %>%
    summarise(
      N_total = n(),
      N_cases = sum(status == 1, na.rm = TRUE),
      N_censored = sum(status == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(outcome = outcome_label) %>%
    select(outcome, PRS_q, N_total, N_cases, N_censored)
  
  list(
    cont_stats = cont_stats,
    cat_stats = cat_stats,
    quintile_counts = quintile_counts
  )
}


## ----------------------------
## Endpoint specification
## ----------------------------
endpoints <- data.frame(
  column_name    = c("Bc", "Pc", "Crc", "Chd0", "Lc", "Rnc", "Pnc", "Stroke1", "Stroke3"),
  sex_stratified = c("yes", "yes", "no", "no", "no", "no", "no", "no", "no"),
  title          = c(
    "Breast cancer (BC)",
    "Prostate cancer (PC)",
    "Colorectal cancer (CRC)",
    "Coronary Heart Disease (CHD)",
    "Lung cancer (LC)",
    "Renal Cell Carcinoma (RCC)",
    "Pancreatic cancer (PCa)",
    "Cerebrovascular accident (CVA)",
    "Ischemic Cerebrovascular accident (I-CVA)"
  ),
  outcome        = c("breast", "prostate", "colorectal", "CHD", "lung", "renal", "pancreatic", "stroke", "stroke"),
  acro           = c("BC", "PC", "CRC", "CHD", "LC", "RCC", "PCa", "CVA", "I-CVA"),
  bs_col         = c(
    "breast_Bs_final",
    "prostate_Bs_final",
    "colorectal_Bs_final",
    "Chd_Bs_New_final",
    "lung_Bs_final",
    "renal_Bs_final",
    "pancreatic_Bs_final",
    "Cerebro_Bs_New_final",
    "Cerebro_Bs_New_final"
  ),
  stringsAsFactors = FALSE
)


## ----------------------------
## Run survival models across endpoints
## ----------------------------
case_counts <- list()
results_cont_list <- vector("list", nrow(endpoints))
results_cat_list  <- vector("list", nrow(endpoints))
quintile_count_list <- vector("list", nrow(endpoints))

for (i in seq_len(nrow(endpoints))) {
  
  # Endpoint-specific metadata
  pat  <- endpoints$column_name[i]
  outc <- endpoints$outcome[i]
  ttl  <- endpoints$title[i]
  acr  <- endpoints$acro[i]
  bsc  <- endpoints$bs_col[i]
  
  # Restrict to the appropriate sex for breast and prostate cancer
  sex_f <- if (endpoints$sex_stratified[i] == "yes") {
    if (outc == "breast") {
      "female"
    } else if (outc == "prostate") {
      "male"
    } else {
      "all"
    }
  } else {
    "all"
  }
  
  # Prepare endpoint-specific survival dataset
  prep <- prep_surv_df(
    df = raw,
    pattern = pat,
    prs_outcome = outc,
    sex_filter = sex_f,
    bs_col = bsc
  )
  
  # Case counts after all endpoint-specific exclusions
  n_total <- nrow(prep$df)
  n_cases <- sum(prep$df$status == 1, na.rm = TRUE)
  n_cens  <- sum(prep$df$status == 0, na.rm = TRUE)
  
  case_counts[[ttl]] <- data.frame(
    outcome = ttl,
    N_total = n_total,
    N_cases = n_cases,
    N_censored = n_cens,
    stringsAsFactors = FALSE
  )
  
  # Fit survival models and extract tables
  fit_out <- fit_survival_models(
    df = prep$df,
    prs_var = prep$prs_var,
    outcome_label = acr,
    sex_filter = sex_f
  )
  
  results_cont_list[[i]]   <- fit_out$cont_stats
  results_cat_list[[i]]    <- fit_out$cat_stats
  quintile_count_list[[i]] <- fit_out$quintile_counts
}


## ----------------------------
## Combine all output tables
## ----------------------------
HR_continuous <- do.call(rbind, results_cont_list)
HR_quintiles  <- do.call(rbind, results_cat_list)
cases_table   <- do.call(rbind, case_counts)
quintile_case_counts <- do.call(rbind, quintile_count_list)


## ----------------------------
## Print case counts
## ----------------------------
print(cases_table)
print(quintile_case_counts)


## ----------------------------
## Create output directory
## ----------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


## ----------------------------
## Save output tables
## ----------------------------

# Continuous PRS Cox models
fwrite(
  HR_continuous,
  file.path(outdir, "HR_continuous_PRS_allOutcomes.txt"),
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# Quintile-based PRS Cox models
fwrite(
  HR_quintiles,
  file.path(outdir, "HR_quintiles_PRS_allOutcomes.txt"),
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# Case counts after removal of baseline-prevalent individuals
fwrite(
  cases_table,
  file.path(outdir, "case_counts_afterRemovingPrevalent.txt"),
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# Event counts within each PRS quintile
fwrite(
  quintile_case_counts,
  file.path(outdir, "case_counts_by_quintile_allOutcomes.txt"),
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)