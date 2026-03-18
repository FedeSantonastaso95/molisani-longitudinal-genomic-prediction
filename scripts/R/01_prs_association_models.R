################################################################################
# 01_prs_association_models.R
#
# Description:
# This script evaluates the association between standardised polygenic risk
# scores (PRSs) and disease outcomes using logistic regression models.
#
# Three outcome definitions are considered for each endpoint:
#   1. All-event       -> prevalent at baseline OR event during follow-up
#   2. Prevalent-only  -> prevalent disease at baseline only
#   3. Incident-only   -> follow-up events among baseline-free participants only
#
# Models are adjusted for:
#   - PRS (standardised to z-score)
#   - Age at baseline
#   - Sex (for non-sex-specific outcomes only)
#   - Genetic principal components PC1-PC5
#
# Notes:
# - For baseline CHD and stroke variables:
#     0 = No
#     1 = Yes, documented
#     2 = Yes, not documented
#   In the analysis, both 1 and 2 are treated as prevalent disease.
#
# - For CHD and stroke, follow-up events occurring in participants already
#   prevalent at baseline are treated as secondary events and are excluded
#   from incident-style handling.
#
# Usage:
#   Rscript 01_prs_association_models.R <input_file> <output_file>
#
################################################################################


# ============================ #
# 0. Load required libraries   #
# ============================ #
library(data.table)
library(dplyr)
library(tibble)


# ============================ #
# 1. Parse command-line input  #
# ============================ #

args <- commandArgs(trailingOnly = TRUE)
script_name <- basename(commandArgs()[1])

if (length(args) < 2) {
  stop(
    paste("Usage: Rscript", script_name, "<input_file> <output_file>")
  )
}

input_file  <- args[1]
output_file <- args[2]


# ============================ #
# 2. Analysis options          #
# ============================ #

# Confidence interval method:
# TRUE  = Wald CI via confint.default() (faster)
# FALSE = profile-likelihood CI via confint() (slower)
use_wald_ci <- TRUE


# ============================ #
# 3. Helper functions          #
# ============================ #

# Detect PC1-PC5 columns using common naming conventions
detect_pc_cols <- function(data, n_pc = 5) {
  
  candidate_sets <- list(
    paste0("PC", 1:n_pc),
    paste0("PC_", 1:n_pc),
    paste0("pc", 1:n_pc),
    paste0("pc_", 1:n_pc),
    paste0("EV", 1:n_pc),
    paste0("EV_", 1:n_pc)
  )
  
  for (candidate in candidate_sets) {
    if (all(candidate %in% names(data))) {
      return(candidate)
    }
  }
  
  stop(
    "Could not detect 5 principal component columns.\n",
    "Searched for: PC1-PC5, PC_1-PC_5, pc1-pc5, pc_1-pc_5, EV1-EV5, EV_1-EV_5.\n",
    "Available columns containing 'PC' or 'EV': ",
    paste(grep("PC|pc|EV|ev", names(data), value = TRUE), collapse = ", ")
  )
}


# Return the follow-up event column corresponding to each endpoint
get_followup_column <- function(endpoint_code) {
  if (endpoint_code == "Chd0")    return("Chd0_Fnof_Fup2020")
  if (endpoint_code == "Stroke1") return("Stroke1_Fnof_Fup2020")
  if (endpoint_code == "Stroke3") return("Stroke3_Fnof_Fup2020")
  return(paste0("Fnf_", endpoint_code, "_Fup2020"))
}


# Return the scaled PRS column corresponding to the endpoint outcome
get_prs_column <- function(outcome_name) {
  paste0("PRED_eurW_", outcome_name, "_scale")
}


# Build the analysis dataset according to the selected scenario
build_analysis_df <- function(data, endpoint_row, scenario) {
  
  followup_col <- get_followup_column(endpoint_row$column_name)
  baseline_col <- endpoint_row$bs_col
  
  d <- data %>%
    mutate(
      .baseline = .data[[baseline_col]],
      .followup = .data[[followup_col]]
    )
  
  if (scenario == "all_event") {
    
    # For selected endpoints (e.g. CHD and stroke), remove participants
    # with both prevalent disease at baseline and a follow-up event flag
    if (isTRUE(endpoint_row$drop_secondary_events)) {
      d <- d %>% filter(!(.baseline == 1 & .followup == 1))
    }
    
    d <- d %>%
      mutate(.y = if_else(.baseline == 1 | .followup == 1, 1L, 0L))
  }
  
  if (scenario == "prevalent") {
    d <- d %>%
      mutate(.y = if_else(.baseline == 1, 1L, 0L))
  }
  
  if (scenario == "incident") {
    d <- d %>%
      filter(.baseline == 0) %>%
      mutate(.y = if_else(.followup == 1, 1L, 0L))
  }
  
  return(d)
}


# Build the logistic regression formula
build_model_formula <- function(prs_col, sex_stratified, pc_cols) {
  
  # Model covariates:
  # - PRS: standardised polygenic risk score (z-score)
  # - Age: age at baseline, modeled as a continuous variable
  # - Sex: biological sex, included only for non-sex-specific outcomes
  # - PC1-PC5: genetic principal components used to adjust for population
  #            structure / ancestry differences
  
  if (sex_stratified) {
    rhs <- paste(c(
      prs_col,
      "Age",
      pc_cols
    ), collapse = " + ")
  } else {
    rhs <- paste(c(
      prs_col,
      "Age",
      "factor(Sex)",
      pc_cols
    ), collapse = " + ")
  }
  
  as.formula(paste0(".y ~ ", rhs))
}


# ============================ #
# 4. Load dataset              #
# ============================ #

df <- fread(input_file, data.table = FALSE) %>%
  filter(!FID %in% excluded_fids)


# ============================ #
# 5. Recode baseline variables #
#    and standardize PRSs      #
# ============================ #

analysis_df <- df %>%
  mutate(
    
    # ----------------------------------------------------------------------
    # Baseline CHD recoding
    #
    # Original coding:
    #   0 = No
    #   1 = Yes, documented
    #   2 = Yes, not documented
    #
    # For analysis purposes:
    #   - 0 is coded as no baseline disease
    #   - 1 and 2 are both coded as prevalent disease
    #   - any other value is set to NA
    # ----------------------------------------------------------------------
    chd_baseline = case_when(
      Chd_Bs_New == 0 ~ 0,
      Chd_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # ----------------------------------------------------------------------
    # Baseline cerebrovascular disease (stroke) recoding
    #
    # Original coding:
    #   0 = No
    #   1 = Yes, documented
    #   2 = Yes, not documented
    #
    # For analysis purposes:
    #   - 0 is coded as no baseline disease
    #   - 1 and 2 are both coded as prevalent disease
    #   - any other value is set to NA
    # ----------------------------------------------------------------------
    stroke_baseline = case_when(
      Cerebro_Bs_New == 0 ~ 0,
      Cerebro_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # ----------------------------------------------------------------------
    # Baseline cancer definitions based on ICD-9 codes
    #
    # A participant is classified as prevalent at baseline if:
    #   - Mt_Diagnosis_New == 1, and
    #   - at least one ICD-9 diagnosis field matches the disease code
    # ----------------------------------------------------------------------
    breast_baseline = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 174 | Mt2_Icd9_New == 174 | Mt3_Icd9_New == 174),
      1, 0
    ),
    
    prostate_baseline = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 185 | Mt2_Icd9_New == 185 | Mt3_Icd9_New == 185),
      1, 0
    ),
    
    colorectal_baseline = ifelse(
      Mt_Diagnosis_New == 1 &
        (
          Mt1_Icd9 %in% c(153, 154) |
            Mt2_Icd9_New %in% c(153, 154) |
            Mt3_Icd9_New %in% c(153, 154)
        ),
      1, 0
    ),
    
    lung_baseline = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 162 | Mt2_Icd9_New == 162 | Mt3_Icd9_New == 162),
      1, 0
    ),
    
    renal_baseline = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 189 | Mt2_Icd9_New == 189 | Mt3_Icd9_New == 189),
      1, 0
    ),
    
    pancreatic_baseline = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 157 | Mt2_Icd9_New == 157 | Mt3_Icd9_New == 157),
      1, 0
    ),
    
    # ----------------------------------------------------------------------
    # Standardize PRSs to z-scores
    #
    # Odds ratios from the regression models can therefore be interpreted as
    # the change in odds associated with a 1-standard-deviation increase in PRS.
    # ----------------------------------------------------------------------
    PRED_eurW_CHD_scale        = as.numeric(scale(PRED_eurW_CHD)),
    PRED_eurW_stroke_scale     = as.numeric(scale(PRED_eurW_stroke)),
    PRED_eurW_breast_scale     = as.numeric(scale(PRED_eurW_breast)),
    PRED_eurW_prostate_scale   = as.numeric(scale(PRED_eurW_prostate)),
    PRED_eurW_colorectal_scale = as.numeric(scale(PRED_eurW_colorectal)),
    PRED_eurW_lung_scale       = as.numeric(scale(PRED_eurW_lung)),
    PRED_eurW_renal_scale      = as.numeric(scale(PRED_eurW_renal)),
    PRED_eurW_pancreatic_scale = as.numeric(scale(PRED_eurW_pancreatic))
  ) %>%
  mutate(
    
    # ----------------------------------------------------------------------
    # Secondary CHD event handling
    #
    # If a participant is already prevalent for CHD at baseline and also has
    # a CHD follow-up event flag, the follow-up flag is set to NA so that
    # this individual does not contribute as an incident case.
    # ----------------------------------------------------------------------
    Chd0_Fnof_Fup2020 = if_else(
      chd_baseline == 1 & Chd0_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Chd0_Fnof_Fup2020)
    ),
    
    # ----------------------------------------------------------------------
    # Secondary stroke event handling
    #
    # If a participant is already prevalent for stroke at baseline and also
    # has a stroke follow-up event flag, the follow-up flag is set to NA so
    # that this individual does not contribute as an incident case.
    # ----------------------------------------------------------------------
    Stroke1_Fnof_Fup2020 = if_else(
      stroke_baseline == 1 & Stroke1_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Stroke1_Fnof_Fup2020)
    )
  )


# ============================ #
# 6. Detect PC1-PC5 columns    #
# ============================ #

pc_cols <- detect_pc_cols(analysis_df, n_pc = 5)


# ============================ #
# 7. Define endpoints          #
# ============================ #

endpoints <- data.frame(
  column_name = c("Bc", "Pc", "Crc", "Chd0", "Lc", "Rnc", "Pnc", "Stroke1", "Stroke3"),
  sex_stratified = c("yes", "yes", "no", "no", "no", "no", "no", "no", "no"),
  title = c(
    "Breast cancer (BC)",
    "Prostate cancer (PC)",
    "Colorectal cancer (CRC)",
    "Coronary Heart Disease (CHD)",
    "Lung cancer (LC)",
    "Renal Cell Carcinoma (RCC)",
    "Pancreatic cancer (PCa)",
    "Stroke",
    "Ischemic stroke (I-Stroke)"
  ),
  outcome = c("breast", "prostate", "colorectal", "CHD", "lung", "renal", "pancreatic", "stroke", "stroke"),
  acro = c("BC", "PC", "CRC", "CHD", "LC", "RNC", "PNC", "Stroke", "I-Stroke"),
  bs_col = c(
    "breast_baseline",
    "prostate_baseline",
    "colorectal_baseline",
    "chd_baseline",
    "lung_baseline",
    "renal_baseline",
    "pancreatic_baseline",
    "stroke_baseline",
    "stroke_baseline"
  ),
  drop_secondary_events = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

scenario_map <- data.frame(
  scenario = c("all_event", "prevalent", "incident"),
  scenario_label = c("All-event", "Prevalent-only", "Incident-only"),
  stringsAsFactors = FALSE
)


# ============================ #
# 8. Run association models    #
# ============================ #

final_results <- data.frame()

for (i in seq_len(nrow(endpoints))) {
  
  endpoint_row <- endpoints[i, ]
  prs_col <- get_prs_column(endpoint_row$outcome)
  
  for (j in seq_len(nrow(scenario_map))) {
    
    scenario_name  <- scenario_map$scenario[j]
    scenario_label <- scenario_map$scenario_label[j]
    
    # Build analysis dataset for the selected endpoint and scenario
    df_scenario <- build_analysis_df(analysis_df, endpoint_row, scenario_name)
    
    # Build regression formula
    model_formula <- build_model_formula(
      prs_col = prs_col,
      sex_stratified = endpoint_row$sex_stratified == "yes",
      pc_cols = pc_cols
    )
    
    # Select variables required by the model
    if (endpoint_row$sex_stratified == "yes") {
      required_vars <- c(".y", prs_col, "Age", pc_cols)
    } else {
      required_vars <- c(".y", prs_col, "Age", "Sex", pc_cols)
    }
    
    # Restrict to complete cases
    df_model <- df_scenario %>%
      select(all_of(required_vars)) %>%
      filter(if_all(everything(), ~ !is.na(.)))
    
    # Count cases and controls
    n_cases <- sum(df_model$.y == 1)
    n_controls <- sum(df_model$.y == 0)
    
    # Skip models with no variation in outcome
    if (n_cases == 0 || n_controls == 0) {
      final_results <- bind_rows(
        final_results,
        data.frame(
          endpoint = endpoint_row$acro,
          column_name = endpoint_row$column_name,
          title = endpoint_row$title,
          outcome = endpoint_row$outcome,
          scenario = scenario_name,
          scenario_label = scenario_label,
          n_cases = n_cases,
          n_controls = n_controls,
          variable = prs_col,
          beta = NA_real_,
          SE = NA_real_,
          pvalue = NA_character_,
          OR = NA_real_,
          CI_lower = NA_real_,
          CI_upper = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    
    # Fit logistic regression model
    fit <- glm(
      formula = model_formula,
      data = df_model,
      family = binomial(link = "logit")
    )
    
    fit_summary <- summary(fit)
    
    # Extract only the PRS coefficient
    prs_result <- as.data.frame(fit_summary$coefficients) %>%
      rownames_to_column("variable") %>%
      filter(variable == prs_col)
    
    # Compute confidence intervals
    ci_matrix <- if (use_wald_ci) {
      suppressMessages(confint.default(fit))
    } else {
      suppressMessages(confint(fit))
    }
    
    ci_prs <- ci_matrix[prs_col, , drop = FALSE]
    
    # Store formatted results
    result_row <- prs_result %>%
      transmute(
        endpoint = endpoint_row$acro,
        column_name = endpoint_row$column_name,
        title = endpoint_row$title,
        outcome = endpoint_row$outcome,
        scenario = scenario_name,
        scenario_label = scenario_label,
        n_cases = n_cases,
        n_controls = n_controls,
        variable = variable,
        beta = Estimate,
        SE = `Std. Error`,
        pvalue = `Pr(>|z|)`,
        OR = exp(beta),
        CI_lower = exp(ci_prs[1]),
        CI_upper = exp(ci_prs[2])
      ) %>%
      mutate(
        across(c(beta, SE, OR, CI_lower, CI_upper), ~ round(., 2)),
        pvalue = format(pvalue, digits = 2, scientific = TRUE)
      )
    
    final_results <- bind_rows(final_results, result_row)
  }
}


# ============================ #
# 9. Export results            #
# ============================ #

fwrite(
  final_results,
  file = output_file,
  row.names = FALSE,
  na = "NA",
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)