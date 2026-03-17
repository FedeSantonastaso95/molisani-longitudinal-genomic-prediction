################################################################################
# PRS Logistic Regression Models (Age, Sex, PC1–PC5 adjusted)
#
# This script performs association analyses between standardized PRSs and
# disease outcomes under three scenarios:
#   1. All-event
#   2. Prevalent-only
#   3. Incident-only
#
# Models are adjusted for:
#   - Age
#   - Sex (when applicable)
#   - Genetic principal components (PC1–PC5)
#
# INPUT:
#   - Phenotype + covariate + PRS dataset
#
# OUTPUT:
#   - Table of ORs, 95% CI, and p-values
#
# USAGE (recommended):
#   Rscript 01_prs_association_models.R input_file output_file
#
################################################################################

# ============================ #
# 0. Libraries                 #
# ============================ #
library(data.table)
library(dplyr)
library(tibble)

# ============================ #
# 1. Input arguments           #
# ============================ #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 01_prs_association_models.R <input_file> <output_file>")
}

input_file  <- args[1]
output_file <- args[2]


# CI method
use_wald_ci <- TRUE


# ============================ #
# 2. Helper functions          #
# ============================ #

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
    if (all(candidate %in% names(data))) return(candidate)
  }
  
  stop("Could not detect PC1–PC5 columns.")
}

get_followup_column <- function(endpoint_code) {
  if (endpoint_code == "Chd0")    return("Chd0_Fnof_Fup2020")
  if (endpoint_code == "Stroke1") return("Stroke1_Fnof_Fup2020")
  if (endpoint_code == "Stroke3") return("Stroke3_Fnof_Fup2020")
  paste0("Fnf_", endpoint_code, "_Fup2020")
}

get_prs_column <- function(outcome_name) {
  paste0("PRED_eurW_", outcome_name, "_scale")
}

build_analysis_df <- function(data, endpoint_row, scenario) {
  
  fu_col <- get_followup_column(endpoint_row$column_name)
  bs_col <- endpoint_row$bs_col
  
  d <- data %>%
    mutate(
      .baseline = .data[[bs_col]],
      .followup = .data[[fu_col]]
    )
  
  if (scenario == "all_event") {
    
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


# ============================ #
# 3. Load data                 #
# ============================ #
df <- fread(input_file, data.table = FALSE) %>%
  filter(!FID %in% excluded_fids)


# ============================ #
# 4. Recode + PRS scaling      #
# ============================ #
analysis_df <- df %>%
  mutate(
    
    # Baseline recoding
    chd_baseline = case_when(
      Chd_Bs_New == 0 ~ 0,
      Chd_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    stroke_baseline = case_when(
      Cerebro_Bs_New == 0 ~ 0,
      Cerebro_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Cancer baseline definitions
    breast_baseline = ifelse(Mt_Diagnosis_New == 1 &
                               (Mt1_Icd9 == 174 | Mt2_Icd9_New == 174 | Mt3_Icd9_New == 174), 1, 0),
    
    prostate_baseline = ifelse(Mt_Diagnosis_New == 1 &
                                 (Mt1_Icd9 == 185 | Mt2_Icd9_New == 185 | Mt3_Icd9_New == 185), 1, 0),
    
    colorectal_baseline = ifelse(Mt_Diagnosis_New == 1 &
                                   (Mt1_Icd9 %in% c(153,154) |
                                      Mt2_Icd9_New %in% c(153,154) |
                                      Mt3_Icd9_New %in% c(153,154)), 1, 0),
    
    lung_baseline = ifelse(Mt_Diagnosis_New == 1 &
                             (Mt1_Icd9 == 162 | Mt2_Icd9_New == 162 | Mt3_Icd9_New == 162), 1, 0),
    
    renal_baseline = ifelse(Mt_Diagnosis_New == 1 &
                              (Mt1_Icd9 == 189 | Mt2_Icd9_New == 189 | Mt3_Icd9_New == 189), 1, 0),
    
    pancreatic_baseline = ifelse(Mt_Diagnosis_New == 1 &
                                   (Mt1_Icd9 == 157 | Mt2_Icd9_New == 157 | Mt3_Icd9_New == 157), 1, 0),
    
    # PRS scaling
    across(starts_with("PRED_eurW_"), ~ as.numeric(scale(.)), .names = "{.col}_scale")
  ) %>%
  mutate(
    # Remove secondary events
    Chd0_Fnof_Fup2020 = if_else(
      chd_baseline == 1 & Chd0_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Chd0_Fnof_Fup2020)
    ),
    Stroke1_Fnof_Fup2020 = if_else(
      stroke_baseline == 1 & Stroke1_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Stroke1_Fnof_Fup2020)
    )
  )


# ============================ #
# 5. Detect PCs                #
# ============================ #
pc_cols <- detect_pc_cols(analysis_df)


# ============================ #
# 6. Endpoint definitions      #
# ============================ #
endpoints <- data.frame(
  column_name = c("Bc","Pc","Crc","Chd0","Lc","Rnc","Pnc","Stroke1","Stroke3"),
  sex_stratified = c("yes","yes","no","no","no","no","no","no","no"),
  outcome = c("breast","prostate","colorectal","CHD","lung","renal","pancreatic","stroke","stroke"),
  acro = c("BC","PC","CRC","CHD","LC","RNC","PNC","Stroke","I-Stroke"),
  bs_col = c("breast_baseline","prostate_baseline","colorectal_baseline","chd_baseline",
             "lung_baseline","renal_baseline","pancreatic_baseline","stroke_baseline","stroke_baseline"),
  drop_secondary_events = c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE),
  stringsAsFactors = FALSE
)


scenario_map <- data.frame(
  scenario = c("all_event","prevalent","incident"),
  stringsAsFactors = FALSE
)


# ============================ #
# 7. Run models                #
# ============================ #
final_results <- data.frame()

for (i in seq_len(nrow(endpoints))) {
  
  ep <- endpoints[i, ]
  prs_col <- get_prs_column(ep$outcome)
  
  for (sc in scenario_map$scenario) {
    
    df2 <- build_analysis_df(analysis_df, ep, sc)
    
    if (ep$sex_stratified == "yes") {
      needed <- c(".y", prs_col, "Age", pc_cols)
      rhs <- paste(c(prs_col, "Age", pc_cols), collapse = " + ")
    } else {
      needed <- c(".y", prs_col, "Age", "Sex", pc_cols)
      rhs <- paste(c(prs_col, "Age", "factor(Sex)", pc_cols), collapse = " + ")
    }
    
    df_model <- df2 %>%
      select(all_of(needed)) %>%
      filter(if_all(everything(), ~ !is.na(.)))
    
    if (sum(df_model$.y == 1) == 0 || sum(df_model$.y == 0) == 0) next
    
    fit <- glm(as.formula(paste0(".y ~ ", rhs)),
               data = df_model,
               family = binomial())
    
    co <- summary(fit)$coefficients
    ci <- if (use_wald_ci) confint.default(fit) else confint(fit)
    
    final_results <- bind_rows(
      final_results,
      data.frame(
        endpoint = ep$acro,
        scenario = sc,
        OR = round(exp(co[prs_col, "Estimate"]), 2),
        CI_lower = round(exp(ci[prs_col, 1]), 2),
        CI_upper = round(exp(ci[prs_col, 2]), 2),
        pvalue = format(co[prs_col, "Pr(>|z|)"], scientific = TRUE, digits = 2)
      )
    )
  }
}


# ============================ #
# 8. Save results              #
# ============================ #
fwrite(final_results, output_file, sep = "\t", na = "NA", quote = FALSE)