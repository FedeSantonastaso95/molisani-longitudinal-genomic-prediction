############################################################
# CHD survival analysis using PRS and polygenic scores
#
# Description:
# This script performs Cox proportional hazards analyses for
# incident Coronary Heart Disease (CHD).
#
# Predictors tested:
#   - the CHD-specific PRS
#   - polygenic scores (PGS) for intermediate endophenotypes
#     / risk factors, as defined in the external legend file
#
# Analysis design:
#   - incident-only risk set
#   - baseline-prevalent CHD cases are excluded
#   - secondary CHD events are handled upstream by setting the
#     follow-up CHD event flag to NA when prevalent CHD is already
#     present at baseline (thus preventing their inclusion as incident events)
#
# Model adjustment:
#   - Age
#   - Sex
#   - PC1–PC5
#
# Event coding note:
# The Cox proportional hazards model (survival::coxph) is invariant
# to the specific numeric coding of the event indicator, as long as
# the ordering is preserved (i.e. higher value = event).
#
# Therefore, both of the following encodings are valid:
#   - 0 = no event (censored), 1 = event
#   - 1 = no event (censored), 2 = event
#
# In this script:
#   - the event indicator is explicitly recoded to a binary variable:
#       status = 1 (event), 0 (no event)
#   - this ensures consistency across all endpoints and models
#
# Output:
#   A results table with HR, 95% CI, p-value, case counts, and
#   model details for each PRS/PGS tested.
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
## Input files
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

legend_file <- "path/to/2025.01.30_SurvivalCurves_riskFactors_diseases.txt"

outfile <- "path/to/output/CHD_cox_incident_vs_nonincident_AgeSexPC1to5.txt"


## ----------------------------
## 1. Read data and define baseline flags
## ----------------------------
df <- fread(infile, data.table = FALSE) %>%
  mutate(
    # Baseline coding for CHD:
    #   0 = No
    #   1 = Yes, documented
    #   2 = Yes, not documented
    #
    # In the analysis, both 1 and 2 are treated as prevalent CHD.
    Chd_Bs_New_final = case_when(
      Chd_Bs_New == 0 ~ 0,
      Chd_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    # Secondary CHD events are not counted as incident events:
    # if CHD is already prevalent at baseline and the follow-up event flag is 1,
    # the follow-up flag is set to NA.
    Chd0_Fnof_Fup2020 = if_else(
      Chd_Bs_New_final == 1 & Chd0_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Chd0_Fnof_Fup2020)
    )
  )


## ----------------------------
## 2. Resolve CHD event/date columns
## ----------------------------
resolve_outcome_columns <- function(dat, column_name) {
  
  cols <- grep(column_name, names(dat), value = TRUE)
  
  # Prefer date columns containing "Dat" or "Date"
  col_dat <- cols[grepl("Dat", cols)]
  if (length(col_dat) > 1 && any(grepl("Date", col_dat))) {
    col_dat <- col_dat[grepl("Date", col_dat)]
  }
  
  # Prefer event flag columns containing "Fnof"
  col_fnf <- cols[grepl("Fnof", cols)]
  col_fnf <- col_fnf[!grepl("Date|Annipers|Dat", col_fnf)]
  
  # Fallback to broader "Fn" search if needed
  if (length(col_fnf) == 0) {
    col_fnf <- cols[grepl("Fn", cols)]
    col_fnf <- col_fnf[!grepl("Dat|Date|Annipers", col_fnf)]
  }
  
  if (length(col_dat) != 1 || length(col_fnf) != 1) {
    stop(
      paste0(
        "Ambiguous CHD columns.",
        " Dat: ", paste(col_dat, collapse = ", "),
        " | Fnf/Fnof: ", paste(col_fnf, collapse = ", ")
      )
    )
  }
  
  list(
    column_dataexit = col_dat,
    column_fnf = col_fnf
  )
}


## ----------------------------
## 3. Restrict to the incident-only CHD risk set
## ----------------------------
#
# Baseline-prevalent CHD cases are excluded so that the comparison is
# incident CHD vs no incident CHD among those free of CHD at baseline.
filter_prevalent_only <- function(dat, bs_final) {
  
  if (!(bs_final %in% names(dat))) {
    stop(paste0("Baseline prevalence column not found: ", bs_final))
  }
  
  dat %>% filter(.data[[bs_final]] == 0L)
}


## ----------------------------
## 4. Fit CHD Cox models for PRS and PGS predictors
## ----------------------------
cox_table_function_chd <- function(dat,
                                   legend_path,
                                   bs_col_final = "Chd_Bs_New_final") {
  
  # Identify CHD follow-up date and event columns
  cols_res <- resolve_outcome_columns(dat, "Chd0")
  column_dataexit <- cols_res$column_dataexit
  column_fnf      <- cols_res$column_fnf
  
  # Build binary event status:
  # - 1 = incident CHD
  # - 0 = no incident CHD
  #
  # Secondary CHD events already recoded to NA upstream will not contribute
  # as incident events here.
  df1 <- dat %>%
    mutate(
      status = ifelse(.data[[column_fnf]] == 1, 1L, 0L)
    )
  
  # Build follow-up time from recruitment to CHD event/censoring
  df1$Recruitment_date <- as.Date(
    gsub("/", "-", df1$Recruitment_date),
    "%d-%m-%Y"
  )
  
  df2 <- df1 %>%
    mutate(
      Dataexit = as.Date(as.character(.data[[column_dataexit]])),
      Dataexit = ifelse(is.na(Dataexit), as.Date("2020-12-31"), Dataexit),
      Dataexit = as.Date(Dataexit, origin = "1970-01-01"),
      time = as.numeric(difftime(Dataexit, Recruitment_date, units = "days")) / 365.25,
      Sex = as.factor(Sex)
    )
  
  # Keep only baseline-free individuals
  df3 <- filter_prevalent_only(dat = df2, bs_final = bs_col_final)
  
  # Read legend of predictors to test
  df_legend <- fread(legend_path, data.table = FALSE) %>%
    # Keep:
    # - all non-disease PGS
    # - disease-specific PRS only if matched to CHD
    filter(!(disease == "yes" & trait != "CHD"))
  
  pcs <- paste0("PC", 1:5)
  results <- list()
  
  for (i in seq_len(nrow(df_legend))) {
    
    prs_raw <- paste0("PRED_eurW_", df_legend$trait[i])
    if (!(prs_raw %in% colnames(df3))) next
    
    prs_sc <- paste0(prs_raw, "_sc")
    
    df3s <- df3 %>%
      mutate(
        !!prs_sc := as.numeric(scale(.data[[prs_raw]]))
      )
    
    # Counts:
    # require time, status, and the predictor to be non-missing
    df_count <- df3s %>%
      filter(
        !is.na(time),
        !is.na(status),
        !is.na(.data[[prs_sc]])
      )
    
    N_event <- sum(df_count$status == 1, na.rm = TRUE)
    N_total <- nrow(df_count)
    N_cens  <- N_total - N_event
    
    # Model dataset:
    # complete-case on all adjustment covariates + predictor
    dfm <- df3s %>%
      filter(
        !is.na(time),
        !is.na(status),
        !is.na(Age),
        !is.na(Sex),
        !is.na(PC1), !is.na(PC2), !is.na(PC3), !is.na(PC4), !is.na(PC5),
        !is.na(.data[[prs_sc]])
      )
    
    # CHD models are explicitly adjusted for:
    # Age + Sex + PC1-PC5
    covs <- c("Age", "Sex", pcs)
    
    form <- as.formula(
      paste0("Surv(time, status) ~ ", paste(c(covs, prs_sc), collapse = " + "))
    )
    
    fit <- coxph(form, data = dfm, ties = "efron")
    s   <- summary(fit)
    
    hr <- s$conf.int[prs_sc, "exp(coef)"]
    lo <- s$conf.int[prs_sc, "lower .95"]
    up <- s$conf.int[prs_sc, "upper .95"]
    pv <- s$coefficients[prs_sc, "Pr(>|z|)"]
    
    pvf <- gsub("e\\+?(-?\\d+)", "×10\\1", formatC(pv, 2, format = "e"))
    
    prefix <- ifelse(df_legend$disease[i] == "yes", "PRS_", "PGS_")
    
    results[[length(results) + 1]] <- tibble(
      predictor         = paste0(prefix, df_legend$trait[i], "(sd)"),
      HR_CI95           = paste0(round(hr, 2), " (", round(lo, 2), "-", round(up, 2), ")"),
      p                 = pvf,
      N_event           = as.integer(N_event),
      N_censored        = as.integer(N_cens),
      N_total           = as.integer(N_total),
      N_variants        = df_legend$N_variants[i],
      study             = df_legend$study[i],
      covariates        = "Age+Sex+PC1+PC2+PC3+PC4+PC5",
      bs_col_used       = bs_col_final,
      fnf_col_used      = column_fnf,
      dataexit_col_used = column_dataexit,
      N_model           = as.integer(nrow(dfm))
    )
  }
  
  final <- bind_rows(results) %>%
    mutate(
      N_event = as.integer(N_event),
      N_censored = as.integer(N_censored),
      N_total = as.integer(N_total)
    ) %>%
    relocate(N_event, N_censored, N_total, .after = p)
  
  final
}


## ----------------------------
## 5. Run CHD analysis
## ----------------------------
results_chd <- cox_table_function_chd(
  dat = df,
  legend_path = legend_file,
  bs_col_final = "Chd_Bs_New_final"
)


## ----------------------------
## 6. Inspect results
## ----------------------------
print(results_chd)


## ----------------------------
## 7. Save results
## ----------------------------
fwrite(
  results_chd,
  outfile,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)