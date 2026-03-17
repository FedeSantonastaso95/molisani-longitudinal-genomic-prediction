############################################################
# Sex interaction and sex-stratified CHD survival analyses
#
# Description:
# This script performs two sets of analyses for incident CHD:
#
#   1. Interaction models:
#      Cox models including Sex * PRS interaction terms
#
#   2. Sex-stratified models:
#      Separate Cox models in women and men for each PRS
#
# Analysis design:
#   - incident-only CHD risk set
#   - baseline-prevalent CHD cases are excluded
#   - secondary CHD events are handled upstream by setting the follow-up
#     CHD event flag to NA when prevalent CHD is already present at baseline
#
# Model adjustment:
#   - Interaction models: Age + Sex + PRS + Sex*PRS
#   - Sex-stratified models: Age + PC1-PC5 + PRS
#
# Output:
#   - one table for interaction results
#   - one table for sex-stratified results
#   - counts of CHD cases and censored individuals in women and men
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
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

outfile_interaction <- "path/to/output/2026.02.12_interactionSex_cox_results.txt"
outfile_sex_stratified <- "path/to/output/2026.02.12_SexStratified.txt"


## ----------------------------
## 1. Read and preprocess data
## ----------------------------
df <- fread(
  infile,
  data.table = FALSE
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
## 2. Traits (PRS/PGS) to include
## ----------------------------
traits <- c("CHD", "ALT", "AST", "ApoB", "BMI", "Uric", "WHR", "SBP", "DBP", "HDL", "LDL", "TRI")
traits_pattern <- paste(traits, collapse = "|")


## ----------------------------
## 3. Identify CHD follow-up columns
## ----------------------------
column_name <- "Chd0"

columns_disease <- grep(column_name, colnames(df), value = TRUE)
column_dataexit <- columns_disease[grepl("Dat", columns_disease)]
column_fnf <- columns_disease[grepl("Fn", columns_disease)]
column_fnf <- column_fnf[!grepl("Dat|Annipers", column_fnf)]

if (length(column_dataexit) != 1) {
  stop("Ambiguous CHD date column: ", paste(column_dataexit, collapse = ", "))
}
if (length(column_fnf) != 1) {
  stop("Ambiguous CHD event column: ", paste(column_fnf, collapse = ", "))
}


## ----------------------------
## 4. Build time-to-event dataset
## ----------------------------
df1 <- df %>%
  mutate(
    # Binary event coding for Surv():
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
## 5. Add standardized PRS/PGS variables
## ----------------------------
df3 <- df2

for (tr in traits) {
  prs_col <- paste0("PRED_eurW_", tr)
  prs_sc  <- paste0(prs_col, "_scaled")
  
  if (prs_col %in% colnames(df3)) {
    df3[[prs_sc]] <- as.numeric(scale(df3[[prs_col]]))
  }
}

predictor_columns <- grep(
  paste0("PRED_eurW_(", traits_pattern, ")_scaled$"),
  colnames(df3),
  value = TRUE
)

df4 <- df3


## ----------------------------
## 6. Interaction models: Sex * PRS
## ----------------------------
#
# Interaction terms are tested for each predictor separately.
# Models are adjusted for Age.
#
# If you want, PCs can also be added here, but this version keeps the
# interaction model aligned with your original analysis structure.
pred_short <- c(
  "CHD", "HDL", "UA", "AST", "BMI", "DBP", "ALT",
  "TG", "ApoB", "WHR", "LDL", "SBP"
)

map_pred <- c(
  CHD  = "PRED_eurW_CHD_scaled",
  SBP  = "PRED_eurW_SBP_scaled",
  DBP  = "PRED_eurW_DBP_scaled",
  BMI  = "PRED_eurW_BMI_scaled",
  WHR  = "PRED_eurW_WHR_scaled",
  LDL  = "PRED_eurW_LDL_scaled",
  HDL  = "PRED_eurW_HDL_scaled",
  TG   = "PRED_eurW_TRI_scaled",
  ApoB = "PRED_eurW_ApoB_scaled",
  ALT  = "PRED_eurW_ALT_scaled",
  AST  = "PRED_eurW_AST_scaled",
  UA   = "PRED_eurW_Uric_scaled"
)

pred_cols <- unname(map_pred[pred_short])

if (any(is.na(pred_cols))) {
  stop("Missing mapping for: ", paste(pred_short[is.na(pred_cols)], collapse = ", "))
}

missing <- setdiff(pred_cols, names(df4))
if (length(missing) > 0) {
  stop("These predictor columns are not present in df4: ", paste(missing, collapse = ", "))
}

interaction_results <- vector("list", length(pred_cols))
names(interaction_results) <- pred_cols

for (i in seq_along(pred_cols)) {
  
  pred <- pred_cols[i]
  pred_label <- names(map_pred)[match(pred, map_pred)]
  pred_bt <- paste0("`", pred, "`")
  
  dtmp <- df4 %>%
    filter(
      !is.na(time),
      !is.na(status),
      !is.na(Age),
      !is.na(Sex),
      !is.na(.data[[pred]])
    )
  
  f <- as.formula(paste0("Surv(time, status) ~ Age + Sex * ", pred_bt))
  fit <- coxph(f, data = dtmp)
  
  coefs <- as.data.frame(summary(fit)$coefficients) %>%
    tibble::rownames_to_column("variable")
  colnames(coefs) <- c("variable", "beta", "HR", "SE", "Z", "pval")
  
  cis <- as.data.frame(summary(fit)$conf.int) %>%
    dplyr::select(-`exp(-coef)`) %>%
    tibble::rownames_to_column("variable")
  colnames(cis) <- c("variable", "HR", "CI95_lower", "CI95_upper")
  
  out <- coefs %>%
    left_join(cis, by = c("variable", "HR")) %>%
    select(variable, HR, CI95_lower, CI95_upper, pval) %>%
    mutate(
      predictor = pred_label,
      model = pred
    )
  
  interaction_results[[i]] <- out
}

table_final_interaction <- bind_rows(interaction_results) %>%
  mutate(
    pval = format(pval, scientific = TRUE, digits = 3)
  )

print(table_final_interaction)

fwrite(
  table_final_interaction,
  outfile_interaction,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


## ----------------------------
## 7. Sex-stratified univariate Cox models
## ----------------------------
#
# Separate models are run in women and men.
# Each model is adjusted for:
#   - Age
#   - PC1-PC5
#   - one standardized PRS/PGS predictor
df_female <- df4 %>% filter(Sex == 0)
df_male   <- df4 %>% filter(Sex == 1)

run_univariate_cox <- function(data, predictors) {
  all_results <- list()
  
  for (pred in predictors) {
    
    dtmp <- data %>%
      filter(
        !is.na(time),
        !is.na(status),
        !is.na(Age),
        !is.na(PC1),
        !is.na(PC2),
        !is.na(PC3),
        !is.na(PC4),
        !is.na(PC5),
        !is.na(.data[[pred]])
      )
    
    full_formula <- paste(
      "Surv(time, status) ~ Age + PC1 + PC2 + PC3 + PC4 + PC5 +",
      pred
    )
    
    fit <- coxph(as.formula(full_formula), data = dtmp)
    
    coef_table <- as.data.frame(summary(fit)$coefficients)
    colnames(coef_table) <- c("beta", "HR", "SE", "Z", "pval")
    coef_table <- tibble::rownames_to_column(coef_table, var = "variable")
    
    conf_int <- as.data.frame(summary(fit)$conf.int) %>%
      dplyr::select(-"exp(-coef)")
    colnames(conf_int) <- c("HR", "CI95_lower", "CI95_upper")
    conf_int <- tibble::rownames_to_column(conf_int, var = "variable")
    
    table_combined <- coef_table %>%
      dplyr::select(-beta, -SE, -Z) %>%
      full_join(conf_int, by = c("variable", "HR")) %>%
      dplyr::select(variable, HR, CI95_lower, CI95_upper, pval)
    
    all_results[[pred]] <- table_combined
  }
  
  bind_rows(all_results) %>%
    filter(variable != "Age") %>%
    mutate(
      variable = gsub("PRED_eurW_|_scaled", "", variable),
      variable = ifelse(variable == "Uric", "UA", variable),
      variable = ifelse(variable == "TRI", "TG", variable)
    )
}

uni_female <- run_univariate_cox(df_female, predictor_columns) %>%
  mutate(group = "Women")

uni_male <- run_univariate_cox(df_male, predictor_columns) %>%
  mutate(group = "Men")

final_uni <- bind_rows(uni_female, uni_male)


## ----------------------------
## 8. Add counts of CHD cases and censored individuals
## ----------------------------
counts_by_group <- df4 %>%
  mutate(group = ifelse(Sex == 0, "Women", "Men")) %>%
  group_by(group) %>%
  summarise(
    CHD_cases  = sum(status == 1, na.rm = TRUE),
    N_censored = sum(status == 0, na.rm = TRUE),
    .groups = "drop"
  )

get_val <- function(df_counts, grp, col) {
  v <- df_counts %>% filter(group == grp) %>% pull(!!sym(col))
  if (length(v) == 0) NA_integer_ else as.integer(v[1])
}

CHD_cases_Women  <- get_val(counts_by_group, "Women", "CHD_cases")
N_censored_Women <- get_val(counts_by_group, "Women", "N_censored")
CHD_cases_Men    <- get_val(counts_by_group, "Men", "CHD_cases")
N_censored_Men   <- get_val(counts_by_group, "Men", "N_censored")

final_uni <- final_uni %>%
  mutate(
    CHD_cases_Women  = CHD_cases_Women,
    N_censored_Women = N_censored_Women,
    CHD_cases_Men    = CHD_cases_Men,
    N_censored_Men   = N_censored_Men
  )


## ----------------------------
## 9. Reorder and save sex-stratified results
## ----------------------------
final_uni <- final_uni %>%
  group_by(group) %>%
  arrange(HR, .by_group = TRUE) %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  ungroup()

all_out <- final_uni

fwrite(
  all_out,
  outfile_sex_stratified,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

print(all_out)