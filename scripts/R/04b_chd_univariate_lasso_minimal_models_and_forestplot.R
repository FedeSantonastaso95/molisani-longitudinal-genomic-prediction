############################################################
# 04b_chd_univariate_lasso_minimal_models_and_forestplot.R
#
# Description:
# This script performs incident CHD survival analyses using:
#   1. univariate Cox models
#   2. LASSO penalised Cox regression for variable selection
#   3. a refitted Cox model including predictors retained by LASSO
#   4. a minimal Cox model including selected predictors of interest
#
# Outcome:
#   - incident Coronary Heart Disease (CHD)
#
# Predictors tested:
#   - the CHD-specific polygenic risk score (PRS)
#   - polygenic scores (PGS) for intermediate endophenotypes / risk factors
#
# Predictor selection:
#   The predictors included in this script were pre-selected based on results
#   from the screening analysis:
#     "04a_chd_survival_prs_and_endophenotype_pgs.R"
#
#   Specifically, only predictors associated with incident CHD at:
#     P ≤ 0.001
#   were retained.
#
#   This threshold was derived from the effective number of independent
#   predictors (n = 36), corresponding to a multiple-testing correction.
#
# Analysis design:
#   - incident-only risk set
#   - baseline-prevalent CHD cases are excluded
#   - secondary CHD events are handled upstream by setting the follow-up
#     CHD event flag to NA when prevalent CHD is already present at baseline
#
# Model adjustment:
#   - Univariate Cox models are adjusted for Age, Sex, and PC1-PC5
#   - LASSO is applied to the standardized PRS/PGS predictors only
#   - predictors retained by LASSO are then refitted in a standard Cox model
#     adjusted for Age, Sex, and PC1-PC5
#   - the minimal Cox model is adjusted for Age, Sex, and PC1-PC5
#
# Output:
#   - summary tables for univariate, LASSO-refit, and minimal models
#   - compact result tables for manuscript text
#   - forest plot comparing the three modelling approaches
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(survival)
library(dplyr)
library(rlang)
library(data.table)
library(tibble)
library(glmnet)
library(ggplot2)
library(grid)


## ----------------------------
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

outdir_tables <- "path/to/output_tables"
outdir_figures <- "path/to/output_figures"

dir.create(outdir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_figures, showWarnings = FALSE, recursive = TRUE)


## ----------------------------
## Predictors retained from the prior screening step
## ----------------------------
#
# These traits correspond to PRS/PGS passing:
#   P ≤ 0.001
# in the screening script:
#   "04a_chd_survival_prs_and_endophenotype_pgs.R"
#
# The threshold was based on the effective number of independent
# predictors (n = 36).
traits <- c(
  "CHD",
  "ALT",
  "AST",
  "ApoB",
  "BMI",
  "Uric",
  "WHR",
  "SBP",
  "DBP",
  "HDL",
  "LDL",
  "TRI"
)


## ----------------------------
## 1. Read and preprocess data
## ----------------------------
df <- fread(infile, data.table = FALSE) %>%
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
## 2. Identify CHD follow-up date and event columns
## ----------------------------
columns_disease <- grep("Chd0", colnames(df), value = TRUE)

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
## 3. Build the incident CHD survival dataset
## ----------------------------
#
# Event coding:
#   status = 1 if incident CHD occurred
#   status = 0 otherwise
#
# Follow-up time is computed from recruitment to event/censoring.
df1 <- df %>%
  mutate(
    status = as.numeric(.data[[column_fnf]] == 1)
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
    Sex = as.factor(Sex)
  ) %>%
  filter(!is.na(time), time >= 0)

# Keep an incident-only risk set by excluding baseline-prevalent CHD
df2 <- df2 %>%
  filter(Chd_Bs_New_final == 0)


## ----------------------------
## 4. Keep selected PRS/PGS predictors and standardize them
## ----------------------------
traits_pattern <- paste(traits, collapse = "|")

df3 <- df2 %>%
  select(-matches(paste0("PRED_eurW_(?!(", traits_pattern, "))"), perl = TRUE))

for (tr in traits) {
  pred_name <- paste0("PRED_eurW_", tr)
  if (pred_name %in% colnames(df3)) {
    df3[[paste0(pred_name, "_scaled")]] <- as.numeric(scale(df3[[pred_name]]))
  }
}

predictor_columns <- grep(
  paste0("PRED_eurW_(", traits_pattern, ")_scaled$"),
  colnames(df3),
  value = TRUE
)

# Complete-case dataset for models adjusted for Age + Sex + PC1-PC5
df_model <- df3 %>%
  filter(
    !is.na(time),
    !is.na(status),
    !is.na(Age),
    !is.na(Sex),
    !is.na(PC1),
    !is.na(PC2),
    !is.na(PC3),
    !is.na(PC4),
    !is.na(PC5)
  )


## ----------------------------
## 5. Univariate Cox models
## Each model: Surv(time, status) ~ Age + Sex + PC1-PC5 + predictor
## ----------------------------
univariate_results <- list()

for (pred in predictor_columns) {
  dtmp <- df_model %>% filter(!is.na(.data[[pred]]))
  
  fit <- coxph(
    as.formula(
      paste(
        "Surv(time, status) ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 +",
        pred
      )
    ),
    data = dtmp
  )
  
  s <- summary(fit)
  row_idx <- which(rownames(s$coefficients) == pred)
  
  univariate_results[[pred]] <- tibble(
    variable = pred,
    HR = s$conf.int[row_idx, "exp(coef)"],
    CI95_lower = s$conf.int[row_idx, "lower .95"],
    CI95_upper = s$conf.int[row_idx, "upper .95"],
    pval = s$coefficients[row_idx, "Pr(>|z|)"],
    model = "Univariate"
  )
}

final_table_univariate <- bind_rows(univariate_results)


## ----------------------------
## 6. LASSO penalised Cox regression
## ----------------------------
#
# LASSO is applied to the standardized PRS/PGS predictors only.
# Predictors retained at lambda.min are subsequently refitted in a
# standard Cox model adjusted for Age, Sex, and PC1-PC5.
df_lasso <- df_model %>%
  select(time, status, all_of(predictor_columns)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

x_lasso <- as.matrix(df_lasso[, predictor_columns, drop = FALSE])
y_lasso <- Surv(df_lasso$time, df_lasso$status)

set.seed(123)
lasso_fit <- cv.glmnet(
  x = x_lasso,
  y = y_lasso,
  family = "cox",
  alpha = 1,
  nfolds = 10
)

best_lambda <- lasso_fit$lambda.min

coef_matrix <- coef(lasso_fit, s = "lambda.min")
coef_values <- as.matrix(coef_matrix)

selected_vars <- rownames(coef_values)[coef_values != 0]

if (length(selected_vars) == 0) {
  stop("No predictors retained by LASSO at lambda.min.")
}

# Refit standard Cox model on selected predictors, adjusted for covariates
df_lasso_refit <- df_model %>%
  select(time, status, Age, Sex, PC1, PC2, PC3, PC4, PC5, all_of(selected_vars)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

lasso_refit_formula <- as.formula(
  paste(
    "Surv(time, status) ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 +",
    paste(selected_vars, collapse = " + ")
  )
)

lasso_refit <- coxph(
  lasso_refit_formula,
  data = df_lasso_refit
)

s_lasso <- summary(lasso_refit)
rows_keep_lasso <- rownames(s_lasso$coefficients) %in% selected_vars

final_table_lasso <- tibble(
  variable = rownames(s_lasso$coefficients)[rows_keep_lasso],
  HR = s_lasso$conf.int[rows_keep_lasso, "exp(coef)"],
  CI95_lower = s_lasso$conf.int[rows_keep_lasso, "lower .95"],
  CI95_upper = s_lasso$conf.int[rows_keep_lasso, "upper .95"],
  pval = s_lasso$coefficients[rows_keep_lasso, "Pr(>|z|)"],
  model = "LASSO"
)


## ----------------------------
## 7. Minimal Cox model
## ----------------------------
#
# This model includes the key predictors retained after LASSO and is
# adjusted for Age, Sex, and PC1-PC5.
minimal_predictors <- c(
  "PRED_eurW_CHD_scaled",
  "PRED_eurW_ApoB_scaled",
  "PRED_eurW_SBP_scaled"
)

minimal_predictors <- minimal_predictors[minimal_predictors %in% colnames(df_model)]

df_minimal <- df_model %>%
  select(
    time, status, Age, Sex, PC1, PC2, PC3, PC4, PC5,
    all_of(minimal_predictors)
  ) %>%
  filter(if_all(everything(), ~ !is.na(.)))

minimal_formula <- as.formula(
  paste(
    "Surv(time, status) ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 +",
    paste(minimal_predictors, collapse = " + ")
  )
)

minimal_fit <- coxph(minimal_formula, data = df_minimal)
s_minimal <- summary(minimal_fit)

rows_keep <- rownames(s_minimal$coefficients) %in% minimal_predictors

final_table_minimal <- tibble(
  variable = rownames(s_minimal$coefficients)[rows_keep],
  HR = s_minimal$conf.int[rows_keep, "exp(coef)"],
  CI95_lower = s_minimal$conf.int[rows_keep, "lower .95"],
  CI95_upper = s_minimal$conf.int[rows_keep, "upper .95"],
  pval = s_minimal$coefficients[rows_keep, "Pr(>|z|)"],
  model = "Minimal"
)


## ----------------------------
## 8. Save core model tables
## ----------------------------
fwrite(
  final_table_univariate,
  file.path(outdir_tables, "CHD_univariate_cox_models.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  final_table_lasso,
  file.path(outdir_tables, "CHD_lasso_refit_cox_models.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  final_table_minimal,
  file.path(outdir_tables, "CHD_minimal_cox_model.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Save LASSO tuning and selected variables
lasso_selection_info <- tibble(
  lambda_min = best_lambda,
  selected_variable = selected_vars
)

fwrite(
  lasso_selection_info,
  file.path(outdir_tables, "CHD_lasso_selection_info.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


## ----------------------------
## 9. Extract compact results for manuscript text
## ----------------------------

# LASSO-refit results used in the Results section
lasso_results_for_text <- final_table_lasso %>%
  filter(variable %in% c(
    "PRED_eurW_CHD_scaled",
    "PRED_eurW_SBP_scaled",
    "PRED_eurW_ApoB_scaled"
  )) %>%
  mutate(
    predictor = gsub("PRED_eurW_|_scaled", "", variable),
    predictor = ifelse(predictor == "CHD", "PRSCHD", paste0("PGI", predictor)),
    HR_CI_text = paste0(
      sprintf("%.2f", HR),
      " (95% CI ",
      sprintf("%.2f", CI95_lower),
      "–",
      sprintf("%.2f", CI95_upper),
      ")"
    ),
    p_text = gsub("e", "×10^", formatC(pval, format = "e", digits = 2)),
    adjustment = "Age + Sex + PC1-PC5"
  ) %>%
  select(predictor, HR_CI_text, p_text, adjustment)

fwrite(
  lasso_results_for_text,
  file.path(outdir_tables, "CHD_lasso_results_for_text.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Minimal model results used in the Results section
minimal_results_for_text <- final_table_minimal %>%
  mutate(
    predictor = gsub("PRED_eurW_|_scaled", "", variable),
    predictor = ifelse(predictor == "CHD", "PRSCHD", paste0("PGI", predictor)),
    HR_CI_text = paste0(
      sprintf("%.2f", HR),
      " (95% CI ",
      sprintf("%.2f", CI95_lower),
      "–",
      sprintf("%.2f", CI95_upper),
      ")"
    ),
    p_text = gsub("e", "×10^", formatC(pval, format = "e", digits = 2)),
    adjustment = "Age + Sex + PC1-PC5"
  ) %>%
  select(predictor, HR_CI_text, p_text, adjustment)

fwrite(
  minimal_results_for_text,
  file.path(outdir_tables, "CHD_minimal_results_for_text.txt"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


## ----------------------------
## 10. Combine results for the forest plot
## ----------------------------
final_table <- bind_rows(
  final_table_univariate,
  final_table_lasso,
  final_table_minimal
) %>%
  mutate(
    variable_clean = gsub("PRED_eurW_|_scaled", "", variable),
    variable_clean = ifelse(variable_clean == "Uric", "UA", variable_clean),
    variable_clean = ifelse(variable_clean == "TRI", "TG", variable_clean),
    dot_shape = ifelse(variable_clean == "CHD", "bianco", "nero"),
    model = factor(model, levels = c("Univariate", "LASSO", "Minimal"))
  ) %>%
  group_by(model) %>%
  arrange(HR, .by_group = TRUE) %>%
  mutate(variable_clean = factor(variable_clean, levels = unique(variable_clean))) %>%
  ungroup()

final_table <- final_table %>%
  mutate(variable_clean = factor(variable_clean, levels = rev(unique(variable_clean))))


## ----------------------------
## 11. Forest plot
## ----------------------------
forest_plot <- ggplot(
  final_table,
  aes(y = variable_clean, x = HR, xmin = CI95_lower, xmax = CI95_upper)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "#9a9a9a") +
  geom_errorbarh(
    height = 0.18,
    linewidth = 0.7,
    color = "black"
  ) +
  geom_point(
    aes(shape = dot_shape, fill = dot_shape),
    size = 3.2,
    stroke = 1.0,
    color = "black"
  ) +
  facet_wrap(~ model, ncol = 3, scales = "fixed") +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(bianco = 22, nero = 22),
    labels = c(bianco = "Disease", nero = "Risk factor")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(bianco = "white", nero = "black"),
    labels = c(bianco = "Disease", nero = "Risk factor")
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 6, stroke = 1.2)),
    fill = guide_legend(override.aes = list(size = 6, color = "black", stroke = 1.2))
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.ticks.length = unit(2.5, "mm"),
    axis.ticks.x = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.5, color = "black"),
    axis.text.x = element_text(size = 10, margin = margin(t = 3)),
    axis.text.y = element_text(size = 10, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 13, face = "bold", margin = margin(b = 3)),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.height = unit(6, "mm"),
    legend.key.width = unit(6, "mm"),
    legend.spacing.y = unit(2, "mm"),
    legend.margin = margin(2, 2, 2, 2),
    panel.spacing = unit(10, "pt")
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.08))
  )

ggsave(
  plot = forest_plot,
  filename = file.path(outdir_figures, "CHD_univariate_lasso_minimal_forestplot.png"),
  dpi = 300,
  width = 9,
  height = 5,
  units = "in",
  device = "png"
)