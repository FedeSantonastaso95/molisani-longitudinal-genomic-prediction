############################################################
# Supplementary Figure 2 - Empirical survival curves (2x2)
#
# Outcomes:
#   - Ischemic stroke (Stroke3)
#   - Lung cancer (LC)
#   - Renal cancer (RNC)
#   - Pancreatic cancer (PNC)
#
# This script:
#   1. builds endpoint-specific incident-only survival datasets
#   2. generates empirical cumulative hazard curves by PRS quintile
#      using survfit()
#   3. fits Cox models to obtain:
#        - HR per 1-SD increase in PRS
#        - HR by PRS quintile (Q2-Q5 vs Q1)
#   4. saves the summary tables used for the figure
#
# Notes:
# - Baseline-prevalent cases are removed for each endpoint.
# - For stroke, secondary events are handled upstream by setting the
#   follow-up event flag to NA when prevalent disease is already present
#   at baseline.
# - Models are adjusted for:
#     * Age
#     * PC1-PC5, when available
#     * Sex, only when appropriate
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(lubridate)
library(survival)
library(survminer)
library(ggpubr)
library(ggplot2)
library(grid)
library(scales)


## ----------------------------
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

outdir_tables  <- "path/to/output_tables"
outdir_figures <- "path/to/output_figures"

dir.create(outdir_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_figures, showWarnings = FALSE, recursive = TRUE)


## ----------------------------
## Quintile colour palette
## ----------------------------
pal_quint <- c("#2986cc", "#9fc5e8", "#999999", "#ea9999", "#cb271b")


## ----------------------------
## Robust date parser
## ----------------------------
parse_date_safe <- function(x) {
  x <- as.character(x)
  suppressWarnings(
    parse_date_time(
      x,
      orders = c("Y-m-d", "ymd", "dmy", "d-m-Y", "d/m/Y", "dmY"),
      exact = FALSE,
      quiet = TRUE
    )
  )
}


## ----------------------------
## Helpers for covariates
## ----------------------------

# Return PC1-PC5 if available
get_pc_terms <- function(df, prefix = "PC", n = 5) {
  pcs <- paste0(prefix, 1:n)
  pcs[pcs %in% names(df)]
}

# Return the most frequent non-missing value
get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Build model covariates:
# - Age always
# - PC1-PC5 when present
# - Sex only if the analysis is not sex-restricted and both sexes are present
build_covariates <- function(df, sex_filter = "all") {
  covs <- c("Age")
  
  pcs <- get_pc_terms(df, prefix = "PC", n = 5)
  covs <- c(covs, pcs)
  
  if (sex_filter == "all" && "Sex" %in% names(df)) {
    sx <- df$Sex
    if (is.factor(sx) || is.character(sx)) sx <- as.character(sx)
    if (length(unique(na.omit(sx))) > 1) {
      covs <- c(covs, "Sex")
    }
  }
  
  covs
}


## ----------------------------
## Read data and define baseline flags
## ----------------------------
raw <- fread(infile, data.table = FALSE) %>%
  mutate(
    # Baseline coding for CHD and stroke:
    #   0 = No
    #   1 = Yes, documented
    #   2 = Yes, not documented
    #
    # In the analysis, both 1 and 2 are treated as prevalent disease.
    Chd_Bs_New_final = case_when(
      Chd_Bs_New == 0 ~ 0,
      Chd_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    Cerebro_Bs_New_final = case_when(
      Cerebro_Bs_New == 0 ~ 0,
      Cerebro_Bs_New %in% c(1, 2) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Baseline cancer flags from ICD-9 diagnosis fields
    lung_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 162 | Mt2_Icd9_New == 162 | Mt3_Icd9_New == 162),
      1, 0
    ),
    renal_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 189 | Mt2_Icd9_New == 189 | Mt3_Icd9_New == 189),
      1, 0
    ),
    pancreatic_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 157 | Mt2_Icd9_New == 157 | Mt3_Icd9_New == 157),
      1, 0
    )
  ) %>%
  mutate(
    # Secondary stroke events are not counted as incident events:
    # if stroke is already prevalent at baseline and the follow-up event flag
    # is 1, the follow-up flag is set to NA.
    Stroke1_Fnof_Fup2020 = if_else(
      Cerebro_Bs_New_final == 1 & Stroke1_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Stroke1_Fnof_Fup2020)
    )
  )


## ----------------------------
## Create PRS quintiles
## ----------------------------
make_quintiles <- function(df, prs_var) {
  qcuts <- quantile(df[[prs_var]], probs = c(.2, .4, .6, .8), na.rm = TRUE)
  
  df %>%
    mutate(
      PRS_q = cut(
        .data[[prs_var]],
        breaks = c(-Inf, qcuts, Inf),
        labels = paste0("Q", 1:5),
        include.lowest = TRUE
      ),
      PRS_q = factor(PRS_q, levels = paste0("Q", 1:5))
    )
}


## ----------------------------
## Prepare survival dataset for one endpoint
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
  
  # Identify endpoint-specific date and event columns
  cols <- grep(pattern, names(df), value = TRUE)
  
  col_date <- cols[grepl("(Date|Dat)", cols, ignore.case = TRUE)][1]
  col_fnf  <- cols[
    grepl("Fno?f", cols, ignore.case = TRUE) &
      !grepl("(Date|Dat)", cols, ignore.case = TRUE)
  ][1]
  
  if (is.na(col_date) || is.na(col_fnf)) {
    stop(
      paste(
        "Cannot identify date/event columns for pattern:", pattern,
        "\nCandidates found:", paste(cols, collapse = ", ")
      )
    )
  }
  
  if (!is.null(bs_col) && !(bs_col %in% names(df))) {
    stop(paste("Baseline prevalence column not found:", bs_col))
  }
  
  # Build time-to-event variables
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
  
  # Apply sex restriction if needed
  if ("Sex" %in% names(df1)) df1 <- df1 %>% mutate(Sex = factor(Sex))
  if (sex_filter == "female") df1 <- df1 %>% filter(Sex == 0)
  if (sex_filter == "male")   df1 <- df1 %>% filter(Sex == 1)
  
  # Remove baseline-prevalent cases to define the incident-only risk set
  if (!is.null(bs_col)) {
    df1 <- df1 %>% filter(is.na(.data[[bs_col]]) | .data[[bs_col]] == 0)
  }
  
  list(
    df = df1,
    prs_var = paste0("PRED_eurW_", prs_outcome)
  )
}


## ----------------------------
## Empirical curve + model summary function
## ----------------------------
make_curve_plot_empirical <- function(df,
                                      prs_var,
                                      title_text,
                                      outcome_label = title_text,
                                      sex_filter = "all",
                                      xlim_years = c(0, 16),
                                      ylim_top = 0.015,
                                      y_mode = c("cumhaz", "event"),
                                      y_break_by = NULL) {
  
  y_mode <- match.arg(y_mode)
  
  # Add PRS quintiles
  df_q <- make_quintiles(df, prs_var)
  
  # Build adjustment covariates
  covs <- build_covariates(df_q, sex_filter = sex_filter)
  
  # Ensure PCs are numeric
  pc_terms <- covs[grepl("^PC[0-9]+$", covs)]
  if (length(pc_terms) > 0) {
    df_q <- df_q %>% mutate(across(all_of(pc_terms), as.numeric))
  }
  
  # Ensure Sex is factor if included
  if ("Sex" %in% covs) {
    df_q <- df_q %>% mutate(Sex = factor(Sex))
  }
  
  # Counts by quintile in the analysis dataset
  counts_q <- df_q %>%
    group_by(PRS_q) %>%
    summarise(
      N_total = n(),
      N_cases = sum(status == 1, na.rm = TRUE),
      N_censored = sum(status == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(outcome = outcome_label) %>%
    rename(quintile = PRS_q)
  
  # ----------------------------------------------------------
  # Continuous Cox model: HR per 1-SD increase in PRS
  # ----------------------------------------------------------
  rhs_cont <- paste(c(paste0("scale(", prs_var, ")"), covs), collapse = " + ")
  f_cont   <- as.formula(paste("Surv(time, status) ~", rhs_cont))
  fit_cont <- coxph(f_cont, data = df_q, ties = "efron")
  
  sm_cont   <- summary(fit_cont)
  ci_cont   <- sm_cont$conf.int
  coef_cont <- sm_cont$coefficients
  row_prs   <- grep("scale\\(", rownames(ci_cont))
  
  cont_stats <- data.frame(
    outcome = outcome_label,
    prs_type = "continuous",
    term = rownames(ci_cont)[row_prs],
    HR = ci_cont[row_prs, "exp(coef)"],
    CI_low = ci_cont[row_prs, "lower .95"],
    CI_high = ci_cont[row_prs, "upper .95"],
    p_value = coef_cont[row_prs, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
  # Categorical Cox model: Q2-Q5 vs Q1
  # ----------------------------------------------------------
  df_q$PRS_q <- relevel(df_q$PRS_q, ref = "Q1")
  
  rhs_cat <- paste(c("PRS_q", covs), collapse = " + ")
  f_cat   <- as.formula(paste("Surv(time, status) ~", rhs_cat))
  fit_cat <- coxph(f_cat, data = df_q, ties = "efron")
  
  sm_cat   <- summary(fit_cat)
  ci_cat   <- sm_cat$conf.int
  coef_cat <- sm_cat$coefficients
  rows_q   <- grep("^PRS_q", rownames(ci_cat))
  
  cat_stats <- data.frame(
    outcome = outcome_label,
    term = c("PRS_qQ1", rownames(ci_cat)[rows_q]),
    quintile = c("Q1", gsub("^PRS_q", "", rownames(ci_cat)[rows_q])),
    HR = c(1, ci_cat[rows_q, "exp(coef)"]),
    CI_low = c(1, ci_cat[rows_q, "lower .95"]),
    CI_high = c(1, ci_cat[rows_q, "upper .95"]),
    p_value = c(NA, coef_cat[rows_q, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(quintile = factor(quintile, levels = paste0("Q", 1:5))) %>%
    arrange(quintile) %>%
    left_join(counts_q, by = c("outcome", "quintile"))
  
  # ----------------------------------------------------------
  # Empirical survival curves by PRS quintile
  # ----------------------------------------------------------
  fit_emp <- survfit(Surv(time, status) ~ PRS_q, data = df_q)
  
  # Legend labels: HR only, as requested
  hr_labs <- c(
    "Q1 (HR = 1.00 [ref])",
    paste0("Q2 (HR = ", sprintf("%.2f", cat_stats$HR[cat_stats$quintile == "Q2"]), ")"),
    paste0("Q3 (HR = ", sprintf("%.2f", cat_stats$HR[cat_stats$quintile == "Q3"]), ")"),
    paste0("Q4 (HR = ", sprintf("%.2f", cat_stats$HR[cat_stats$quintile == "Q4"]), ")"),
    paste0("Q5 (HR = ", sprintf("%.2f", cat_stats$HR[cat_stats$quintile == "Q5"]), ")")
  )
  
  plt <- ggsurvplot(
    fit_emp,
    data = df_q,
    fun = if (y_mode == "event") "event" else "cumhaz",
    conf.int = TRUE,
    conf.int.alpha = 0.08,
    xlim = xlim_years,
    palette = pal_quint,
    legend.title = "PRS risk category",
    legend.labs = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    ggtheme = theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18, family = "Arial")
      ),
    title = title_text,
    xlab = "",
    ylab = "",
    font.legend = 12,
    font.x = 15,
    font.y = 15,
    font.tickslab = 14,
    censor = FALSE
  )
  
  if (is.null(y_break_by)) y_break_by <- ylim_top / 4
  
  plt$plot <- plt$plot +
    scale_y_continuous(
      limits = c(0, ylim_top),
      breaks = seq(0, ylim_top, by = y_break_by),
      labels = function(x) ifelse(x == 0, "0", formatC(x, format = "f", digits = 3))
    ) +
    scale_x_continuous(breaks = seq(xlim_years[1], xlim_years[2], by = 4)) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 1)),
      fill = "none"
    ) +
    scale_colour_manual(values = pal_quint, labels = hr_labs) +
    theme(
      text = element_text(family = "Arial"),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16)
    )
  
  s16 <- summary(fit_emp, times = xlim_years[2], extend = TRUE)
  
  if (y_mode == "event") {
    est16 <- setNames(s16$surv, nm = paste0("Q", 1:5))
  } else {
    est16 <- setNames(s16$cumhaz, nm = paste0("Q", 1:5))
  }
  
  list(
    plot = plt$plot,
    est16 = est16,
    cont_stats = cont_stats,
    cat_stats = cat_stats
  )
}


## ----------------------------
## Endpoint specification for Supplementary Figure 2
## ----------------------------
endpoints <- data.frame(
  column_name = c("Stroke3", "Lc", "Rnc", "Pnc"),
  sex_stratified = c("no", "no", "no", "no"),
  title = c(
    "Ischemic stroke",
    "Lung cancer",
    "Renal cancer",
    "Pancreatic cancer"
  ),
  outcome = c("stroke", "lung", "renal", "pancreatic"),
  acro = c("Stroke3", "LC", "RNC", "PNC"),
  bs_col = c(
    "Cerebro_Bs_New_final",
    "lung_Bs_final",
    "renal_Bs_final",
    "pancreatic_Bs_final"
  ),
  stringsAsFactors = FALSE
)

# Figure order
endpoints_4b <- endpoints %>%
  mutate(column_name = factor(column_name, levels = c("Stroke3", "Lc", "Rnc", "Pnc"))) %>%
  arrange(column_name)


## ----------------------------
## Run analysis and collect outputs
## ----------------------------
plot_list_4b        <- vector("list", nrow(endpoints_4b))
results_cont_list_b <- vector("list", nrow(endpoints_4b))
results_cat_list_b  <- vector("list", nrow(endpoints_4b))

Y_MAX  <- 0.015
Y_STEP <- 0.005

for (i in seq_len(nrow(endpoints_4b))) {
  
  pat      <- as.character(endpoints_4b$column_name[i])
  outc     <- endpoints_4b$outcome[i]
  acr      <- endpoints_4b$acro[i]
  ttl_full <- endpoints_4b$title[i]
  bsc      <- endpoints_4b$bs_col[i]
  
  sex_f <- "all"
  
  prep <- prep_surv_df(
    df = raw,
    pattern = pat,
    prs_outcome = outc,
    sex_filter = sex_f,
    bs_col = bsc
  )
  
  out <- make_curve_plot_empirical(
    prep$df,
    prep$prs_var,
    title_text = ttl_full,
    outcome_label = acr,
    sex_filter = sex_f,
    y_mode = "cumhaz",
    ylim_top = Y_MAX,
    y_break_by = Y_STEP
  )
  
  plot_list_4b[[i]]      <- out$plot
  results_cont_list_b[[i]] <- out$cont_stats
  results_cat_list_b[[i]]  <- out$cat_stats
}

HR_continuous_b <- do.call(rbind, results_cont_list_b)
HR_quintiles_b  <- do.call(rbind, results_cat_list_b)

print(HR_quintiles_b)


## ----------------------------
## Save summary tables
## ----------------------------
outfile_hr_quint <- file.path(
  outdir_tables,
  paste0("SuppFig2_HR_CI_p_cases_censored_QUINTILES_Stroke3_LC_RNC_PNC_ymax", Y_MAX, ".txt")
)

fwrite(
  HR_quintiles_b,
  outfile_hr_quint,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

outfile_hr_cont <- file.path(
  outdir_tables,
  paste0("SuppFig2_HR_continuous_perSD_Stroke3_LC_RNC_PNC_ymax", Y_MAX, ".txt")
)

fwrite(
  HR_continuous_b,
  outfile_hr_cont,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

cat("Saved HR quintile table:", outfile_hr_quint, "\n")
cat("Saved HR continuous table:", outfile_hr_cont, "\n")


## ----------------------------
## Assemble Supplementary Figure 2
## ----------------------------
final_plot_b <- ggarrange(
  plotlist = plot_list_4b,
  ncol = 2,
  nrow = 2
)

final_plot_b <- annotate_figure(
  final_plot_b,
  bottom = text_grob("Follow-up time (years)", size = 18, face = "bold", family = "Arial"),
  left   = text_grob("Cumulative hazard", size = 18, face = "bold", family = "Arial", rot = 90)
)

print(final_plot_b)


## ----------------------------
## Save Supplementary Figure 2
## ----------------------------
png(
  file.path(
    outdir_figures,
    paste0("SupplementaryFigure2_SurvivalPlots_Stroke3_LC_RNC_PNC_ymax", Y_MAX, ".png")
  ),
  width = 8.3,
  height = 8.3,
  units = "in",
  res = 600
)
print(final_plot_b)
dev.off()