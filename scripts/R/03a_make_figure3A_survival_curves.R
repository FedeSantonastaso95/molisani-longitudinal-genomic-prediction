############################################################
# 03a_make_figure3A_survival_curves.R
#
# Description:
# This script generates the empirical cumulative incidence curves used in
# Figure 3A of the manuscript for four outcomes:
#   - Coronary Heart Disease (CHD)
#   - Breast cancer (BC)
#   - Prostate cancer (PC)
#   - Colorectal cancer (CRC)
#
# Curves are drawn by PRS quintiles using empirical survfit() estimates.
# The plotted quantity is cumulative incidence, obtained as:
#   1 - S(t)
# by setting fun = "event" in ggsurvplot().
#
# In addition to the figure, the script exports summary tables containing:
#   - categorical HRs by PRS quintile (Q2-Q5 vs Q1)
#   - continuous HRs per 1-SD increase in PRS
#   - trend HRs per one-quintile increase
#   - rate tables (N, cases, person-years, rate)
#
# Notes:
# - Baseline-prevalent cases are excluded from each endpoint-specific risk set.
# - For CHD, secondary events are handled upstream by recoding the follow-up
#   event flag to NA when prevalent disease is already present at baseline.
# - Models are adjusted for:
#     * Age
#     * Sex, when appropriate
#     * PC1-PC5
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


## ----------------------------
## Input / output paths
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"

outdir_tables <- "path/to/output_tables"
outdir_plots  <- "path/to/output_figures"

dir.create(outdir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_plots,  showWarnings = FALSE, recursive = TRUE)


## ----------------------------
## PRS quintile colour palette
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
## Covariate helpers
## ----------------------------

# Return PC1-PC5 if present in the dataset
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

# Build the adjustment set:
# - Age always
# - PC1-PC5 when available
# - Sex only when the analysis is not sex-restricted and Sex varies
build_covariates <- function(df, sex_filter = "all") {
  covs <- c("Age")
  
  pcs <- get_pc_terms(df, "PC", 5)
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
    
    # Baseline cancer status from ICD-9 diagnosis fields
    breast_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 174 | Mt2_Icd9_New == 174 | Mt3_Icd9_New == 174),
      1, 0
    ),
    prostate_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 &
        (Mt1_Icd9 == 185 | Mt2_Icd9_New == 185 | Mt3_Icd9_New == 185),
      1, 0
    ),
    colorectal_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 &
        (
          Mt1_Icd9 %in% c(153, 154) |
            Mt2_Icd9_New %in% c(153, 154) |
            Mt3_Icd9_New %in% c(153, 154)
        ),
      1, 0
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
  
  # Identify endpoint-specific date and event columns
  cols <- grep(pattern, names(df), value = TRUE)
  
  col_date <- cols[grepl("(Date|Dat)", cols, ignore.case = TRUE)][1]
  col_fnf  <- cols[
    grepl("Fno?f", cols, ignore.case = TRUE) &
      !grepl("(Date|Dat)", cols, ignore.case = TRUE)
  ][1]
  
  # Explicit fallback for CHD naming
  if ((is.na(col_date) || is.na(col_fnf)) && pattern == "Chd0") {
    if ("Chd0_Fnof_Date_Fup2020" %in% names(df)) col_date <- "Chd0_Fnof_Date_Fup2020"
    if ("Chd0_Fnof_Fup2020" %in% names(df))      col_fnf  <- "Chd0_Fnof_Fup2020"
  }
  
  if (is.na(col_date) || is.na(col_fnf)) {
    stop(paste("Cannot identify date/event columns for pattern:", pattern))
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
  
  # Apply sex restriction for sex-specific outcomes
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
## Fit empirical curves and extract model summaries
## ----------------------------
make_curve_plot_empirical <- function(df,
                                      prs_var,
                                      title_text,
                                      outcome_label = title_text,
                                      sex_filter = "all",
                                      xlim_years = c(0, 16),
                                      ylim_top = 0.08,
                                      y_mode = c("event", "cumhaz")) {
  
  y_mode <- match.arg(y_mode)
  
  # Add PRS quintiles
  df_q <- make_quintiles(df, prs_var)
  
  # Build adjustment covariates
  covs <- build_covariates(df_q, sex_filter = sex_filter)
  
  # Ensure numeric PCs
  pc_terms <- covs[grepl("^PC[0-9]+$", covs)]
  if (length(pc_terms) > 0) {
    df_q <- df_q %>% mutate(across(all_of(pc_terms), as.numeric))
  }
  
  # Ensure Sex is factor if included
  if ("Sex" %in% covs) {
    df_q <- df_q %>% mutate(Sex = factor(Sex))
  }
  
  # ==========================================================
  # Counts by quintile
  # ==========================================================
  count_table <- df_q %>%
    group_by(PRS_q) %>%
    summarise(
      N_total = n(),
      N_cases = sum(status == 1, na.rm = TRUE),
      N_censored = sum(status == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(quintile = PRS_q) %>%
    mutate(outcome = outcome_label)
  
  # ==========================================================
  # Continuous Cox model: HR per 1-SD increase in PRS
  # ==========================================================
  rhs_cont <- paste(c(paste0("scale(", prs_var, ")"), covs), collapse = " + ")
  fit_cont <- coxph(as.formula(paste("Surv(time, status) ~", rhs_cont)), data = df_q, ties = "efron")
  sm_cont  <- summary(fit_cont)
  
  ci_cont   <- sm_cont$conf.int
  coef_cont <- sm_cont$coefficients
  row_prs   <- grep("scale\\(", rownames(ci_cont))
  
  HR_continuous <- data.frame(
    outcome = outcome_label,
    term = rownames(ci_cont)[row_prs],
    HR = ci_cont[row_prs, "exp(coef)"],
    CI_low = ci_cont[row_prs, "lower .95"],
    CI_high = ci_cont[row_prs, "upper .95"],
    p_value = coef_cont[row_prs, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  # ==========================================================
  # Trend model: HR per one-quintile increase
  # ==========================================================
  df_q <- df_q %>% mutate(PRS_q_num = as.numeric(PRS_q))
  rhs_trend <- paste(c("PRS_q_num", covs), collapse = " + ")
  fit_trend <- coxph(as.formula(paste("Surv(time, status) ~", rhs_trend)), data = df_q, ties = "efron")
  sm_trend  <- summary(fit_trend)
  
  HR_trend <- data.frame(
    outcome = outcome_label,
    term = "PRS_q_num (per quintile increase)",
    HR = exp(coef(fit_trend)["PRS_q_num"]),
    CI_low = exp(confint(fit_trend)["PRS_q_num", 1]),
    CI_high = exp(confint(fit_trend)["PRS_q_num", 2]),
    p_value = sm_trend$coefficients["PRS_q_num", "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  # ==========================================================
  # Rate table
  # ==========================================================
  rate_table <- df_q %>%
    group_by(PRS_q) %>%
    summarise(
      N = n(),
      cases = sum(status == 1, na.rm = TRUE),
      py = sum(time, na.rm = TRUE),
      rate = cases / py,
      .groups = "drop"
    ) %>%
    mutate(outcome = outcome_label) %>%
    select(outcome, PRS_q, N, cases, py, rate)
  
  # ==========================================================
  # Categorical Cox model: Q2-Q5 vs Q1
  # ==========================================================
  rhs_cat <- paste(c("PRS_q", covs), collapse = " + ")
  fit_cat <- coxph(as.formula(paste("Surv(time, status) ~", rhs_cat)), data = df_q, ties = "efron")
  sm_cat  <- summary(fit_cat)
  
  coef_cat <- sm_cat$coefficients
  ci_cat   <- sm_cat$conf.int
  rows_q   <- grep("^PRS_q", rownames(coef_cat))
  
  HR_categorical <- data.frame(
    outcome = outcome_label,
    quintile = c("Q1", gsub("^PRS_q", "", rownames(coef_cat)[rows_q])),
    HR = c(1, ci_cat[rows_q, "exp(coef)"]),
    CI_low = c(NA, ci_cat[rows_q, "lower .95"]),
    CI_high = c(NA, ci_cat[rows_q, "upper .95"]),
    p_value_quintile = c(NA, coef_cat[rows_q, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  ) %>%
    left_join(count_table, by = c("outcome", "quintile")) %>%
    mutate(quintile = factor(quintile, levels = paste0("Q", 1:5))) %>%
    arrange(outcome, quintile)
  
  # ==========================================================
  # Empirical curves by quintile
  # ==========================================================
  fit_emp <- survfit(Surv(time, status) ~ PRS_q, data = df_q)
  
  # Legend labels based on categorical Cox HRs
  hr_labs <- c(
    "Q1 (HR = 1.00 [ref])",
    paste0("Q2 (HR = ", sprintf("%.2f", HR_categorical$HR[HR_categorical$quintile == "Q2"]), ")"),
    paste0("Q3 (HR = ", sprintf("%.2f", HR_categorical$HR[HR_categorical$quintile == "Q3"]), ")"),
    paste0("Q4 (HR = ", sprintf("%.2f", HR_categorical$HR[HR_categorical$quintile == "Q4"]), ")"),
    paste0("Q5 (HR = ", sprintf("%.2f", HR_categorical$HR[HR_categorical$quintile == "Q5"]), ")")
  )
  
  # Build empirical cumulative incidence plot
  plt <- ggsurvplot(
    fit_emp,
    data = df_q,
    fun = if (y_mode == "event") "event" else "cumhaz",
    conf.int = TRUE,
    conf.int.alpha = 0.08,
    xlim = xlim_years,
    palette = pal_quint,
    legend.title = "PRS risk category",
    legend.labs = paste0("Q", 1:5),
    ggtheme = theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, family = "Arial")
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
  
  plt$plot <- plt$plot +
    scale_y_continuous(
      limits = c(0, ylim_top),
      breaks = seq(0, ylim_top, by = 0.02),
      labels = function(x) paste0(sprintf("%.0f", x * 100), "%")
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
  
  # Extract curve estimate at the end of follow-up
  s16 <- summary(fit_emp, times = xlim_years[2], extend = TRUE)
  
  if (y_mode == "event") {
    est16 <- setNames(1 - s16$surv, nm = paste0("Q", 1:5))
  } else {
    est16 <- setNames(s16$cumhaz, nm = paste0("Q", 1:5))
  }
  
  list(
    plot = plt$plot,
    est16 = est16,
    HR_continuous = HR_continuous,
    HR_trend = HR_trend,
    rate_table = rate_table,
    HR_categorical = HR_categorical
  )
}


## ----------------------------
## Endpoint specification
## ----------------------------
endpoints <- data.frame(
  column_name = c("Bc", "Pc", "Crc", "Chd0"),
  sex_stratified = c("yes", "yes", "no", "no"),
  title = c(
    "Breast cancer (BC)",
    "Prostate cancer (PC)",
    "Colorectal cancer (CRC)",
    "Coronary Heart Disease (CHD)"
  ),
  outcome = c("breast", "prostate", "colorectal", "CHD"),
  acro = c("BC", "PC", "CRC", "CHD"),
  bs_col = c(
    "breast_Bs_final",
    "prostate_Bs_final",
    "colorectal_Bs_final",
    "Chd_Bs_New_final"
  ),
  stringsAsFactors = FALSE
)

# Figure order
endpoints_4 <- endpoints %>%
  mutate(acro = factor(acro, levels = c("CHD", "BC", "PC", "CRC"))) %>%
  arrange(acro)


## ----------------------------
## Run analysis across the four outcomes
## ----------------------------
plot_list_4   <- vector("list", nrow(endpoints_4))
HR_cont_list  <- list()
HR_trend_list <- list()
rates_list    <- list()
HR_cat_list   <- list()

for (i in seq_len(nrow(endpoints_4))) {
  
  pat  <- endpoints_4$column_name[i]
  outc <- endpoints_4$outcome[i]
  ttl  <- endpoints_4$title[i]
  acr  <- as.character(endpoints_4$acro[i])
  bsc  <- endpoints_4$bs_col[i]
  
  # Restrict breast and prostate cancer to the relevant sex
  sex_f <- if (endpoints_4$sex_stratified[i] == "yes") {
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
  
  # Prepare endpoint-specific dataset
  prep <- prep_surv_df(
    df = raw,
    pattern = pat,
    prs_outcome = outc,
    sex_filter = sex_f,
    bs_col = bsc
  )
  
  # Fit models and build empirical curves
  out <- make_curve_plot_empirical(
    prep$df,
    prep$prs_var,
    title_text = ttl,
    outcome_label = acr,
    sex_filter = sex_f,
    y_mode = "event"
  )
  
  plot_list_4[[i]] <- out$plot
  HR_cont_list[[acr]]  <- out$HR_continuous
  HR_trend_list[[acr]] <- out$HR_trend
  rates_list[[acr]]    <- out$rate_table
  HR_cat_list[[acr]]   <- out$HR_categorical
}


## ----------------------------
## Combine summary tables
## ----------------------------
HR_continuous      <- bind_rows(HR_cont_list)
HR_trend           <- bind_rows(HR_trend_list)
rate_table_all     <- bind_rows(rates_list)
HR_categorical_all <- bind_rows(HR_cat_list)

print(HR_categorical_all)


## ----------------------------
## Save summary tables
## ----------------------------
write.csv(
  HR_categorical_all,
  file.path(outdir_tables, "HR_categorical_pvaluePerQuintile_cases_censored_CHD_BC_PC_CRC.csv"),
  row.names = FALSE
)

write.csv(
  HR_continuous,
  file.path(outdir_tables, "HR_continuous_perSD_CHD_BC_PC_CRC.csv"),
  row.names = FALSE
)

write.csv(
  HR_trend,
  file.path(outdir_tables, "HR_trend_ordinalQ1to5_CHD_BC_PC_CRC.csv"),
  row.names = FALSE
)

write.csv(
  rate_table_all,
  file.path(outdir_tables, "RateTable_N_cases_PY_rate_CHD_BC_PC_CRC.csv"),
  row.names = FALSE
)


## ----------------------------
## Assemble Figure 3A
## ----------------------------
final_plot <- ggarrange(plotlist = plot_list_4, ncol = 2, nrow = 2)

final_plot <- annotate_figure(
  final_plot,
  top    = text_grob("A", size = 16, face = "bold", family = "Arial", x = 0.02, hjust = 0),
  bottom = text_grob("Follow-up time (years)", size = 18, face = "bold", family = "Arial"),
  left   = text_grob("Cumulative incidence (%)", size = 18, face = "bold", family = "Arial", rot = 90)
)


## ----------------------------
## Save Figure 3A
## ----------------------------
png(
  file.path(outdir_plots, "Figure3A_SurvivalPlots_CHD_BC_PC_CRC.png"),
  width = 8.3,
  height = 8.3,
  units = "in",
  res = 600
)
print(final_plot)
dev.off()