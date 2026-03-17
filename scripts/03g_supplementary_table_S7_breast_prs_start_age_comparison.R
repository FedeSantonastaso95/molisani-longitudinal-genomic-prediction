############################################################
# 03g_supplementary_table_S7_breast_prs_start_age_comparison.R
#
# Description:
# This script generates a Supplementary Table S7-style summary for
# breast cancer using the same model-based predicted cumulative hazard
# curves shown in Figure 3B.
#
# The table compares:
#   - cases observed by age 50
#   - cases observed by the PRS-based start age
#
# The PRS-based start age is defined as the age at which the predicted
# cumulative hazard curve for each PRS quintile reaches a target value
# (here: 0.010), using the same Cox model and the same predicted curves
# as in Figure 3B.
#
# Important notes:
# - The analysis dataset is restricted to women only.
# - The analysis dataset is restricted to Age_breast <= 80 years.
# - Start ages are derived from model-based predicted curves.
# - Case counts are computed on the observed data.
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(survival)


## ----------------------------
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"
outfile_table <- "path/to/output/Supplementary_Table_S7_breast_prs_start_age_comparison.txt"


## ----------------------------
## 1. Read data
## ----------------------------
df <- fread(
  infile,
  data.table = FALSE
) %>%
  filter(!FID %in% c(
    201629, 201630, 212535, 216032,
    217178, 217489, 400372, 410250
  ))


## ----------------------------
## 2. Recode breast cancer outcome and age at event/censoring
## ----------------------------
#
# Outcome coding:
#   - 2 = breast cancer case
#   - 1 = censored
#
# event_bc is a binary event indicator used for counting cases.
# Age_breast is the age at event or censoring.
df2 <- df %>%
  mutate(
    breast_Bs_final = case_when(
      Mt_Diagnosis_New == 1 & Mt1_Icd9 == 174 ~ 2,
      Fnf_Bc_Fup2020 == 1                     ~ 2,
      TRUE                                    ~ 1
    ),
    event_bc = as.integer(breast_Bs_final == 2),
    Age_breast = case_when(
      breast_Bs_final == 2 &
        Mt_Diagnosis_New == 1 &
        Mt1_Icd9 == 174 ~ Age_Mt1_New,
      
      TRUE ~ (
        as.numeric(
          as.Date(Dataexit_Bc_Fup2020, format = "%Y-%m-%d") -
            as.Date(Recruitment_date, format = "%d/%m/%Y")
        ) / 365.25
      ) + Age
    )
  )


## ----------------------------
## 3. Restrict to women and define standardized PRS
## ----------------------------
pred_column <- "PRED_eurW_breast_noBRCA1_noBRCA2"

df3 <- df2 %>%
  filter(Sex == 0) %>%
  mutate(
    PRS_z = as.numeric(scale(.data[[pred_column]]))
  ) %>%
  filter(
    !is.na(PRS_z),
    !is.na(Age_breast),
    Age_breast <= 80
  )


## ----------------------------
## 4. Define PRS quintiles
## ----------------------------
#
# Quintiles are defined on the standardized PRS distribution.
quint_labels <- c("Q1", "Q2", "Q3", "Q4", "Q5")

quint_breaks <- quantile(
  df3$PRS_z,
  probs = seq(0, 1, by = 0.20),
  na.rm = TRUE
)

df3 <- df3 %>%
  mutate(
    prs_quintile = cut(
      PRS_z,
      breaks = quint_breaks,
      include.lowest = TRUE,
      labels = quint_labels
    ),
    prs_quintile = factor(prs_quintile, levels = quint_labels)
  )


## ----------------------------
## 5. Representative PRS value per quintile
## ----------------------------
#
# The median standardized PRS within each quintile is used to generate
# the predicted curves, exactly as in Figure 3B.
quintile_table <- df3 %>%
  group_by(prs_quintile) %>%
  summarise(
    N = n(),
    PRS_z = median(PRS_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prs_quintile = factor(prs_quintile, levels = quint_labels)
  ) %>%
  arrange(prs_quintile)

print(quintile_table)


## ----------------------------
## 6. Fit the same Cox model used for Figure 3B
## ----------------------------
#
# Age is already the time scale and is therefore not included as a covariate.
# The analysis is restricted to women only, so Sex is not included either.
fit <- coxph(
  Surv(Age_breast, breast_Bs_final) ~ PRS_z + PC1 + PC2 + PC3 + PC4 + PC5,
  data = df3
)


## ----------------------------
## 7. Build the same newdata used for Figure 3B
## ----------------------------
newdata_q <- data.frame(
  PRS_z = quintile_table$PRS_z,
  PC1   = mean(df3$PC1, na.rm = TRUE),
  PC2   = mean(df3$PC2, na.rm = TRUE),
  PC3   = mean(df3$PC3, na.rm = TRUE),
  PC4   = mean(df3$PC4, na.rm = TRUE),
  PC5   = mean(df3$PC5, na.rm = TRUE)
)

sf <- survfit(fit, newdata = newdata_q)


## ----------------------------
## 8. Parameters shared with Figure 3B
## ----------------------------
targets <- c(0.010, 0.013)
target_risk <- 0.010
x_min <- 30


## ----------------------------
## 9. Convert survfit output into a long dataframe
## ----------------------------
curve_ids <- quint_labels
time_vec <- sf$time

surv_mat <- sf$surv
if (is.null(dim(surv_mat))) surv_mat <- matrix(surv_mat, ncol = 1)

cumhaz_mat <- sf$cumhaz
if (is.null(dim(cumhaz_mat))) cumhaz_mat <- matrix(cumhaz_mat, ncol = 1)

lower_mat <- sf$lower
if (!is.null(lower_mat) && is.null(dim(lower_mat))) {
  lower_mat <- matrix(lower_mat, ncol = 1)
}

upper_mat <- sf$upper
if (!is.null(upper_mat) && is.null(dim(upper_mat))) {
  upper_mat <- matrix(upper_mat, ncol = 1)
}

std_err_mat <- sf$std.err
if (!is.null(std_err_mat) && is.null(dim(std_err_mat))) {
  std_err_mat <- matrix(std_err_mat, ncol = 1)
}

n_curves <- ncol(surv_mat)

if (n_curves != length(curve_ids)) {
  stop("Number of curves in survfit does not match the expected number of quintiles.")
}

df_list <- vector("list", n_curves)

for (j in seq_len(n_curves)) {
  tmp <- data.frame(
    time   = time_vec,
    strata = curve_ids[j],
    surv   = surv_mat[, j],
    cumhaz = cumhaz_mat[, j]
  )
  
  if (!is.null(lower_mat) && !is.null(upper_mat)) {
    tmp$lower_surv <- pmin(pmax(lower_mat[, j], 1e-10), 1)
    tmp$upper_surv <- pmin(pmax(upper_mat[, j], 1e-10), 1)
    tmp$lower <- pmax(0, -log(tmp$upper_surv))
    tmp$upper <- pmax(0, -log(tmp$lower_surv))
  } else if (!is.null(std_err_mat)) {
    z <- qnorm(.975)
    tmp$se_surv <- std_err_mat[, j]
    tmp$se_ch <- ifelse(tmp$surv > 0, tmp$se_surv / tmp$surv, NA_real_)
    tmp$lower <- pmax(0, tmp$cumhaz - z * tmp$se_ch)
    tmp$upper <- tmp$cumhaz + z * tmp$se_ch
  } else {
    tmp$lower <- NA_real_
    tmp$upper <- NA_real_
  }
  
  df_list[[j]] <- tmp
}

df0 <- bind_rows(df_list) %>%
  mutate(
    strata = factor(strata, levels = quint_labels)
  )


## ----------------------------
## 10. Re-zero cumulative hazard at age 30
## ----------------------------
#
# This follows the same logic used in Figure 3B.
df1 <- df0 %>%
  group_by(strata) %>%
  mutate(
    ch_at_xmin = {
      v <- cumhaz[time <= x_min]
      if (length(v) == 0) 0 else max(v, na.rm = TRUE)
    },
    cumhaz = cumhaz - ch_at_xmin,
    lower  = pmax(0, lower - ch_at_xmin),
    upper  = pmax(0, upper - ch_at_xmin)
  ) %>%
  ungroup()


## ----------------------------
## 11. Function to extract the start age from predicted curves
## ----------------------------
#
# The start age is defined as the first age at which the re-zeroed
# cumulative hazard reaches or exceeds the target value.
get_cross_target <- function(q_label, target_value) {
  d <- df1 %>%
    filter(strata == q_label) %>%
    arrange(time)
  
  if (nrow(d) == 0) return(NA_real_)
  if (max(d$cumhaz, na.rm = TRUE) < target_value) return(NA_real_)
  
  d$time[which(d$cumhaz >= target_value)[1]]
}


## ----------------------------
## 12. Start ages derived from the same predicted curves
## ----------------------------
start_age_df <- data.frame(
  prs_quintile = quint_labels,
  start_age_exact = sapply(quint_labels, get_cross_target, target_value = target_risk)
) %>%
  mutate(
    prs_quintile = factor(prs_quintile, levels = quint_labels),
    start_age_display = round(start_age_exact)
  ) %>%
  arrange(prs_quintile)

print(start_age_df)


## ----------------------------
## 13. Optional: crossings for all target values
## ----------------------------
cross_df_all <- expand.grid(
  strata = quint_labels,
  target = targets,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    age = get_cross_target(strata, target)
  ) %>%
  ungroup()

print(cross_df_all)


## ----------------------------
## 14. Build Supplementary Table S7-style summary
## ----------------------------
#
# Definitions:
#   N                 = all women in each PRS quintile
#   Cases <=50        = observed breast cancer cases by age 50
#   Cases <=PRS start = observed breast cancer cases by the PRS-based start age
#
# Additional columns:
#   Cases 50-PRS start
#     = extra cases covered if the PRS-based start age is later than 50
#
#   Cases <50 captured
#     = cases between PRS-based start age and 50 if the PRS-based start age
#       is earlier than 50
table_quintiles <- df3 %>%
  left_join(start_age_df, by = "prs_quintile") %>%
  group_by(prs_quintile, start_age_exact, start_age_display) %>%
  summarise(
    N = n(),
    
    `Cases <=50` = sum(
      event_bc == 1 & Age_breast <= 50,
      na.rm = TRUE
    ),
    
    `Cases <=PRS start` = sum(
      event_bc == 1 & Age_breast <= start_age_exact[1],
      na.rm = TRUE
    ),
    
    `Cases 50-PRS start` = ifelse(
      is.na(start_age_exact[1]),
      NA_integer_,
      ifelse(
        start_age_exact[1] > 50,
        sum(
          event_bc == 1 &
            Age_breast > 50 &
            Age_breast <= start_age_exact[1],
          na.rm = TRUE
        ),
        0L
      )
    ),
    
    `Cases <50 captured` = ifelse(
      is.na(start_age_exact[1]),
      NA_integer_,
      ifelse(
        start_age_exact[1] < 50,
        sum(
          event_bc == 1 &
            Age_breast > start_age_exact[1] &
            Age_breast <= 50,
          na.rm = TRUE
        ),
        0L
      )
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    `Start age (PRS)` = start_age_display
  ) %>%
  select(
    `PRS quintile` = prs_quintile,
    N,
    `Start age (PRS)`,
    `Cases <=50`,
    `Cases <=PRS start`,
    `Cases 50-PRS start`,
    `Cases <50 captured`
  )

print(table_quintiles)


## ----------------------------
## 15. Optional pretty-print version
## ----------------------------
table_quintiles_pretty <- table_quintiles %>%
  mutate(
    N = format(N, big.mark = ",", scientific = FALSE)
  )

print(table_quintiles_pretty)


## ----------------------------
## 16. Simple checks
## ----------------------------
cat("Total women in table:", sum(table_quintiles$N), "\n")

start_age_df_check <- start_age_df %>%
  mutate(
    start_age_display = round(start_age_exact)
  ) %>%
  arrange(prs_quintile)

print(start_age_df_check)


## ----------------------------
## 17. Save table
## ----------------------------
fwrite(
  table_quintiles,
  outfile_table,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)