############################################################
# Figure 3B
# Breast cancer cumulative incidence by PRS quintiles
#
# This script:
#   1. defines a breast cancer time-to-event dataset in women only
#   2. fits a Cox model with continuous standardized PRS
#   3. generates predicted cumulative incidence curves at representative
#      PRS values (median PRS_z within each quintile)
#   4. identifies the ages at which each curve crosses two predefined
#      cumulative incidence thresholds
#   5. adds points at those intersections
#   6. prints the crossing ages in the console
#
# Notes:
# - Age is used as the time scale.
# - The analysis is restricted to age <= 80 years.
# - PRS quintiles are defined from the standardized PRS distribution.
# - Curves are plotted as cumulative incidence:
#     1 - S(t)
############################################################


## ----------------------------
## Libraries
## ----------------------------
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(grid)
library(scales)


## ----------------------------
## Input / output
## ----------------------------
infile <- "path/to/phenotype_dataset.txt"
outfile_png <- "path/to/output/Figure3B_BC_quintiles.png"


## ----------------------------
## Read data
## ----------------------------
df <- fread(
  infile,
  data.table = FALSE
)


## ----------------------------
## Recode breast cancer outcome
## ----------------------------
#
# Outcome coding for Surv():
#   - 2 = breast cancer case
#   - 1 = censored
#
# Age_breast is the age at event or censoring:
#   - prevalent cases: reported age at onset
#   - all other women: baseline age + observed follow-up time
df2 <- df %>%
  mutate(
    breast_Bs_final = case_when(
      Mt_Diagnosis_New == 1 & Mt1_Icd9 == 174 ~ 2,
      Fnf_Bc_Fup2020 == 1                     ~ 2,
      TRUE                                    ~ 1
    ),
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
## Women only + standardized PRS
## ----------------------------
pred_column <- "PRED_eurW_breast_noBRCA1_noBRCA2"

df3 <- df2 %>%
  filter(Sex == 0) %>%
  mutate(
    PRS_z = as.numeric(scale(.data[[pred_column]]))
  ) %>%
  filter(!is.na(PRS_z), !is.na(Age_breast), Age_breast <= 80)


## ----------------------------
## Define empirical PRS quintiles
## ----------------------------
#
# ntile() is used here to avoid issues with duplicated quantile cut points.
quint_labels <- c("Q1", "Q2", "Q3", "Q4", "Q5")

df3 <- df3 %>%
  mutate(
    prs_quintile_num = ntile(PRS_z, 5),
    prs_quintile = factor(
      paste0("Q", prs_quintile_num),
      levels = quint_labels
    )
  )


## ----------------------------
## Representative PRS value per quintile
## ----------------------------
#
# Each predicted curve is generated using the median PRS_z observed
# within the corresponding empirical quintile.
quintile_table <- df3 %>%
  group_by(prs_quintile) %>%
  summarise(
    PRS_z = median(PRS_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prs_quintile = factor(prs_quintile, levels = quint_labels)
  ) %>%
  arrange(prs_quintile)

print(quintile_table)


## ----------------------------
## Cox model with continuous PRS
## ----------------------------
#
# Age is already the time scale, so it is not included as a covariate.
# Sex is omitted because the analysis is restricted to women only.
fit <- coxph(
  Surv(Age_breast, breast_Bs_final) ~ PRS_z + PC1 + PC2 + PC3 + PC4 + PC5,
  data = df3
)


## ----------------------------
## Newdata for predicted curves
## ----------------------------
#
# Curves are predicted at:
#   - PRS_z = median PRS_z within each quintile
#   - PC1-PC5 = sample means
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
## Plot parameters
## ----------------------------
targets <- c(0.010, 0.013)

x_min <- 30
x_max_base <- 70
y_min <- 0
y_max <- 0.05

cap_to_decade <- function(x) ceiling(x / 10) * 10

pal <- c(
  "Q1" = "#2C7FB8",
  "Q2" = "#9ECAE1",
  "Q3" = "#969696",
  "Q4" = "#F4A6A6",
  "Q5" = "#D7301F"
)

hline_cols <- c(
  "0.01"  = "black",
  "0.013" = "#238B45"
)

y_lab_fun <- function(x) {
  paste0(sprintf("%.0f", x * 100), "%")
}


## ----------------------------
## Convert survfit object to plotting dataframe
## ----------------------------
curve_ids <- quint_labels
time_vec <- sf$time

surv_mat <- sf$surv
if (is.null(dim(surv_mat))) surv_mat <- matrix(surv_mat, ncol = 1)

lower_mat <- sf$lower
if (!is.null(lower_mat) && is.null(dim(lower_mat))) {
  lower_mat <- matrix(lower_mat, ncol = 1)
}

upper_mat <- sf$upper
if (!is.null(upper_mat) && is.null(dim(upper_mat))) {
  upper_mat <- matrix(upper_mat, ncol = 1)
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
    surv   = surv_mat[, j]
  )
  
  tmp$cuminc <- 1 - tmp$surv
  
  # Confidence intervals transformed from survival to cumulative incidence
  if (!is.null(lower_mat) && !is.null(upper_mat)) {
    tmp$lower_surv <- pmin(pmax(lower_mat[, j], 1e-10), 1)
    tmp$upper_surv <- pmin(pmax(upper_mat[, j], 1e-10), 1)
    tmp$lower <- pmax(0, 1 - tmp$upper_surv)
    tmp$upper <- pmin(1, 1 - tmp$lower_surv)
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
## Determine age at crossing for each target
## ----------------------------
#
# For each quintile-specific predicted curve, find the first age at which
# cumulative incidence reaches or exceeds the requested threshold.
get_cross_target <- function(q_label, target_value, data_curves) {
  d <- data_curves %>%
    filter(strata == q_label) %>%
    arrange(time)
  
  if (nrow(d) == 0) return(NA_real_)
  if (max(d$cuminc, na.rm = TRUE) < target_value) return(NA_real_)
  
  d$time[which(d$cuminc >= target_value)[1]]
}

cross_df_all <- expand.grid(
  strata = quint_labels,
  target = targets,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    age = get_cross_target(strata, target, df0)
  ) %>%
  ungroup()

print(cross_df_all)


## ----------------------------
## Print crossing ages in a clean table
## ----------------------------
#
# This table reports, for each PRS quintile, the age at which the
# predicted cumulative incidence curve first reaches each target.
print_crossing_ages <- function(crossing_df) {
  cross_wide <- crossing_df %>%
    mutate(
      target_label = case_when(
        target == 0.010 ~ "age_at_1.0pct",
        target == 0.013 ~ "age_at_1.3pct",
        TRUE ~ paste0("age_at_", target)
      )
    ) %>%
    select(strata, target_label, age) %>%
    tidyr::pivot_wider(
      names_from = target_label,
      values_from = age
    ) %>%
    arrange(strata)
  
  cat("\n============================================================\n")
  cat("Ages at which each PRS quintile reaches the target cumulative incidence\n")
  cat("============================================================\n")
  print(cross_wide)
  cat("============================================================\n\n")
  
  invisible(cross_wide)
}

crossing_age_table <- print_crossing_ages(cross_df_all)


## ----------------------------
## Dynamic x-axis upper limit
## ----------------------------
cross_max <- suppressWarnings(max(cross_df_all$age, na.rm = TRUE))
if (!is.finite(cross_max)) cross_max <- x_max_base
x_max <- if (cross_max > x_max_base) cap_to_decade(cross_max) else x_max_base


## ----------------------------
## Re-zero curves at x_min
## ----------------------------
#
# Curves are displayed from age 30 onward and re-centered so that
# cumulative incidence starts at 0 at age 30.
df1 <- df0 %>%
  group_by(strata) %>%
  mutate(
    ci_at_xmin = {
      v <- cuminc[time <= x_min]
      if (length(v) == 0) 0 else max(v, na.rm = TRUE)
    },
    cuminc = cuminc - ci_at_xmin,
    lower  = pmax(0, lower - ci_at_xmin),
    upper  = pmax(0, upper - ci_at_xmin)
  ) %>%
  ungroup()


## ----------------------------
## Restrict plotting range and force point at x_min
## ----------------------------
df_plot <- df1 %>%
  filter(time <= x_max) %>%
  group_by(strata) %>%
  arrange(time, .by_group = TRUE) %>%
  {
    add0 <- distinct(., strata) %>%
      mutate(
        time = x_min,
        cuminc = 0,
        lower = 0,
        upper = 0,
        surv = NA_real_,
        ci_at_xmin = NA_real_
      )
    bind_rows(add0, .)
  } %>%
  filter(time >= x_min, time <= x_max) %>%
  arrange(strata, time) %>%
  mutate(
    cuminc = ifelse(time == x_min, 0, cuminc),
    lower  = ifelse(time == x_min, 0, lower),
    upper  = ifelse(time == x_min, 0, upper)
  ) %>%
  ungroup()

if (nrow(df_plot) == 0) {
  stop("df_plot is empty after range restriction. Please check x_min and x_max.")
}


## ----------------------------
## Plotting data for intersections
## ----------------------------
df_plot <- df_plot %>%
  mutate(
    strata_lab = factor(strata, levels = quint_labels, labels = quint_labels)
  )

cross_df <- cross_df_all %>%
  filter(!is.na(age), age >= x_min, age <= x_max) %>%
  mutate(
    strata = factor(strata, levels = quint_labels),
    strata_lab = factor(strata, levels = quint_labels, labels = quint_labels),
    point_y = target
  )

print(cross_df)


## ----------------------------
## Horizontal lines data
## ----------------------------
hline_df <- data.frame(
  target = targets,
  col = c("black", "#238B45")
)


## ----------------------------
## Plot
## ----------------------------
p <- ggplot(df_plot, aes(x = time, y = cuminc)) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = strata_lab),
    alpha = 0.045,
    colour = NA,
    show.legend = FALSE
  ) +
  
  geom_step(
    aes(color = strata_lab, linetype = strata_lab),
    linewidth = 0.95,
    show.legend = TRUE
  ) +
  
  geom_segment(
    data = hline_df,
    aes(x = x_min, xend = x_max, y = target, yend = target),
    inherit.aes = FALSE,
    linetype = "dotted",
    color = hline_df$col,
    linewidth = 0.8
  ) +
  
  geom_point(
    data = cross_df,
    aes(x = age, y = point_y, color = strata_lab),
    size = 2.8,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = pal, name = "PRS quintile") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_linetype_manual(values = setNames(rep("solid", 5), quint_labels), name = "PRS quintile") +
  
  scale_x_continuous(
    breaks = seq(x_min, x_max, by = 10),
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  scale_y_continuous(
    breaks = seq(y_min, y_max, by = 0.01),
    labels = y_lab_fun
  ) +
  
  coord_cartesian(
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.1)),
    linetype = "none"
  ) +
  
  labs(
    x = "Age (years)",
    y = "Breast cancer incidence (%)"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 14),
    legend.position = c(0.02, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(p)


## ----------------------------
## Optional export
## ----------------------------
png(
  outfile_png,
  width = 6.5,
  height = 3.8,
  units = "in",
  res = 600
)
print(p)
dev.off()