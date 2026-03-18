############################################################
# Breast cancer cumulative risk at ages 50 and 52
#
# This script estimates the cumulative breast cancer rate from a Cox model
# using age as the time scale and a standardized continuous PRS
# (breast PRS excluding BRCA1/2).
#
# The cumulative rate is extracted at:
#   - age 50 years
#   - age 52 years
#
# The age-52 estimate was added to reflect the next biennial screening
# interval after age 50.
#
# Notes:
# - The analysis is restricted to women only.
# - Follow-up is restricted to age <= 80 years.
# - The plotted curve is model-based and corresponds to the average covariate
#   profile in the study population, with PRS fixed at its mean (z = 0).
############################################################


## ----------------------------
## Libraries
## ----------------------------
conflicted::conflicts_prefer(dplyr::filter)

library(data.table)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(scales)


## ----------------------------
## 1. Read input data
## ----------------------------
df <- fread(
  "path/to/phenotype_dataset.txt",
  data.table = FALSE
) 


## ----------------------------
## 2. Recode breast cancer outcome
## ----------------------------
#
# Outcome coding for Surv():
#   - 2 = breast cancer case
#   - 1 = censored
#
# Age_breast is the age at event or censoring:
#   - prevalent cases: reported age at onset
#   - incident cases and censored individuals: baseline age plus follow-up time
df2 <- df %>%
  mutate(
    breast_Bs_final = ifelse(
      Mt_Diagnosis_New == 1 & Mt1_Icd9 == 174,
      2,
      ifelse(Fnf_Bc_Fup2020 == 1, 2, 1)
    ),
    Age_breast = case_when(
      breast_Bs_final == 2 &
        Mt_Diagnosis_New == 1 &
        Mt1_Icd9 == 174 ~ Age_Mt1_New,
      
      TRUE ~ (
        as.numeric(
          as.Date(Dataexit_Bc_Fup2020, "%Y-%m-%d") -
            as.Date(Recruitment_date, "%d/%m/%Y")
        ) / 365.25
      ) + Age
    )
  )


## ----------------------------
## 3. Restrict to women and define PRS
## ----------------------------
df3 <- df2 %>%
  filter(Sex == 0) %>%
  mutate(
    PRS_z = as.numeric(scale(PRED_eurW_breast_noBRCA1_noBRCA2))
  ) %>%
  filter(Age_breast <= 80)


## ----------------------------
## 4. Fit Cox proportional hazards model
## ----------------------------
#
# The model is adjusted for:
#   - Age
#   - Sex
#   - PC1-PC5
#   - standardized PRS
#
# Note:
# Sex is constant in this women-only dataset, so it is not necessary.
# It is therefore omitted from the rewritten model below.
fit <- coxph(
  Surv(Age_breast, breast_Bs_final) ~
    Age + PC1 + PC2 + PC3 + PC4 + PC5 + PRS_z,
  data = df3
)


## ----------------------------
## 5. Generate predicted cumulative hazard curve
## ----------------------------
#
# The predicted curve is obtained for an "average" participant:
#   - Age set to the sample mean
#   - PC1-PC5 set to their sample means
#   - PRS_z fixed at 0 (mean PRS)
nd <- data.frame(
  Age   = mean(df3$Age, na.rm = TRUE),
  PC1   = mean(df3$PC1, na.rm = TRUE),
  PC2   = mean(df3$PC2, na.rm = TRUE),
  PC3   = mean(df3$PC3, na.rm = TRUE),
  PC4   = mean(df3$PC4, na.rm = TRUE),
  PC5   = mean(df3$PC5, na.rm = TRUE),
  PRS_z = 0
)

sf <- survfit(fit, newdata = nd)


## ----------------------------
## 6. Extract cumulative hazard at ages 50 and 52
## ----------------------------
#
# For each target age, the closest available time point in the survfit object
# is used.
get_val <- function(target_age, survfit_obj) {
  idx <- which.min(abs(survfit_obj$time - target_age))
  list(
    time = survfit_obj$time[idx],
    cumhaz = survfit_obj$cumhaz[idx]
  )
}

v50 <- get_val(50, sf)
v52 <- get_val(52, sf)

lab_50 <- sprintf("%.3f", v50$cumhaz)
lab_52 <- sprintf("%.3f", v52$cumhaz)


## ----------------------------
## 7. Base plot without legend
## ----------------------------
p0 <- ggsurvplot(
  sf,
  data = df3,
  fun = "cumhaz",
  conf.int = TRUE,
  conf.int.alpha = 0.08,
  censor = FALSE,
  legend = "none",
  palette = "#cb271b",
  xlab = "Age (years)",
  ylab = "Cumulative breast cancer rate",
  ylim = c(0, 0.05),
  ggtheme = theme_bw()
)$plot


## ----------------------------
## 8. Final formatted plot
## ----------------------------
p <- p0 +
  coord_cartesian(
    xlim = c(0, 80),
    ylim = c(0, 0.052),
    clip = "off"
  ) +
  scale_x_continuous(
    breaks = seq(0, 80, 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.05, 0.01),
    labels = function(x) sub("\\.?0+$", "", sprintf("%.2f", x)),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  
  # Vertical guide lines at ages 50 and 52
  geom_vline(xintercept = 50, linetype = "dotted", color = "#bcbcbc") +
  geom_vline(xintercept = 52, linetype = "dotted", color = "black") +
  
  # Points marking the cumulative hazard at ages 50 and 52
  geom_point(aes(x = v50$time, y = v50$cumhaz), color = "#cb271b", size = 3) +
  geom_point(aes(x = v52$time, y = v52$cumhaz), color = "#cb271b", size = 3) +
  
  # Horizontal guide lines from y-axis to ages 50 and 52
  geom_segment(
    aes(x = 0, xend = 50, y = v50$cumhaz, yend = v50$cumhaz),
    linetype = "dotted",
    color = "#bcbcbc"
  ) +
  geom_segment(
    aes(x = 0, xend = 52, y = v52$cumhaz, yend = v52$cumhaz),
    linetype = "dotted",
    color = "black"
  ) +
  
  # Labels for the cumulative hazard values
  annotate(
    "text",
    x = 0.8,
    y = v50$cumhaz,
    label = lab_50,
    hjust = 0,
    vjust = -0.35,
    color = "#7f7f7f",
    size = 4
  ) +
  annotate(
    "text",
    x = 0.8,
    y = v52$cumhaz,
    label = lab_52,
    hjust = 0,
    vjust = -0.35,
    color = "black",
    size = 4
  )

print(p)


## ----------------------------
## 9. Optional export
## ----------------------------
png(
  "path/to/output/CumulativeRate_By50_By52_AllWomen_by80.png",
  width = 6,
  height = 4.5,
  units = "in",
  res = 600
)
print(p)
dev.off()