# =========================================================
# In this script, you will find the final Figure 5 with two panels
# and the final supplementary table with the four requested models,
# all adjusted for Age, Sex, and PC1-PC5.
#
# IMPORTANT:
# - In the PRS-only model, PRS quintiles are defined on all participants.
# - In the variant-specific models, PRS quintiles are defined ONLY among
#   non-carriers of the tested rare variant; carriers are treated as a
#   separate group.
#
# Final 4-panel figure:
# - A: PRS only
# - B: LPA (rs3798220_C)
# - C: LDLR (rs730882080_T)
# - D: SLC4A11 (rs58757394_T)
#
# Final supplementary table:
# - PRS only                -> reference = Q1
# - LPA rs3798220           -> reference = Q1
# - LDLR rs730882080        -> reference = Q1
# - SLC4A11 rs58757394      -> reference = Q5
# =========================================================

# =========================================================
# Libraries
# =========================================================
library(survival)
library(survminer)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(rlang)
library(data.table)
library(ggplot2)
library(patchwork)
library(ggtext)
library(cowplot)
library(grid)

theme_set(theme_bw(base_family = "Arial"))

# =========================================================
# Covariate names: EDIT HERE IF NEEDED
# =========================================================
age_var <- "Age"
sex_var <- "Sex"
pc_vars <- c("PC1", "PC2", "PC3", "PC4", "PC5")

# =========================================================
# Load carrier file
# =========================================================
df_carr <- fread("/path/to/carriers_rarVar.raw")

df_carr[, 7:ncol(df_carr)] <- lapply(df_carr[, 7:ncol(df_carr)], round)

# PLINK dosages correspond to the major allele in this dataset.
# Recode so that:
#   0 = non-carrier of the rare allele
#   1 = heterozygous carrier
#   2 = homozygous carrier of the rare allele
df_carr[, 7:ncol(df_carr)] <- lapply(df_carr[, 7:ncol(df_carr)], function(col) {
  if (is.numeric(col)) {
    col_new <- col
    col_new[col == 0] <- 2
    col_new[col == 2] <- 0
    return(col_new)
  } else {
    return(col)
  }
})

old_names <- names(df_carr)[7:ncol(df_carr)]
new_names <- sub("_(?!.*_).*$", "", old_names, perl = TRUE)
names(df_carr)[7:ncol(df_carr)] <- new_names

# =========================================================
# Load rare variants and SNP list
# =========================================================
var <- fread("/path/to/listVariants_Burden_ALL_final.txt")

snplist <- fread("/path/to/imputed_genotypes.snplist")

var1 <- var %>%
  left_join(
    snplist,
    by = c(
      "GENPOS"  = "V4",
      "CHROM"   = "V1",
      "ALLELE0" = "V5",
      "ALLELE1" = "V6"
    )
  ) %>%
  dplyr::rename(bgen_ID = V2) %>%
  dplyr::select(-V3)

setDT(var1)
unique_combos <- unique(var1[, .(gene_name, trait)])

# =========================================================
# Load sample file
# =========================================================
sample <- fread("/path/to/imputed_genotypes.sample")

# Remove the first line of the PLINK/BGEN sample file if it contains the
# standard dummy row
sample <- sample[-1, ]

# =========================================================
# Load phenotype data
# =========================================================
df <- fread(
  "/path/to/phenotypes_covariates_prs.txt",
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
  mutate(Fnf_cancer_FUP2020 = ifelse(is.na(Fnf_cancer_FUP2020), 0, Fnf_cancer_FUP2020)) %>%
  left_join(df_carr, by = c("FID" = "FID", "IID" = "IID"))

df1 <- sample %>%
  left_join(df, by = c("ID_1" = "FID"))

df1 <- df1 %>%
  select(
    IID,
    Chd0_Fnof_Fup2020,
    Annipers_Chd0_Fnof_Fup2020,
    rs3798220,
    rs58757394,
    rs730882080,
    PRED_eurW_CHD,
    Bio_Lp_A,
    SBP,
    LDL,
    Chd_Bs_New,
    all_of(age_var),
    all_of(sex_var),
    all_of(pc_vars)
  ) %>%
  rename(
    Age = all_of(age_var),
    Sex = all_of(sex_var),
    PC1 = all_of(pc_vars[1]),
    PC2 = all_of(pc_vars[2]),
    PC3 = all_of(pc_vars[3]),
    PC4 = all_of(pc_vars[4]),
    PC5 = all_of(pc_vars[5])
  ) %>%
  mutate(Chd_Bs_New_final = ifelse(Chd_Bs_New == 0, 0, 1)) %>%
  mutate(
    Chd0_Fnof_Fup2020 = if_else(
      Chd_Bs_New_final == 1 & Chd0_Fnof_Fup2020 == 1,
      as.numeric(NA),
      as.numeric(Chd0_Fnof_Fup2020)
    )
  ) %>%
  filter(Chd_Bs_New_final != 1)

# =========================================================
# Helpers
# =========================================================
y_label_formatter <- function(x) {
  ifelse(x == 0, "0", sprintf("%.2f", x))
}

mm_to_in <- function(mm) mm / 25.4

x_lab_common <- "Time (years)"
y_lab_common <- "Cumulative CHD hazard"

q_labels <- paste0("Q", 1:5)
legend_levels <- c(q_labels, "Carrier")

global_y_max <- 0.08
y_breaks_common <- seq(0, global_y_max, by = 0.02)

base_panel_theme <- function() {
  theme_classic(base_family = "Arial") +
    theme(
      axis.title = element_text(size = 11, family = "Arial", face = "bold"),
      axis.text = element_text(size = 10, family = "Arial"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.4),
      panel.grid = element_blank(),
      plot.margin = margin(12, 8, 8, 8),
      plot.title = ggtext::element_markdown(
        family = "Arial",
        face = "bold",
        hjust = 0.5,
        size = 11.5
      )
    )
}

get_mode_value <- function(x) {
  x_nonmiss <- x[!is.na(x)]
  ux <- unique(x_nonmiss)
  ux[which.max(tabulate(match(x_nonmiss, ux)))]
}

build_reference_covariates <- function(df_model) {
  data.frame(
    Age = median(df_model$Age, na.rm = TRUE),
    Sex = get_mode_value(df_model$Sex),
    PC1 = median(df_model$PC1, na.rm = TRUE),
    PC2 = median(df_model$PC2, na.rm = TRUE),
    PC3 = median(df_model$PC3, na.rm = TRUE),
    PC4 = median(df_model$PC4, na.rm = TRUE),
    PC5 = median(df_model$PC5, na.rm = TRUE)
  )
}

get_prs_quintile_medians <- function(df, prs_var = "PRED_eurW_CHD") {
  df %>%
    mutate(
      CHD_quint = ntile(.data[[prs_var]], 5),
      CHD_quint = factor(CHD_quint, levels = 1:5, labels = paste0("Q", 1:5))
    ) %>%
    group_by(CHD_quint) %>%
    summarise(
      prs_value = median(.data[[prs_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(CHD_quint)
}

get_noncarrier_prs_quintile_medians <- function(df, prs_var = "PRED_eurW_CHD") {
  df %>%
    mutate(
      CHD_quint = ntile(.data[[prs_var]], 5),
      CHD_quint = factor(CHD_quint, levels = 1:5, labels = paste0("Q", 1:5))
    ) %>%
    group_by(CHD_quint) %>%
    summarise(
      prs_value = median(.data[[prs_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(CHD_quint)
}

# =========================================================
# Extract predicted cumulative hazard + CI from coxph
# =========================================================
get_predicted_curves <- function(fit, newdata, legend_var = "legend_group") {
  
  sf <- survfit(fit, newdata = newdata, conf.type = "log-log")
  n_curves <- nrow(newdata)
  
  if (n_curves == 1) {
    surv_mat  <- matrix(sf$surv,  ncol = 1)
    lower_mat <- matrix(sf$lower, ncol = 1)
    upper_mat <- matrix(sf$upper, ncol = 1)
  } else {
    surv_mat  <- sf$surv
    lower_mat <- sf$lower
    upper_mat <- sf$upper
  }
  
  out <- lapply(seq_len(n_curves), function(i) {
    data.frame(
      time = sf$time,
      surv = surv_mat[, i],
      lower_surv = lower_mat[, i],
      upper_surv = upper_mat[, i],
      legend_group = as.character(newdata[[legend_var]][i])
    )
  }) %>%
    bind_rows() %>%
    mutate(
      legend_group = factor(legend_group, levels = legend_levels),
      cumhaz = -log(pmax(surv, 1e-12)),
      cumhaz_low = -log(pmax(upper_surv, 1e-12)),
      cumhaz_high = -log(pmax(lower_surv, 1e-12))
    )
  
  out0 <- newdata %>%
    transmute(
      time = 0,
      surv = 1,
      lower_surv = 1,
      upper_surv = 1,
      legend_group = factor(.data[[legend_var]], levels = legend_levels),
      cumhaz = 0,
      cumhaz_low = 0,
      cumhaz_high = 0
    )
  
  bind_rows(out0, out) %>%
    arrange(legend_group, time)
}

# =========================================================
# Build PRS-only prediction dataframe
# =========================================================
build_prs_only_sf <- function(df) {
  
  df_prs <- df %>%
    filter(
      !is.na(PRED_eurW_CHD),
      !is.na(Annipers_Chd0_Fnof_Fup2020),
      !is.na(Chd0_Fnof_Fup2020),
      !is.na(Age),
      !is.na(Sex),
      !is.na(PC1),
      !is.na(PC2),
      !is.na(PC3),
      !is.na(PC4),
      !is.na(PC5)
    )
  
  fit <- coxph(
    Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
      PRED_eurW_CHD + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
    data = df_prs,
    ties = "efron",
    x = TRUE
  )
  
  ref_cov <- build_reference_covariates(df_prs)
  prs_medians <- get_prs_quintile_medians(df_prs, prs_var = "PRED_eurW_CHD")
  
  newdata <- prs_medians %>%
    mutate(
      Age = ref_cov$Age[1],
      Sex = ref_cov$Sex[1],
      PC1 = ref_cov$PC1[1],
      PC2 = ref_cov$PC2[1],
      PC3 = ref_cov$PC3[1],
      PC4 = ref_cov$PC4[1],
      PC5 = ref_cov$PC5[1],
      legend_group = factor(as.character(CHD_quint), levels = legend_levels)
    ) %>%
    transmute(
      PRED_eurW_CHD = prs_value,
      Age, Sex, PC1, PC2, PC3, PC4, PC5,
      legend_group
    )
  
  get_predicted_curves(fit, newdata = newdata, legend_var = "legend_group")
}

# =========================================================
# Build PRS + variant prediction dataframe
# =========================================================
build_prs_variant_sf <- function(df, variant_col) {
  
  var_sym <- rlang::ensym(variant_col)
  
  df_all <- df %>%
    mutate(
      carrier = ifelse(!!var_sym == 1, 1, 0)
    ) %>%
    filter(
      !is.na(PRED_eurW_CHD),
      !is.na(Annipers_Chd0_Fnof_Fup2020),
      !is.na(Chd0_Fnof_Fup2020),
      !is.na(carrier),
      !is.na(Age),
      !is.na(Sex),
      !is.na(PC1),
      !is.na(PC2),
      !is.na(PC3),
      !is.na(PC4),
      !is.na(PC5)
    )
  
  fit <- coxph(
    Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
      PRED_eurW_CHD + carrier + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
    data = df_all,
    ties = "efron",
    x = TRUE
  )
  
  ref_cov <- build_reference_covariates(df_all)
  
  prs_medians_noncarrier <- df_all %>%
    filter(carrier == 0) %>%
    get_noncarrier_prs_quintile_medians(prs_var = "PRED_eurW_CHD") %>%
    mutate(
      carrier = 0,
      Age = ref_cov$Age[1],
      Sex = ref_cov$Sex[1],
      PC1 = ref_cov$PC1[1],
      PC2 = ref_cov$PC2[1],
      PC3 = ref_cov$PC3[1],
      PC4 = ref_cov$PC4[1],
      PC5 = ref_cov$PC5[1],
      legend_group = factor(as.character(CHD_quint), levels = legend_levels)
    ) %>%
    select(PRED_eurW_CHD = prs_value, carrier, Age, Sex, PC1, PC2, PC3, PC4, PC5, legend_group)
  
  carrier_prs_value <- df_all %>%
    filter(carrier == 1) %>%
    summarise(prs_value = median(PRED_eurW_CHD, na.rm = TRUE)) %>%
    pull(prs_value)
  
  if (length(carrier_prs_value) == 0 || is.na(carrier_prs_value)) {
    carrier_prs_value <- median(df_all$PRED_eurW_CHD, na.rm = TRUE)
  }
  
  carrier_row <- data.frame(
    PRED_eurW_CHD = carrier_prs_value,
    carrier = 1,
    Age = ref_cov$Age[1],
    Sex = ref_cov$Sex[1],
    PC1 = ref_cov$PC1[1],
    PC2 = ref_cov$PC2[1],
    PC3 = ref_cov$PC3[1],
    PC4 = ref_cov$PC4[1],
    PC5 = ref_cov$PC5[1],
    legend_group = factor("Carrier", levels = legend_levels)
  )
  
  newdata <- bind_rows(prs_medians_noncarrier, carrier_row)
  
  get_predicted_curves(fit, newdata = newdata, legend_var = "legend_group")
}

# =========================================================
# Plotting function
# =========================================================
make_panel_plot <- function(sf,
                            plot_title,
                            y_max = 0.12,
                            show_carrier_ci = TRUE,
                            carrier_ci_alpha = 0.04,
                            legend_position = c(0.03, 0.985),
                            legend_justification = c(0, 1)) {
  
  col_vec <- setNames(
    c("#2986cc", "#9fc5e8", "#999999", "#ea9999", "#cb271b", "#000000"),
    legend_levels
  )
  ltype_vec <- setNames(c(rep("solid", 5), "dotted"), legend_levels)
  
  fill_vec <- col_vec
  lwd_vec  <- setNames(c(rep(0.8, 5), 0.8), legend_levels)
  
  sf_noncarrier <- sf %>% filter(legend_group != "Carrier")
  sf_carrier    <- sf %>% filter(legend_group == "Carrier")
  
  p <- ggplot() +
    geom_ribbon(
      data = sf_noncarrier,
      aes(
        x = time,
        ymin = pmax(cumhaz_low, 0),
        ymax = pmin(cumhaz_high, y_max),
        fill = legend_group,
        group = legend_group
      ),
      alpha = 0.10,
      color = NA,
      show.legend = FALSE
    )
  
  if (show_carrier_ci) {
    p <- p +
      geom_ribbon(
        data = sf_carrier,
        aes(
          x = time,
          ymin = pmax(cumhaz_low, 0),
          ymax = pmin(cumhaz_high, y_max),
          fill = legend_group,
          group = legend_group
        ),
        alpha = carrier_ci_alpha,
        color = NA,
        show.legend = FALSE
      )
  }
  
  p +
    geom_line(
      data = sf_noncarrier,
      aes(
        x = time,
        y = cumhaz,
        color = legend_group,
        linewidth = legend_group,
        linetype = legend_group,
        group = legend_group
      )
    ) +
    geom_line(
      data = sf_carrier,
      aes(
        x = time,
        y = cumhaz,
        color = legend_group,
        linewidth = legend_group,
        linetype = legend_group,
        group = legend_group
      )
    ) +
    scale_color_manual(
      values = col_vec,
      breaks = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    scale_fill_manual(
      values = fill_vec,
      breaks = legend_levels,
      drop = FALSE,
      guide = "none"
    ) +
    scale_linewidth_manual(
      values = lwd_vec,
      breaks = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    scale_linetype_manual(
      values = ltype_vec,
      breaks = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 16),
      breaks = c(0, 4, 8, 12, 16),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, y_max, by = 0.02),
      labels = y_label_formatter,
      expand = c(0, 0),
      limits = c(0, y_max)
    ) +
    labs(
      x = x_lab_common,
      y = y_lab_common,
      title = plot_title
    ) +
    base_panel_theme() +
    theme(
      legend.position = legend_position,
      legend.justification = legend_justification,
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.5, family = "Arial"),
      legend.key.width = unit(1.1, "lines"),
      legend.key.height = unit(0.75, "lines"),
      legend.spacing.y = unit(0.05, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA)
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          linewidth = 0.8,
          linetype = c("solid", "solid", "solid", "solid", "solid", "dotted")
        )
      ),
      linewidth = "none"
    ) +
    theme(
      plot.margin = margin(8, 8, 8, 8)
    )
}

# =========================================================
# Build prediction data for final 4-panel figure
# =========================================================
sf_prs_only    <- build_prs_only_sf(df1)
sf_rs3798220   <- build_prs_variant_sf(df1, rs3798220)
sf_rs730882080 <- build_prs_variant_sf(df1, rs730882080)
sf_rs58757394  <- build_prs_variant_sf(df1, rs58757394)

# =========================================================
# Build final 4 figure panels
# =========================================================
p_A <- make_panel_plot(
  sf = sf_prs_only,
  plot_title = "<b>PRS only</b>",
  y_max = 0.08,
  show_carrier_ci = FALSE,
  carrier_ci_alpha = 0.04,
  legend_position = c(0.03, 0.97),
  legend_justification = c(0, 1)
)

p_B <- make_panel_plot(
  sf = sf_rs3798220,
  plot_title = "<b><i>LPA</i> (rs3798220_C)</b>",
  y_max = 0.08,
  show_carrier_ci = TRUE,
  carrier_ci_alpha = 0.04,
  legend_position = c(0.03, 0.97),
  legend_justification = c(0, 1)
)

p_C <- make_panel_plot(
  sf = sf_rs730882080,
  plot_title = "<b><i>LDLR</i> (rs730882080_T)</b>",
  y_max = 0.08,
  show_carrier_ci = TRUE,
  carrier_ci_alpha = 0.04,
  legend_position = c(0.03, 0.97),
  legend_justification = c(0, 1)
)

p_D <- make_panel_plot(
  sf = sf_rs58757394,
  plot_title = "<b><i>SLC4A11</i> (rs58757394_T)</b>",
  y_max = 0.08,
  show_carrier_ci = TRUE,
  carrier_ci_alpha = 0.04,
  legend_position = c(0.03, 0.97),
  legend_justification = c(0, 1)
)

# =========================================================
# Combine final 4-panel figure
# =========================================================
final_plot <- ((p_A | p_B) / (p_C | p_D)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 11, family = "Arial")
  )

final_plot

# =========================================================
# Save figure
# =========================================================
ggsave(
  "/path/to/figures/CHD_PRS_four_panels.pdf",
  final_plot,
  width = mm_to_in(180),
  height = mm_to_in(180),
  units = "in",
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  "/path/to/figures/CHD_PRS_four_panels.png",
  final_plot,
  width = mm_to_in(180),
  height = mm_to_in(180),
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

# =========================================================
# MODEL TABLES
# =========================================================
build_quintile_hr_table_prs_only <- function(df, model_name = "PRS_only") {
  
  df_model <- df %>%
    filter(
      !is.na(PRED_eurW_CHD),
      !is.na(Annipers_Chd0_Fnof_Fup2020),
      !is.na(Chd0_Fnof_Fup2020),
      !is.na(Age),
      !is.na(Sex),
      !is.na(PC1),
      !is.na(PC2),
      !is.na(PC3),
      !is.na(PC4),
      !is.na(PC5)
    ) %>%
    mutate(
      CHD_group = ntile(PRED_eurW_CHD, 5),
      CHD_group = factor(CHD_group, levels = 1:5, labels = paste0("Q", 1:5))
    )
  
  fit <- coxph(
    Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
      CHD_group + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
    data = df_model,
    ties = "efron",
    x = TRUE
  )
  
  sm <- summary(fit)
  
  coef_tab <- as.data.frame(sm$coefficients)
  coef_tab$term <- rownames(coef_tab)
  
  ci_tab <- as.data.frame(sm$conf.int)
  ci_tab$term <- rownames(ci_tab)
  
  coef_tab <- coef_tab %>%
    dplyr::rename(
      beta = coef,
      se_beta = `se(coef)`,
      z_value = z,
      P_value = `Pr(>|z|)`
    )
  
  ci_tab <- ci_tab %>%
    dplyr::rename(
      HR = `exp(coef)`,
      CI_low = `lower .95`,
      CI_high = `upper .95`
    )
  
  model_terms_tab <- coef_tab %>%
    left_join(ci_tab %>% dplyr::select(term, HR, CI_low, CI_high), by = "term") %>%
    mutate(HR_95CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high))
  
  hr_quint_tab <- model_terms_tab %>%
    filter(grepl("^CHD_group", term)) %>%
    mutate(group = sub("^CHD_group", "", term)) %>%
    dplyr::select(group, HR_95CI, P_value)
  
  ref_row <- data.frame(
    group = "Q1",
    HR_95CI = "1.00 (ref)",
    P_value = NA_real_
  )
  
  hr_quint_tab <- bind_rows(ref_row, hr_quint_tab) %>%
    mutate(group = factor(group, levels = paste0("Q", 1:5))) %>%
    arrange(group) %>%
    mutate(group = as.character(group))
  
  count_tab <- df_model %>%
    group_by(CHD_group) %>%
    summarise(
      N_total = n(),
      N_cases = sum(Chd0_Fnof_Fup2020 == 1, na.rm = TRUE),
      N_censored = sum(Chd0_Fnof_Fup2020 == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = CHD_group) %>%
    mutate(group = as.character(group))
  
  count_tab %>%
    left_join(hr_quint_tab, by = "group") %>%
    mutate(model = model_name) %>%
    dplyr::select(model, group, HR_95CI, P_value, N_cases, N_censored, N_total)
}

build_quintile_hr_table_variant <- function(df, variant_col, model_name, reference_group = "Q1") {
  
  var_sym <- rlang::sym(variant_col)
  
  df_model <- df %>%
    filter(
      !is.na(PRED_eurW_CHD),
      !is.na(Annipers_Chd0_Fnof_Fup2020),
      !is.na(Chd0_Fnof_Fup2020),
      !is.na(Age),
      !is.na(Sex),
      !is.na(PC1),
      !is.na(PC2),
      !is.na(PC3),
      !is.na(PC4),
      !is.na(PC5)
    ) %>%
    mutate(
      carrier = ifelse(!!var_sym == 1, 1, 0)
    ) %>%
    filter(!is.na(carrier))
  
  noncarrier_cut <- df_model %>%
    filter(carrier == 0) %>%
    mutate(
      CHD_quint_noncarrier = ntile(PRED_eurW_CHD, 5),
      CHD_quint_noncarrier = factor(
        CHD_quint_noncarrier,
        levels = 1:5,
        labels = paste0("Q", 1:5)
      )
    ) %>%
    select(IID, CHD_quint_noncarrier)
  
  df_model <- df_model %>%
    left_join(noncarrier_cut, by = "IID") %>%
    mutate(
      risk_group = case_when(
        carrier == 1 ~ "Carrier",
        carrier == 0 ~ as.character(CHD_quint_noncarrier),
        TRUE ~ NA_character_
      )
    )
  
  if (reference_group == "Q5") {
    fit_levels <- c("Q5", "Q1", "Q2", "Q3", "Q4", "Carrier")
  } else {
    fit_levels <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Carrier")
  }
  
  df_model <- df_model %>%
    mutate(
      risk_group = factor(risk_group, levels = fit_levels)
    )
  
  fit <- coxph(
    Surv(Annipers_Chd0_Fnof_Fup2020, Chd0_Fnof_Fup2020) ~
      risk_group + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5,
    data = df_model,
    ties = "efron",
    x = TRUE
  )
  
  sm <- summary(fit)
  
  coef_tab <- as.data.frame(sm$coefficients)
  coef_tab$term <- rownames(coef_tab)
  
  ci_tab <- as.data.frame(sm$conf.int)
  ci_tab$term <- rownames(ci_tab)
  
  coef_tab <- coef_tab %>%
    dplyr::rename(
      beta = coef,
      se_beta = `se(coef)`,
      z_value = z,
      P_value = `Pr(>|z|)`
    )
  
  ci_tab <- ci_tab %>%
    dplyr::rename(
      HR = `exp(coef)`,
      CI_low = `lower .95`,
      CI_high = `upper .95`
    )
  
  model_terms_tab <- coef_tab %>%
    left_join(ci_tab %>% dplyr::select(term, HR, CI_low, CI_high), by = "term") %>%
    mutate(HR_95CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high))
  
  hr_tab <- model_terms_tab %>%
    filter(grepl("^risk_group", term)) %>%
    mutate(group = sub("^risk_group", "", term)) %>%
    dplyr::select(group, HR_95CI, P_value)
  
  ref_row <- data.frame(
    group = reference_group,
    HR_95CI = "1.00 (ref)",
    P_value = NA_real_
  )
  
  hr_tab <- bind_rows(ref_row, hr_tab) %>%
    mutate(group = factor(group, levels = c(paste0("Q", 1:5), "Carrier"))) %>%
    arrange(group) %>%
    mutate(group = as.character(group))
  
  count_tab <- df_model %>%
    mutate(
      risk_group = factor(as.character(risk_group), levels = c(paste0("Q", 1:5), "Carrier"))
    ) %>%
    group_by(risk_group) %>%
    summarise(
      N_total = n(),
      N_cases = sum(Chd0_Fnof_Fup2020 == 1, na.rm = TRUE),
      N_censored = sum(Chd0_Fnof_Fup2020 == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = risk_group) %>%
    mutate(group = as.character(group))
  
  count_tab %>%
    left_join(hr_tab, by = "group") %>%
    mutate(model = model_name) %>%
    dplyr::select(model, group, HR_95CI, P_value, N_cases, N_censored, N_total)
}

# =========================================================
# Build final tables for all requested models
# =========================================================
tab_prs_only <- build_quintile_hr_table_prs_only(
  df = df1,
  model_name = "PRS_only"
)

tab_rs3798220 <- build_quintile_hr_table_variant(
  df = df1,
  variant_col = "rs3798220",
  model_name = "LPA_rs3798220",
  reference_group = "Q1"
)

tab_rs730882080 <- build_quintile_hr_table_variant(
  df = df1,
  variant_col = "rs730882080",
  model_name = "LDLR_rs730882080",
  reference_group = "Q1"
)

tab_rs58757394 <- build_quintile_hr_table_variant(
  df = df1,
  variant_col = "rs58757394",
  model_name = "SLC4A11_rs58757394",
  reference_group = "Q5"
)

tab_all_models <- bind_rows(
  tab_prs_only,
  tab_rs3798220,
  tab_rs730882080,
  tab_rs58757394
)

fwrite(
  tab_all_models,
  "/path/to/survival_files/CHD_all_models_quintile_HR_table_rareVar.tsv",
  sep = "\t"
)