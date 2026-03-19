# =========================================================
# 09_rare_variant_carrier_phenotype_distribution.R
#
# Description:
# This script generates a supplementary figure showing the distribution
# of selected quantitative phenotypes in carriers versus non-carriers
# of rare variants associated with CHD-related traits.
#
# Variants shown:
# - rs730882080_T  -> LDL cholesterol
# - rs58757394_T   -> systolic blood pressure
# - rs3798220_C    -> lipoprotein(a)
#
# IMPORTANT:
# - Genotype dosages were exported with PLINK2 using `--export A`.
# - In this dataset, the exported dosage corresponds to the major allele.
# - Therefore, dosages are recoded so that:
#     0 = non-carrier of the rare allele
#     1 = heterozygous carrier
#     2 = homozygous carrier of the rare allele
#
# Output:
# - A supplementary violin plot figure with one panel per variant/phenotype pair
# - Descriptive statistics printed to the console
# =========================================================

# =========================================================
# Libraries
# =========================================================
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(
  theme_classic(base_size = 14) +
    theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
)

# =========================================================
# Define input/output paths
# =========================================================
carriers_file   <- "/path/to/carriers_rarVar.raw"
sample_file     <- "/path/to/imputed_genotypes.sample"
phenotype_file  <- "/path/to/phenotypes_covariates_prs.txt"
figure_file_png <- "/path/to/figures/supp_rare_variant_violin_panels.png"
figure_file_pdf <- "/path/to/figures/supp_rare_variant_violin_panels.pdf"

# =========================================================
# Load carrier file
# =========================================================
df_carr <- fread(carriers_file)

# Round all dosage columns
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

# Clean variant column names by removing the suffix after the last underscore
old_names <- names(df_carr)[7:ncol(df_carr)]
new_names <- sub("_(?!.*_).*$", "", old_names, perl = TRUE)
names(df_carr)[7:ncol(df_carr)] <- new_names

# =========================================================
# Load sample file
# =========================================================
sample <- fread(sample_file)

# Remove the first line of the PLINK/BGEN sample file if it contains
# the standard dummy row
sample <- sample[-1, ]

# =========================================================
# Load phenotype data
# =========================================================
df <- fread(
  phenotype_file,
  data.table = FALSE
) %>%
  left_join(df_carr, by = c("FID" = "FID", "IID" = "IID"))

# Restrict to sequence-imputed samples
df1 <- sample %>%
  left_join(df, by = c("ID_1" = "FID"))

# Keep only variables needed for the supplementary figure
df1 <- df1 %>%
  select(
    IID,
    rs3798220,
    rs58757394,
    rs730882080,
    Bio_Lp_A,
    SBP,
    LDL
  )

# =========================================================
# Recode carrier status for plotting
# =========================================================
# Here we collapse genotypes as:
#   0 = non-carrier
#   1 or 2 = carrier
df1 <- df1 %>%
  mutate(
    carrier_rs730882080 = factor(
      ifelse(rs730882080 >= 1, 1, 0),
      levels = c(0, 1),
      labels = c("Non-carrier", "Carrier")
    ),
    carrier_rs58757394 = factor(
      ifelse(rs58757394 >= 1, 1, 0),
      levels = c(0, 1),
      labels = c("Non-carrier", "Carrier")
    ),
    carrier_rs3798220 = factor(
      ifelse(rs3798220 >= 1, 1, 0),
      levels = c(0, 1),
      labels = c("Non-carrier", "Carrier")
    )
  )

# =========================================================
# Plotting helpers
# =========================================================
carrier_cols <- c(
  "Non-carrier" = "grey70",
  "Carrier" = "#d73027"
)

ci_bar <- stat_summary(
  fun.data = mean_cl_boot,
  geom = "errorbar",
  width = 0.18,
  linewidth = 1.1,
  color = "black"
)

ci_point <- stat_summary(
  fun = mean,
  geom = "point",
  size = 3.5,
  color = "black"
)

# =========================================================
# Descriptive statistics helper
# =========================================================
report_stats <- function(data, genotype, phenotype, phenotype_label) {
  carrier_var <- ifelse(data[[genotype]] >= 1, "Carrier", "Non-carrier")
  
  out <- data.frame(
    group = carrier_var,
    phenotype = data[[phenotype]]
  ) %>%
    filter(!is.na(group), !is.na(phenotype)) %>%
    group_by(group) %>%
    summarise(
      N = n(),
      Mean = mean(phenotype, na.rm = TRUE),
      SD = sd(phenotype, na.rm = TRUE),
      Median = median(phenotype, na.rm = TRUE),
      IQR = IQR(phenotype, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n====================================================\n")
  cat("Variant:", genotype, "\n")
  cat("Phenotype:", phenotype_label, "\n")
  print(out)
  cat("====================================================\n")
}

# =========================================================
# Build violin plots
# =========================================================

# Panel A: LDL ~ rs730882080
dA <- df1 %>%
  filter(!is.na(LDL), !is.na(carrier_rs730882080))

pA <- ggplot(dA, aes(carrier_rs730882080, LDL, fill = carrier_rs730882080)) +
  geom_violin(trim = FALSE, color = "grey30", linewidth = 0.3, alpha = 0.95) +
  ci_bar +
  ci_point +
  scale_fill_manual(values = carrier_cols) +
  labs(
    x = "rs730882080_T",
    y = "LDL cholesterol (mg/dL)"
  )

# Panel B: SBP ~ rs58757394
dB <- df1 %>%
  filter(!is.na(SBP), !is.na(carrier_rs58757394))

pB <- ggplot(dB, aes(carrier_rs58757394, SBP, fill = carrier_rs58757394)) +
  geom_violin(trim = FALSE, color = "grey30", linewidth = 0.3, alpha = 0.95) +
  ci_bar +
  ci_point +
  scale_fill_manual(values = carrier_cols) +
  labs(
    x = "rs58757394_T",
    y = "Systolic blood pressure (mmHg)"
  )

# Panel C: Lp(a) ~ rs3798220
dC <- df1 %>%
  filter(!is.na(Bio_Lp_A), !is.na(carrier_rs3798220))

pC <- ggplot(dC, aes(carrier_rs3798220, Bio_Lp_A, fill = carrier_rs3798220)) +
  geom_violin(trim = FALSE, color = "grey30", linewidth = 0.3, alpha = 0.95) +
  ci_bar +
  ci_point +
  scale_fill_manual(values = carrier_cols) +
  labs(
    x = "rs3798220_C",
    y = "Lipoprotein(a) (mg/dL)"
  )

# =========================================================
# Combine panels into final supplementary figure
# =========================================================
final_plot <- (pA | pB | pC) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18, family = "Arial"),
    plot.tag.position = c(0.02, 0.98)
  )

print(final_plot)

# =========================================================
# Save figure
# =========================================================
ggsave(
  figure_file_png,
  final_plot,
  width = 8,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  figure_file_pdf,
  final_plot,
  width = 8,
  height = 4.5,
  bg = "white"
)

# =========================================================
# Print descriptive statistics to console
# =========================================================
report_stats(
  df1,
  genotype = "rs730882080",
  phenotype = "LDL",
  phenotype_label = "LDL cholesterol (mg/dL)"
)

report_stats(
  df1,
  genotype = "rs58757394",
  phenotype = "SBP",
  phenotype_label = "Systolic blood pressure (mmHg)"
)

report_stats(
  df1,
  genotype = "rs3798220",
  phenotype = "Bio_Lp_A",
  phenotype_label = "Lipoprotein(a) (mg/dL)"
)