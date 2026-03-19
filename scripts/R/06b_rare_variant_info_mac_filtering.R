################################################################################
# 06b_rare_variant_info_mac_filtering.R
#
# Description:
# This script annotates candidate rare variants with imputation quality (INFO)
# and allele frequency (A1FREQ) using association summary statistic files
# generated with REGENIE from sequence-imputed genotype data.
#
# The input legend file contains, for each trait, the path to the corresponding
# REGENIE summary statistics file (one file per trait).
#
# REGENIE summary statistics files are expected to include at least the
# following columns:
#   - CHROM      : chromosome
#   - GENPOS     : genomic position (GRCh38)
#   - ALLELE0    : reference allele
#   - ALLELE1    : alternate allele
#   - A1FREQ     : allele frequency of ALLELE1
#   - INFO       : imputation quality score
#
# For each variant-trait combination, the script:
#   1. matches the variant to the corresponding REGENIE summary statistics file
#   2. extracts INFO and A1FREQ
#   3. computes minor allele frequency (MAF) as:
#         MAF = min(A1FREQ, 1 - A1FREQ)
#   4. computes minor allele count (MAC) as:
#         MAC = 2 * N * MAF
#   5. applies quality filters:
#         INFO >= 0.8
#         MAC  >= 5
#
# Notes:
# - INFO and A1FREQ are derived from REGENIE association analyses performed on
#   the same imputed genotype dataset used for rare variant analyses.
# - MAF is obtained by folding A1FREQ to the minor allele frequency scale.
# - Paths shown below are placeholders for reproducibility.
#
# Inputs:
#   - Candidate rare variant table
#   - Legend file mapping traits to REGENIE summary statistics paths
#
# Outputs:
#   - Annotated table with INFO, A1FREQ, MAF, and MAC
#   - Filtered table retaining variants with INFO >= 0.8 and MAC >= 5
#
# Usage:
#   Rscript 06b_rare_variant_info_mac_filtering.R
#
################################################################################

# =========================
# Load libraries
# =========================
library(data.table)
library(dplyr)

# =========================
# Define input/output paths
# =========================
candidate_file        <- "/path/to/rare_variants_candidates.tsv"
legend_file           <- "/path/to/regenie_sumstats_legend_MAC1.tsv"
annotated_output_file <- "/path/to/rare_variants_with_info_maf_mac.tsv"
filtered_output_file  <- "/path/to/rare_variants_filtered_info08_mac5.tsv"

# If N is not already present in the candidate table, define a fixed sample size
# here. Otherwise, the script will use the N column from the candidate table.
study_n <- NA_real_

# =========================
# Load inputs
# =========================
dc <- fread(candidate_file)
legend_imp <- fread(legend_file, data.table = FALSE)

# The legend file is expected to already contain paths to REGENIE summary
# statistics computed on sequence-imputed genotypes (MAC1 dataset).

setDT(dc)
setDT(legend_imp)

# =========================
# Basic checks
# =========================
required_dc_cols <- c("trait", "variant")
missing_dc_cols <- setdiff(required_dc_cols, names(dc))
if (length(missing_dc_cols) > 0) {
  stop(
    "Missing required columns in candidate_file: ",
    paste(missing_dc_cols, collapse = ", ")
  )
}

required_legend_cols <- c("trait", "original_folder")
missing_legend_cols <- setdiff(required_legend_cols, names(legend_imp))
if (length(missing_legend_cols) > 0) {
  stop(
    "Missing required columns in legend_file: ",
    paste(missing_legend_cols, collapse = ", ")
  )
}

# =========================
# Restrict to overlapping traits
# =========================
traits_use <- intersect(unique(dc$trait), unique(legend_imp$trait))
legend_use <- legend_imp[trait %chin% traits_use, .(trait, original_folder)]

if (length(traits_use) == 0) {
  stop("No overlapping traits found between candidate table and legend file.")
}

# =========================
# Initialise annotation columns
# =========================
dc[, INFO := as.numeric(NA)]
dc[, A1FREQ := as.numeric(NA)]

# =========================
# Helper function to create variant key
# =========================
# Expected key format:
#   CHR-POS-ALLELE0-ALLELE1
make_key <- function(chr, pos, a0, a1) {
  paste0(chr, "-", pos, "-", toupper(a0), "-", toupper(a1))
}

# =========================
# Annotate INFO and A1FREQ by trait
# =========================
for (tr in legend_use$trait) {
  
  f <- legend_use[trait == tr, original_folder][1]
  
  if (!file.exists(f)) {
    warning("REGENIE file not found for trait ", tr, ": ", f)
    next
  }
  
  ss <- fread(f)
  
  required_ss_cols <- c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "INFO", "A1FREQ")
  missing_ss_cols <- setdiff(required_ss_cols, names(ss))
  
  if (length(missing_ss_cols) > 0) {
    stop(
      "Missing required columns in REGENIE summary statistics file for trait ",
      tr, ": ", paste(missing_ss_cols, collapse = ", ")
    )
  }
  
  ss <- ss[, ..required_ss_cols]
  
  # Build matching key in REGENIE summary statistics
  ss[, key := make_key(CHROM, GENPOS, ALLELE0, ALLELE1)]
  ss <- ss[, .(key, INFO, A1FREQ)]
  setkey(ss, key)
  
  # Subset candidate rows for current trait
  idx <- which(dc$trait == tr)
  
  # Match candidate variants to REGENIE summary statistics
  m <- ss[.(dc$variant[idx]), on = "key", nomatch = NA]
  
  # Fill INFO and A1FREQ in the candidate table
  dc$INFO[idx]   <- m$INFO
  dc$A1FREQ[idx] <- m$A1FREQ
}

# =========================
# Compute MAF and MAC
# =========================
# MAF is derived from A1FREQ by folding to the minor allele frequency scale.
dc[, MAF := fifelse(
  is.na(A1FREQ),
  NA_real_,
  pmin(A1FREQ, 1 - A1FREQ)
)]

# Use N from candidate table if available; otherwise use study_n
if ("N" %in% names(dc)) {
  dc[, N_for_MAC := as.numeric(N)]
} else {
  dc[, N_for_MAC := study_n]
}

dc[, MAC := 2 * N_for_MAC * MAF]

# =========================
# Save annotated table
# =========================
fwrite(
  dc,
  file = annotated_output_file,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# =========================
# Apply INFO and MAC filters
# =========================
filtered_dc <- dc %>%
  filter(!is.na(INFO), !is.na(MAF), !is.na(MAC)) %>%
  filter(INFO >= 0.8, MAC >= 5)

# =========================
# Save filtered table
# =========================
fwrite(
  filtered_dc,
  file = filtered_output_file,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)