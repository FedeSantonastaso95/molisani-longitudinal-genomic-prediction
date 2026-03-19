################################################################################
# 06a_rare_variant_lovo_processing.R
#
# Description:
# This script processes leave-one-variant-out (LOVO) rare variant association
# results generated with REGENIE step 2.
#
# The script performs two main tasks:
#   1. Reads all REGENIE LOVO output files from a specified directory and
#      merges them into a single table, extracting metadata from file names
#      (Gene, Mask, Threshold, Trait, Model).
#   2. Post-processes the merged LOVO results to quantify the contribution of
#      individual variants within each associated gene-based signal.
#
# Specifically, the script:
#   - parses REGENIE output files
#   - combines all LOVO result files into a unified table
#   - retains only the burden test corresponding to the top gene-level signal
#     identified in the reference burden results table
#   - computes the variance explained proxy as:
#         r2 = CHISQ / N
#   - identifies the reference (non-LOVO) result within each group, defined as
#     the gene-based burden test including all rare variants in the mask
#     (i.e. without leaving any variant out)
#   - computes the relative reduction in variance explained after excluding
#     each variant:
#         delta_r2 = (ref_r2 - r2) / ref_r2
#
# LOVO results are evaluated in two steps:
#   1. Technical validity:
#      results are considered technically valid if they have:
#      - a unique reference row
#      - a significant reference association signal (LOG10P >= 5.6)
#      - a non-negative delta_r2
#
#   2. Variant contribution:
#      among technically valid LOVO results, variants are considered
#      contributory when their exclusion leads to a reduction in variance
#      explained of at least 10% relative to the corresponding full gene-based
#      burden test:
#          delta_r2 >= 0.10
#
# Status labels:
#   - "valid"                 -> LOVO result passing technical validity checks
#   - "ref"                   -> reference (non-LOVO) row
#   - "no_ref_or_non_unique"  -> missing or non-unique reference row
#   - "low_log10p"            -> reference row does not meet significance threshold
#   - "delta_r2_le_0"         -> invalid or non-informative delta_r2
#
# Contribution labels:
#   - "reference"             -> reference (non-LOVO) row
#   - "contributory"          -> technically valid LOVO result with delta_r2 >= 0.10
#   - "non_contributory"      -> technically valid LOVO result with delta_r2 < 0.10
#   - NA                      -> row not eligible for contribution assessment
#
# Inputs:
#   - A directory containing REGENIE LOVO output files (.regenie)
#   - A reference table listing the highest LOG10P burden-test result for each
#     gene/mask/threshold/trait combination
#
# Outputs:
#   - A merged LOVO results table
#   - A post-processed LOVO summary table including r2, core_id, delta_r2,
#     technical validity status, and contribution status
#   - A filtered table containing contributory variants only
#
# Notes:
# - REGENIE output file names are expected to follow the pattern:
#     GENE-MASK-THRESHOLD-TRAIT-MODEL.regenie
#
# - The reference (non-LOVO) result corresponds to the gene-based burden test
#   including all variants in the mask. It is identified as the row where
#   ID_main matches core_id.
#
# - The significance threshold for the reference signal is currently set to:
#     LOG10P >= 5.6
#
# - Paths shown below are placeholders for reproducibility. Original analyses
#   were run on protected institutional infrastructure.
#
# Usage:
#   Rscript 06a_rare_variant_lovo_processing.R
#
################################################################################

# =========================
# Load libraries
# =========================
library(data.table)
library(dplyr)
library(tidyr)

# =========================
# Define input/output paths
# =========================
lovo_results_dir            <- "/path/to/lovo_results"
reference_file              <- "/path/to/highest_log10p_reference.tsv"
merged_output_file          <- "/path/to/lovo_merged_filtered.tsv"
final_output_file           <- "/path/to/lovo_postprocessed.tsv"
contributory_output_file    <- "/path/to/lovo_contributory_variants.tsv"

# =========================
# Helper function to read and annotate REGENIE output files
# =========================
read_and_annotate <- function(file_path) {
  
  file_name <- basename(file_path)
  file_name <- sub("\\.regenie$", "", file_name)
  
  # Expected format:
  # GENE-MASK-THRESHOLD-TRAIT-MODEL.regenie
  parts <- unlist(strsplit(file_name, "-"))
  
  if (length(parts) < 5) {
    stop(paste("Unexpected filename format:", file_name))
  }
  
  gene <- parts[1]
  mask <- parts[2]
  threshold <- parts[3]
  trait <- parts[4]
  model <- paste(parts[5:length(parts)], collapse = "-")
  
  df <- tryCatch(
    fread(
      file_path,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      data.table = FALSE
    ),
    error = function(e) NULL
  )
  
  # If file is empty or unreadable, create a 1-row NA data frame with header
  if (is.null(df) || nrow(df) == 0) {
    header_only <- fread(
      file_path,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      nrows = 0,
      data.table = FALSE
    )
    
    df <- as.data.frame(matrix(NA, ncol = ncol(header_only), nrow = 1))
    colnames(df) <- colnames(header_only)
  }
  
  df <- df %>%
    mutate(across(everything(), as.character)) %>%
    mutate(
      Gene = gene,
      Mask = mask,
      Threshold = threshold,
      Trait = trait,
      Model = model
    )
  
  return(df)
}

# =========================
# Read and merge LOVO REGENIE output files
# =========================
file_list <- list.files(
  lovo_results_dir,
  pattern = "MAC1.*\\.regenie$",
  full.names = TRUE
)

if (length(file_list) == 0) {
  stop("No REGENIE LOVO output files found in lovo_results_dir.")
}

data_list <- lapply(file_list, read_and_annotate)

# Keep only columns shared across all files
common_cols <- Reduce(intersect, lapply(data_list, colnames))
data_list <- lapply(data_list, function(df) df[, common_cols, drop = FALSE])

merged_data <- bind_rows(data_list)

# =========================
# Parse variant ID fields
# =========================
merged_data_parsed <- merged_data %>%
  separate(
    ID,
    into = c("ID_main", "ChrPos", "ALLELE0", "ALLELE1"),
    sep = ":",
    extra = "drop",
    fill = "right"
  ) %>%
  dplyr::select(-ChrPos)

# =========================
# Load reference burden results
# =========================
reference_dt <- fread(reference_file)

reference_dt <- reference_dt %>%
  dplyr::select(trait, GENE, MASK, THRESHOLD, TEST, gene_name)

# =========================
# Retain only the burden test of interest
# =========================
filtered_data <- merged_data_parsed %>%
  inner_join(
    reference_dt,
    by = c(
      "Gene" = "GENE",
      "Mask" = "MASK",
      "Threshold" = "THRESHOLD",
      "Trait" = "trait"
    )
  ) %>%
  filter(TEST.x == TEST.y) %>%
  dplyr::select(-TEST.y) %>%
  rename(TEST = TEST.x)

# Save merged and filtered LOVO table
fwrite(
  filtered_data,
  merged_output_file,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# =========================
# Compute variance explained and post-process LOVO results
# =========================
df <- copy(filtered_data)
setDT(df)

# Ensure numeric columns are numeric
numeric_cols <- c("CHISQ", "N", "LOG10P")
for (col in intersect(numeric_cols, names(df))) {
  set(df, j = col, value = as.numeric(df[[col]]))
}

# Compute variance explained proxy
df[, r2 := CHISQ / N]

# Extract core ID by removing LOVO suffix
df[, core_id := sub("_chr.*", "", ID_main)]

# Identify groups with at least one exact reference row
core_ids <- df[ID_main == core_id, unique(core_id)]

split_list <- list()

for (cid in core_ids) {
  
  cid_group <- df[core_id == cid]
  
  if (nrow(cid_group) == 0) next
  
  for (trait_name in unique(cid_group$Trait)) {
    
    trait_group <- cid_group[Trait == trait_name]
    
    for (test_name in unique(trait_group$TEST)) {
      
      group_dt <- copy(trait_group[TEST == test_name])
      ref_row <- group_dt[ID_main == cid]
      
      # Default technical status
      group_dt[, status := "valid"]
      
      # Default contribution status
      group_dt[, contribution_status := NA_character_]
      
      # Check that the reference row exists and is unique
      if (nrow(ref_row) != 1) {
        group_dt[, status := "no_ref_or_non_unique"]
        split_list[[paste(cid, trait_name, test_name, sep = "|")]] <- group_dt
        next
      }
      
      # Require significant reference signal
      if (is.na(ref_row$LOG10P) || ref_row$LOG10P < 5.6) {
        group_dt[, status := "low_log10p"]
        split_list[[paste(cid, trait_name, test_name, sep = "|")]] <- group_dt
        next
      }
      
      # Mark reference row
      group_dt[ID_main == cid, `:=`(
        status = "ref",
        contribution_status = "reference"
      )]
      
      # Compute delta_r2 relative to the full burden test
      ref_r2 <- ref_row$r2
      
      if (is.na(ref_r2) || ref_r2 <= 0) {
        group_dt[ID_main != cid, status := "delta_r2_le_0"]
        split_list[[paste(cid, trait_name, test_name, sep = "|")]] <- group_dt
        next
      }
      
      group_dt[, delta_r2 := (ref_r2 - r2) / ref_r2]
      
      # Flag invalid or non-informative LOVO rows
      group_dt[ID_main != cid & (is.na(delta_r2) | delta_r2 < 0), status := "delta_r2_le_0"]
      
      # Assign contribution status only to technically valid LOVO rows
      group_dt[ID_main != cid & status == "valid" & delta_r2 >= 0.10,
               contribution_status := "contributory"]
      
      group_dt[ID_main != cid & status == "valid" & delta_r2 < 0.10,
               contribution_status := "non_contributory"]
      
      split_list[[paste(cid, trait_name, test_name, sep = "|")]] <- group_dt
    }
  }
}

# Combine all groups into final table
combined_df <- rbindlist(split_list, use.names = TRUE, fill = TRUE)

# =========================
# Save final post-processed LOVO table
# =========================
fwrite(
  combined_df,
  final_output_file,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)

# =========================
# Save contributory variants only
# =========================
contributory_df <- combined_df %>%
  filter(contribution_status == "contributory")

fwrite(
  contributory_df,
  contributory_output_file,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)