#!/bin/bash

################################################################################
# extract_rare_variant_dosages.sh
#
# Description:
# This script extracts genotype dosages for a set of selected rare variants
# from sequence-imputed genotype data using PLINK2.
#
# Genotypes are exported using the PLINK2 option `--export A`, which outputs
# the dosage of the counted allele (A1) for each variant and individual.
#
# IMPORTANT:
# The exported dosage corresponds to the A1 allele in PLINK and does not
# necessarily represent the minor (rare) allele. In this dataset, A1 often
# corresponds to the common allele. Therefore, rare-allele carrier status
# must be derived downstream by combining genotype dosages with allele
# frequency information (e.g. A1FREQ from REGENIE summary statistics).
#
# Inputs:
#   - Imputed genotype dataset (BGEN + SAMPLE)
#   - List of variants to extract (one variant ID per line)
#
# Outputs:
#   - PLINK2 .raw file containing allele dosages
#
# Usage:
#   bash extract_rare_variant_dosages.sh
#
################################################################################

# =========================
# Load software
# =========================
module load plink/2.00_20211217

# =========================
# Define input/output paths
# =========================
BGEN_FILE="/path/to/imputed_genotypes.bgen"
SAMPLE_FILE="/path/to/imputed_genotypes.sample"
VARIANT_LIST="/path/to/rare_variants_to_extract.txt"
OUTPUT_PREFIX="/path/to/carriers_rare_variants"

# =========================
# Run PLINK2
# =========================
plink2 \
  --bgen "${BGEN_FILE}" \
  ref-first \
  --sample "${SAMPLE_FILE}" \
  --extract "${VARIANT_LIST}" \
  --export A \
  --out "${OUTPUT_PREFIX}"

# =========================
# Done
# =========================
echo "Extraction completed. Output saved to: ${OUTPUT_PREFIX}.raw"
