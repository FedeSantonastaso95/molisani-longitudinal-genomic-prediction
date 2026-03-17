############################################################
# Compute SCORE2 for primary prevention participants
#
# Description:
# This script computes SCORE2 using the reference implementation
# provided by the RiskScorescvd package.
#
# Eligibility criteria:
#   - QC-passed participants only
#   - no prevalent CHD at baseline
#   - no prevalent stroke at baseline
#   - no diabetes
#   - age between 40 and 69 years
#
# Required SCORE2 inputs:
#   - Age (years)
#   - Gender ("male" / "female")
#   - smoker (0/1; current smoker no/yes)
#   - systolic.bp (mmHg)
#   - diabetes (0/1)
#   - total.chol (mmol/L)
#   - total.hdl (mmol/L)
#
# Cholesterol conversion:
#   The source dataset is assumed to store cholesterol in mg/dL.
#   Conversion to mmol/L is performed using:
#     mmol/L = mg/dL * 0.02586
#
# Output:
#   A tab-delimited file containing:
#   - Idth_Ms2022
#   - SCORE2_score
#   - SCORE2_strat
#   - FID
############################################################

## ----------------------------
## Libraries
## ----------------------------
library(RiskScorescvd)
library(data.table)
library(dplyr)
library(tibble)

## ----------------------------
## Constants
## ----------------------------
mmol_per_mgdl_chol <- 0.02586

## ----------------------------
## Input / output paths
## ----------------------------
infile_main <- "path/to/main_longitudinal_dataset.txt"
infile_aux  <- "path/to/auxiliary_phenotype_dataset.txt"
outfile     <- "path/to/output/2026.02.18_SCORE2_computed.txt"

## ----------------------------
## 1. Load main dataset
## ----------------------------
df_main <- fread(
  infile_main,
  data.table = FALSE
) 

## ----------------------------
## 2. Load auxiliary dataset
## ----------------------------
# Smoking and diabetes information are taken from the auxiliary file.
df_aux <- fread(
  infile_aux,
  data.table = FALSE
) %>%
  select(
    Idth_Ms2022,
    Smoking_Tabacco,
    T2d_Arw
  )

## ----------------------------
## 3. Merge datasets
## ----------------------------
df <- df_main %>%
  left_join(df_aux, by = "Idth_Ms2022")

## ----------------------------
## 4. Define SCORE2-eligible sample
## ----------------------------
df_score <- df %>%
  select(
    Idth_Ms2022,
    Smoking_Tabacco,
    SBP,
    HDL,
    Age,
    Sex,
    Cholesterol_TOT,
    Chd_Bs_New,
    Cerebro_Bs_New,
    T2d_Arw
  ) %>%
  # Exclude diabetes
  filter(T2d_Arw != 1) %>%
  # Exclude prevalent CHD
  filter(!Chd_Bs_New %in% c(1, 2)) %>%
  # Exclude prevalent stroke
  filter(!Cerebro_Bs_New %in% c(1, 2)) %>%
  # Keep SCORE2 age range
  filter(Age >= 40, Age <= 69)

## ----------------------------
## 5. Prepare SCORE2 input variables
## ----------------------------
df_score <- df_score %>%
  mutate(
    Gender = ifelse(Sex == 0, "female", "male"),
    smoker = ifelse(Smoking_Tabacco == 1, 1, 0),
    diabetes = ifelse(T2d_Arw == 1, 1, 0),
    total.hdl = HDL * mmol_per_mgdl_chol,
    total.chol = Cholesterol_TOT * mmol_per_mgdl_chol,
    systolic.bp = SBP
  )

## ----------------------------
## 6. Keep only variables required by SCORE2
## ----------------------------
df_input <- df_score %>%
  select(
    Idth_Ms2022,
    Age,
    Gender,
    smoker,
    systolic.bp,
    diabetes,
    total.chol,
    total.hdl
  )

## ----------------------------
## 7. Compute continuous SCORE2
## ----------------------------
df_score2_value <- df_input %>%
  rowwise() %>%
  mutate(
    SCORE2_score = RiskScorescvd::SCORE2(
      Risk.region = "Moderate",
      Age = Age,
      Gender = Gender,
      smoker = smoker,
      systolic.bp = systolic.bp,
      diabetes = diabetes,
      total.chol = total.chol,
      total.hdl = total.hdl,
      classify = FALSE
    )
  ) %>%
  ungroup()

## ----------------------------
## 8. Compute SCORE2 category
## ----------------------------
df_score2_class <- df_input %>%
  rowwise() %>%
  mutate(
    SCORE2_strat = RiskScorescvd::SCORE2(
      Risk.region = "Moderate",
      Age = Age,
      Gender = Gender,
      smoker = smoker,
      systolic.bp = systolic.bp,
      diabetes = diabetes,
      total.chol = total.chol,
      total.hdl = total.hdl,
      classify = TRUE
    )
  ) %>%
  ungroup()

## ----------------------------
## 9. Merge SCORE2 score and category
## ----------------------------
df_score2 <- df_score2_value %>%
  select(Idth_Ms2022, SCORE2_score) %>%
  left_join(
    df_score2_class %>% select(Idth_Ms2022, SCORE2_strat),
    by = "Idth_Ms2022"
  )

## ----------------------------
## 10. Merge back to all QC-passed participants
## ----------------------------
# This keeps all QC-passed individuals in the final output and leaves
# SCORE2 variables as NA for those not eligible for SCORE2 calculation.
df_ids <- df_main %>%
  select(Idth_Ms2022)

df_output <- df_ids %>%
  left_join(df_score2, by = "Idth_Ms2022")

## ----------------------------
## 11. Add FID
## ----------------------------
# FID is added back using the same QC-passed main dataset.
df_fid <- df_main %>%
  mutate(
    Chd_Bs_New_final = ifelse(Chd_Bs_New == 0, 0, 1)
  ) %>%
  filter(!(Chd_Bs_New_final == 1 & Chd0_Fnof_Fup2020 %in% c(0, 1))) %>%
  select(Idth_Ms2022, FID)

df_output <- df_output %>%
  left_join(df_fid, by = "Idth_Ms2022")

## ----------------------------
## 12. Save output
## ----------------------------
fwrite(
  df_output,
  file = outfile,
  na = "NA",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE
)