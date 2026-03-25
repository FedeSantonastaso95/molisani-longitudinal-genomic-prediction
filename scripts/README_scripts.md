### `01_prs_association_models.R`

Evaluates the association between standardised polygenic risk scores (PRSs) and disease outcomes using logistic regression models.

Analyses are performed using three outcome definitions:
- all events (prevalent or incident)
- prevalent cases only
- incident cases among baseline-free individuals

Models are adjusted for age, sex (when applicable), and genetic principal components (PC1–PC5).

**Used for**
- Main PRS association analyses


### `02_figure2_and_supplementary_figureS1.R`

Generates plots for Figure 2 and Supplementary Figure S1.

The script includes:
- forest plots comparing PRS effect estimates from the Moli-sani cohort with external European reference studies
- PRS quintile-based association plots used in Figure 2B

**Used for**
- Figure 2
- Supplementary Figure S1

**Notes**
- Requires previously generated PRS association results
- File paths should be adapted locally

### `03_survival_analysis_prs.R`

Performs longitudinal survival analyses to evaluate the association between polygenic risk scores (PRS) and incident disease risk.

The script:
- prepares endpoint-specific survival datasets
- fits Cox proportional hazards models using:
  - continuous standardized PRS
  - PRS quintiles (Q1 as reference)
- estimates hazard ratios and confidence intervals
- exports results tables for continuous and quintile-based models
- reports case counts after exclusion of baseline-prevalent individuals

**Used for**
- Main longitudinal analyses

**Notes**
- Baseline-prevalent cases are excluded from all analyses
- Models are adjusted for age, sex (when applicable), and genetic principal components (PC1–PC5)


### `03_survival_analysis_prs.R`

Performs longitudinal survival analyses to evaluate the association between polygenic risk scores (PRS) and incident disease risk.

The script:
- prepares endpoint-specific survival datasets
- fits Cox proportional hazards models using:
  - continuous standardized PRS
  - PRS quintiles (Q1 as reference)
- estimates hazard ratios and confidence intervals
- exports results tables for continuous and quintile-based models
- reports case counts after exclusion of baseline-prevalent individuals

**Used for**
- Main longitudinal analyses

**Notes**
- Baseline-prevalent cases are excluded from all analyses
- Models are adjusted for age, sex (when applicable), and genetic principal components (PC1–PC5)


### `03a_make_figure3A_survival_curves.R`

Generates empirical survival curves for Figure 3A across PRS quintiles for:
- coronary heart disease (CHD)
- breast cancer (BC)
- prostate cancer (PC)
- colorectal cancer (CRC)

Curves are estimated using empirical `survfit()` models and plotted as cumulative hazard.

The script also:
- estimates hazard ratios across PRS quintiles (Q2–Q5 vs Q1)
- computes hazard ratios per 1-SD increase in PRS
- evaluates trend effects across quintiles
- exports summary tables including case counts, person-years, and incidence rates

**Used for**
- Figure 3A

**Notes**
- Baseline-prevalent cases are excluded from all analyses
- Models are adjusted for age, sex (when applicable), and genetic principal components (PC1–PC5)


### `03b_make_supplFig2_survival_curves.R`

Generates empirical survival curves for Supplementary Figure 2 across PRS quintiles for:
- ischemic stroke
- lung cancer
- renal cancer
- pancreatic cancer

Curves are estimated using empirical `survfit()` models and plotted as cumulative hazard.

The script:
- builds endpoint-specific survival datasets including incident events only
- fits Cox proportional hazards models to estimate:
  - hazard ratios per 1-SD increase in PRS
  - hazard ratios across PRS quintiles (Q2–Q5 vs Q1)
- exports summary tables used for the figure

**Used for**
- Supplementary Figure 2

**Notes**
- Baseline-prevalent cases are excluded from all analyses
- Models are adjusted for age, sex (when applicable), and genetic principal components (PC1–PC5)


### `03c_compare_breast_prs_with_without_BRCA1_BRCA2.R`

Performs breast cancer risk analyses comparing polygenic risk scores (PRS) computed with and without variants in the BRCA1 and BRCA2 loci.

The script:
- evaluates the impact of excluding high-penetrance loci on PRS performance
- estimates risk across PRS strata
- supports sensitivity analyses for breast cancer

The analysis is restricted to events observed up to age 80 and uses age as the time scale.

**Used for**
- Breast cancer sensitivity analyses

**Notes**
- The number of breast cancer cases differs from the full dataset because:
  - only events occurring up to age 80 are included
  - individuals without age-at-onset information cannot be included
- The analysis dataset includes 435 cases and 9,191 censored individuals


### `03d_breast_prs_noBRCA1_2_prevalence_by_age.R`

Estimates age-specific cumulative breast cancer risk across PRS quintiles using a PRS constructed excluding BRCA1 and BRCA2 loci.

The script:
- defines the breast cancer survival dataset in women only
- constructs PRS quintiles based on the BRCA1/2-excluded PRS
- estimates the cumulative proportion of breast cancer cases observed by age 40 and age 50 within each PRS quintile

**Used for**
- Breast cancer age-specific risk analyses

**Notes**
- The analysis is restricted to individuals with evaluable age at event or censoring (up to age 80)
- The number of cases differs from the full dataset because:
  - events occurring after age 80 are excluded
  - individuals without age-at-onset information cannot be included
  
  
### `03e_breast_prs_noBRCA1_2_cumulative_risk_50_52y.R`

Estimates cumulative breast cancer risk at ages 50 and 52 using Cox proportional hazards models with age as the time scale and a standardized continuous PRS excluding BRCA1/2 loci.

The script:
- fits Cox models using continuous PRS
- derives model-based cumulative risk estimates
- extracts cumulative risk at age 50 and age 52
- provides estimates aligned with screening-relevant time points

**Used for**
- Breast cancer absolute risk estimation

**Notes**
- Analysis restricted to women
- Follow-up limited to age ≤ 80 years
- Estimates are based on model-predicted risk at average covariate values, with PRS fixed at its mean


### `03f_make_figure3B_BreastCancer_screening.R`

Generates Figure 3B, showing cumulative breast cancer risk across PRS quintiles and identifying age-specific risk equivalence.

The script:
- defines a breast cancer time-to-event dataset in women only
- fits Cox proportional hazards models using continuous standardized PRS
- generates model-based cumulative hazard curves for representative PRS values (median within each quintile)
- identifies the ages at which each PRS stratum reaches predefined cumulative risk thresholds


**Used for**
- Figure 3B

**Notes**
- Age is used as the time scale
- Analysis restricted to age ≤ 80 years
- PRS quintiles are derived from the standardized PRS distribution


### `03g_supplementary_table_S7_breast_prs_start_age_comparison.R`

Generates Supplementary Table S7 summarizing breast cancer risk across PRS strata using model-based cumulative hazard estimates.

The script:
- uses predicted cumulative hazard curves derived from Cox models (as in Figure 3B)
- defines PRS-specific starting ages based on a target cumulative risk threshold
- compares the number of cases observed by age 50 with those observed by the PRS-based starting age
- summarizes results across PRS quintiles

**Used for**
- Supplementary Table S7

**Notes**
- Analysis restricted to women
- Follow-up limited to age ≤ 80 years
- Starting ages are derived from model-based predicted curves
- Case counts are based on observed data

### `04a_chd_survival_prs_and_endophenotype_pgs.R`

Performs longitudinal survival analyses for incident coronary heart disease (CHD) using Cox proportional hazards models.

The script:
- evaluates the association between CHD-specific PRS and disease risk
- tests additional polygenic scores (PGS) for intermediate risk factors
- uses an incident-only design excluding baseline-prevalent CHD cases
- estimates hazard ratios, confidence intervals, and p-values for each score
- exports summary tables including case counts and model results

**Used for**
- CHD longitudinal risk analyses
- Evaluation of PRS and endophenotype-based PGS

**Notes**
- Baseline-prevalent CHD cases are excluded
- Secondary events are handled upstream to avoid inclusion as incident events
- Models are adjusted for age, sex, and genetic principal components (PC1–PC5)
- Event indicators are recoded to a binary format for consistency across analyses


### `04b_chd_univariate_lasso_minimal_models_and_forestplot.R`

Performs incident CHD survival analyses combining univariate, penalised, and multivariable modelling approaches.

The script:
- fits univariate Cox models for each PRS/PGS predictor
- applies LASSO penalised Cox regression for variable selection
- refits a multivariable Cox model including predictors retained by LASSO
- fits a minimal Cox model with selected predictors of interest
- generates forest plots comparing effect estimates across modelling strategies
- exports summary tables for all models

**Used for**
- CHD multivariable modelling and predictor selection
- Forest plot visualization of model results

**Notes**
- Predictors are pre-selected based on screening results from `04a_chd_survival_prs_and_endophenotype_pgs.R`
- Only predictors with P ≤ 0.001 are included (multiple-testing threshold)
- Analysis restricted to incident CHD events (baseline-prevalent cases excluded)
- Models are adjusted for age, sex, and genetic principal components (PC1–PC5)

### `04c_chd_prs_and_baseline_traits_models.R`

Fits Cox proportional hazards models for incident coronary heart disease (CHD) to compare genetic and phenotypic risk contributions.

The script:
- fits a traits-only model including selected polygenic predictors (PRS and PGS)
- fits a full model combining polygenic predictors with corresponding baseline phenotypic measurements
- evaluates whether baseline traits capture part of the genetic risk
- exports summary tables for both models

**Used for**
- CHD risk modelling combining genetic and baseline phenotypic factors

**Notes**
- Predictors include CHD PRS and selected polygenic scores (e.g. ApoB, SBP)
- Corresponding baseline traits (e.g. ApoB, SBP) are included in the full model
- Analysis restricted to incident CHD events (baseline-prevalent cases excluded)
- Models are adjusted for age, sex, and genetic principal components (PC1–PC5)


### `04d_chd_sex_interaction_and_sex_stratified_models.R`

Performs sex-specific analyses for incident coronary heart disease (CHD) to assess differences in genetic risk effects between women and men.

The script:
- fits Cox models including PRS × sex interaction terms
- performs sex-stratified Cox analyses in women and men separately
- estimates hazard ratios for PRS within each sex
- exports summary tables for interaction and stratified analyses
- reports case counts and sample sizes by sex

**Used for**
- Sex-specific CHD risk analyses

**Notes**
- Analysis restricted to incident CHD events (baseline-prevalent cases excluded)
- Interaction models include age, sex, PRS, and PRS × sex terms
- Sex-stratified models are adjusted for age and genetic principal components (PC1–PC5)


### `04e_chd_composite_score.R`

Builds and evaluates a composite polygenic score for incident coronary heart disease (CHD) based on multivariable Cox model coefficients.

The script:
- fits a multivariable Cox model including selected polygenic predictors
- derives a composite score using model beta coefficients
- standardizes and tests the composite score in a Cox model
- exports model coefficients and summary results

**Used for**
- Composite genetic risk modelling for CHD

**Notes**
- Composite score integrates CHD PRS and selected polygenic predictors (e.g. ApoB, SBP)
- Analysis restricted to incident CHD events (baseline-prevalent cases excluded)
- Secondary events are handled upstream to avoid inclusion as incident events


### `05a_compute_SCORE2_primary_prevention.R`

Computes SCORE2 cardiovascular risk estimates for primary prevention participants using the reference implementation from the RiskScorescvd package.

The script:
- selects eligible individuals based on clinical criteria (no prevalent CHD or stroke, no diabetes, age 40–69 years)
- prepares required inputs (age, sex, smoking status, blood pressure, cholesterol levels)
- converts cholesterol measurements from mg/dL to mmol/L
- calculates SCORE2 risk estimates and risk categories
- exports a table with SCORE2 scores and classifications

**Used for**
- Clinical risk estimation using SCORE2
- Comparison between genetic and clinical risk models

**Notes**
- Analysis restricted to QC-passed participants without baseline cardiovascular disease
- SCORE2 inputs follow standard definitions from the RiskScorescvd implementation
- Cholesterol values are converted to mmol/L prior to computation


### `05b_incident_chd_prediction_model_evaluation.R`

Evaluates predictive performance of incident CHD risk models combining genetic and clinical predictors.

The script:
- builds the analysis dataset for incident CHD prediction
- computes composite genetic risk scores
- fits multiple Cox models including genetic and clinical predictors
- evaluates model performance using:
  - Harrell’s C-index with bootstrap confidence intervals
  - change in C-index relative to the SCORE2-only model
  - time-dependent AUC at 10 years with bootstrap confidence intervals
- exports a summary table of model performance metrics

**Used for**
- Comparison of predictive performance across CHD risk models
- Evaluation of genetic vs clinical risk prediction

**Notes**
- Analysis restricted to primary prevention population (no baseline CHD or stroke, no diabetes, age 40–69 years)
- SCORE2 is also residualized on age and sex and evaluated in alternative models


### `05c_chd_model_performance_fig4_dca_calibration_auc.R`

Generates Figure 4 and supporting analyses summarizing predictive performance of CHD risk models.

The script:
- produces forest plots of 10-year AUC estimates
- performs decision curve analysis (DCA)
- generates calibration plots for uncalibrated and calibrated models
- computes additional summaries for model performance and calibration

It also:
- uses model evaluation results from `05b_incident_chd_prediction_model_evaluation.R`
- rebuilds the analysis dataset for decision curve and calibration analyses
- exports numerical summaries for main text and supplementary tables

**Used for**
- Figure 4
- Supplementary calibration and performance analyses

**Notes**
- Outcome: incident coronary heart disease (CHD)
- Integrates outputs from model evaluation and extends them with calibration and clinical utility analyses

### `05d_chd_cox_models_prs_score2_family_history.R`

Fits Cox proportional hazards models to evaluate the independent and combined associations of CHD PRS, SCORE2, and family history with incident coronary heart disease, adjusting for age, sex, and genetic principal components.

The script:
- fits sequential Cox models including PRS alone and in combination with SCORE2 and family history
- adjusts all models for Age, Sex, and PC1–PC5
- uses a common complete-case dataset across models to ensure direct comparability
- estimates hazard ratios and 95% confidence intervals for the main predictors
- extracts results only for PRS, SCORE2, and family history
- exports results formatted for Supplementary Table S12

**Used for**
- Supplementary Table S12
- Joint modelling of genetic, clinical, and familial risk

**Notes**
- Analysis restricted to the primary prevention population (no baseline CHD, no baseline stroke, no diabetes, age 40–69 years)
- PRS and SCORE2 are standardized to 1 standard deviation
- Family history is modelled as a binary predictor
- All reported estimates come from covariate-adjusted Cox models
- SCORE2 input is derived from `05a_compute_SCORE2_primary_prevention.R`

## Rare variant analyses

Rare variant gene-based association analyses were performed with **REGENIE** using a custom workflow maintained separately at:

`https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie`

This repository does **not** reimplement the full burden-testing pipeline. Instead, it documents the study-specific downstream steps used in this project, including:

- SLURM-based leave-one-variant-out (LOVO) analyses with REGENIE
- post-processing of LOVO output files
- annotation and filtering of selected rare variants
- extraction of variant dosages from sequence-imputed genotypes
- downstream association, survival, and interaction analyses for selected rare variants


### `scripts/slurm/regenie_lovo_submit.sh`

Submits **leave-one-variant-out (LOVO)** rare variant analyses with REGENIE step 2 on an HPC cluster using SLURM.

The script:
- reads a tab-delimited job table with one row per gene/mask/trait analysis
- builds the corresponding `--mask-lovo` argument for each test
- submits one SLURM job per row
- runs REGENIE step 2 with the same burden-test configuration used in the primary analysis

**Used for**
- LOVO analyses of exome-wide significant gene-based associations

**Notes**
- Runs downstream of the primary REGENIE gene-based association workflow
- Uses the REGENIE `--mask-lovo` option
- Requires a job table specifying gene, mask, threshold, phenotype file, covariate file, prediction file, and input paths
- Paths in the public version are placeholders and should be adapted locally


### `config/lovo_jobs_example.tsv`

Example tab-delimited input table used by `scripts/slurm/regenie_lovo_submit.sh`.

Each row defines one LOVO analysis and includes:
- gene identifier
- variant mask
- allele-frequency threshold
- trait name
- phenotype and covariate files
- REGENIE step 1 predictions
- input directories mounted inside the container

**Used for**
- Example configuration for LOVO job submission

**Notes**
- This file is provided as a template only
- Real paths and protected institutional file locations are not distributed in the public repository


### `06a_rare_variant_lovo_processing.R`

Processes REGENIE LOVO output files and quantifies the contribution of individual variants to significant gene-based association signals.

The script:
- reads and merges LOVO REGENIE output files across genes, masks, and traits
- extracts metadata from file names
- retains the burden test corresponding to the top gene-level signal from the reference results table
- computes a variance-explained proxy as `R2 = CHISQ / N`
- defines the reference result as the full gene-based burden test including all rare variants in the mask
- computes the relative contribution of each variant as:

`delta_r2 = (R2_full - R2_reduced) / R2_full`

- labels LOVO results according to technical validity
- identifies variants considered contributory when `delta_r2 >= 0.10`

**Used for**
- Post-processing of LOVO analyses
- Identification of candidate driver variants within significant gene-based signals

**Notes**
- LOVO results are first filtered for technical validity
- The public version uses placeholder paths and anonymized file structure
- The reference result corresponds to the full burden test without excluding any variant


### `06b_rare_variant_info_mac_filtering.R`

Annotates candidate rare variants with imputation quality and allele-frequency metrics derived from REGENIE summary statistics generated from the same sequence-imputed genotype dataset.

The script:
- matches selected variants to trait-specific REGENIE summary statistics files
- extracts:
  - `INFO`
  - `A1FREQ`
  - `N`
- computes minor allele frequency as:

`MAF = min(A1FREQ, 1 - A1FREQ)`

- computes minor allele count as:

`MAC = 2 * N * MAF`

- exports an annotated table with INFO, A1FREQ, MAF, and MAC

**Used for**
- Variant-level annotation of selected rare variants
- Documentation of upstream filtering criteria

**Notes**
- This script is intended to document how INFO and MAC were derived
- In the study workflow, the rare variant table used downstream is assumed to already satisfy the predefined INFO and MAC criteria


### `scripts/bash/extract_rare_variant_dosages.sh`

Extracts genotype dosages for selected rare variants from sequence-imputed genotype data using PLINK2.

The script:
- reads an imputed genotype dataset in BGEN format
- extracts a list of selected rare variants
- exports allele dosages using `PLINK2 --export A`

**Used for**
- Extraction of genotype dosages for selected rare variants

**Notes**
- In this dataset, the exported dosage corresponds to the **major allele**
- Rare-allele carrier status is therefore derived downstream by recoding the exported dosages
- The public version uses placeholder paths and anonymized inputs


### `06c_rare_variant_chd_association_analysis.R`

Performs downstream association analyses for selected rare variants in relation to coronary heart disease (CHD).

The script:
- loads PLINK2 genotype dosage output for selected rare variants
- recodes dosages from the major allele to rare-allele carrier dosage
- links variants to their imputed-genotype identifiers
- prepares the phenotype dataset
- performs cross-sectional logistic regression analyses for CHD
- performs survival analyses for incident CHD using Cox proportional hazards models
- restricts survival analyses to variants showing nominal evidence of association in the cross-sectional analysis (`P < 0.05`)

**Used for**
- Cross-sectional and longitudinal CHD analyses for selected rare variants

**Notes**
- The input rare variant table is assumed to already contain the final set of variants retained after upstream filtering
- Baseline-prevalent CHD cases are excluded from survival analyses
- Models are adjusted for age, sex, and genetic principal components (PC1–PC5)


### `06d_chd_prs_rare_variant_survival_figure_and_tables.R`

Generates the final CHD cumulative hazard figures combining PRS stratification and rare variant carrier status, and produces the corresponding supplementary tables.

The script:
- fits Cox models for:
  - PRS only
  - PRS plus carrier status for selected rare variants
- defines PRS quintiles:
  - on all participants in the PRS-only model
  - only among non-carriers in variant-specific models
- treats rare variant carriers as a separate group
- generates:
  - the main two-panel cumulative hazard figure (LPA and LDLR)
  - a supplementary cumulative hazard figure for SLC4A11
- exports model-derived hazard ratio tables for the supplementary material

**Used for**
- Figure 5 (LPA and LDLR panels)
- Supplementary Figure (SLC4A11)
- Supplementary survival model table

**Notes**
- Variant-specific PRS strata are defined among non-carriers only
- Carrier groups are shown separately
- The SLC4A11 analysis is presented as supplementary due to lack of independent replication
- Models are adjusted for age, sex, and PC1–PC5


### `06e_chd_prs_rare_variant_interaction_models.R`

Performs supplementary analyses to evaluate the relationship between CHD polygenic risk score (PRS) and rare variant carrier status in incident coronary heart disease (CHD).

The script:
- loads genotype dosages for selected rare variants exported with PLINK2
- recodes genotype dosages so that:
  - `0` = non-carrier of the rare allele
  - `1` = heterozygous carrier
  - `2` = homozygous carrier of the rare allele
- restricts the analysis to sequence-imputed samples
- prepares the CHD survival dataset
- standardizes the CHD PRS
- defines carrier status for the three selected rare variants:
  - `rs3798220` (`LPA`, Lp(a))
  - `rs58757394` (`SLC4A11`, SBP)
  - `rs730882080` (`LDLR`, LDL)

For each rare variant, the script fits:
- a linear model testing the association between PRS and rare variant carrier status
- a logistic model testing the association between carrier status and PRS
- Cox proportional hazards models for incident CHD including:
  - rare variant carrier status only
  - PRS only
  - rare variant carrier status and PRS jointly
  - rare variant carrier status × PRS interaction

The script also:
- compares additive and interaction Cox models using likelihood ratio tests
- quantifies changes in effect estimates before and after mutual adjustment
- exports formatted supplementary tables for:
  - PRS ~ carrier analyses
  - carrier ~ PRS analyses
  - Cox model comparison before and after mutual adjustment
  - PRS × rare variant interaction models

**Used for**
- Supplementary analyses of PRS and rare variant interplay in incident CHD
- Supplementary interaction tables

**Notes**
- Genotype dosages were exported with `PLINK2 --export A`
- In this dataset, exported dosages correspond to the major allele and are recoded downstream to define rare allele carrier status
- Analysis is restricted to sequence-imputed samples with non-missing CHD follow-up, PRS, and covariates
- Cox models are adjusted for age, sex, and genetic principal components (PC1–PC5)


### `06f_carrier_vs_noncarrier_phenotypes.R`

Generates a supplementary figure showing the distribution of quantitative phenotypes in carriers versus non-carriers of selected rare variants.

The script:
- loads recoded rare variant carrier data
- merges carrier status with phenotype data
- generates violin plots comparing phenotype distributions in carriers and non-carriers
- reports descriptive statistics for each variant–phenotype pair

**Used for**
- Supplementary figure on phenotype distributions by rare variant carrier status

**Notes**
- Carrier status is defined from recoded dosages, where `1` or `2` indicates carriage of the rare allele
- Phenotypes shown correspond to traits linked to the selected rare variants