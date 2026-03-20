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
- annotates curves with age labels at these threshold intersections

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