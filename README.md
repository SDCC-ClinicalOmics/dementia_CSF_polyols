# Systems Medicine Approach Identifies CSF Polyol Dysregulation Across the Dementia Continuums

## Overview

This repository contains the code necessary to replicate the findings from our manuscript titled **"Systems medicine approach identifies CSF polyol dysregulation across the dementia continuums"** along with the associated datasets.

## Repository Structure

### `/code` — Analysis and Modeling Scripts

#### `metabolomics_evaluation/`
Quality control, exploratory data analysis, and univariate statistics for CSF and plasma targeted and untargeted metabolomics.

For each tissue–targeted/untargeted combination, we provide:

- **`00_tissue_targeted/untargeted_QC.Rmd`** — Data quality evaluation including:
  - Peak filtering and blank-signal removal
  - Missingness and duplication assessment
  - Data normalization evaluation

- **`01_tissue_targeted/untargeted_EDA.Rmd`** — Exploratory data analysis including:
  - PCA analyses to identify outliers
  - Potential confounding factor identification

- **`03_tissue_targeted/untargeted_Univariate.Rmd`** — Statistical testing:
  - Univariate statistical approaches for continuous clinical biomarkers
  - Bayesian modeling results

- **`03A_bayesian_modelling.R`** — R script for Bayesian models applied to binary clinical group contrasts

#### `modelling/`
Python scripts and Quarto notebooks for ElasticNet regression modeling and Structural Equation Modeling (SEM) approaches:

- **`000_dementia_ElasticNet_modelruns.py`** — ElasticNet model discovery
- **`01_dementiabiomarkers_elasticnet_discovery.qmd`** — Discovery analysis notebook
- **`001_dementia_ElasticNet_modelbootstrapping.py`** — Model bootstrapping
- **`02_dementiabiomarkers_elasticnet_bootstraping.qmd`** — Bootstrapping notebook
- **`002_dementia_ElasticNet_validation.py`** — Model validation
- **`03_dementiabiomarkers_elasticnet_validation.qmd`** — Validation notebook
- **`003_dementia_SEM_mediation.py`** — SEM mediation analysis
- **`04_dementiabiomarkers_SEM_mediation.qmd`** — SEM mediation notebook
- **`005_dementia_SEM_validation_genetics.py`** — Genetic validation of SEM results
- **`05_dementiabiomarkers_SEM_mediation_validation.qmd`** — Validation notebook
- **`004_SNPs_extraction.py`** & **`004_SNPs_extraction.sh`** — SNP extraction pipeline

#### Root Analysis Scripts

- **`metabolomics_validation_CruchagaCohort.R`** — Validation of metabolomics and proteomics findings in the Cruchaga (ADNI) cohort
- **`polyols_focused_analysis.R`** — Main manuscript analysis and figure generation focused on CSF polyols

### `/data` — Datasets

Raw and processed data files organized by sample type and external cohorts

### `/results` — Output Directory

Analysis results and generated outputs