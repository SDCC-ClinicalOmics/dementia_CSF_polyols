# Systems medicine approach identifies CSF polyol dysregulation across the dementia continuums

In this repository, we have submitted the code necessary to replicate our manuscript titled "Systems medicine approach identifies CSF polyol dysregulation across the dementia continuums" with our datasets. The repository is organized as follows:

/code
    - ../metabolomics_evaluation: QC, EDA and univariate statistics RMarkdowns for both CSF and plasma tissues targeted and untargeted metabolomics. For each tissue-targeted/untargeted GC/MS combination, we include:
        - 00_tissue_targeted/untargeted_QC.Rmd: includes general data quality evaluation, peak-filtering, blank-signal removal, missingness and duplications evaluation and data normalization assessment.
        - 01_tissue_targeted/untargeted_EDA.Rmd: exploratory data analysis, with PCA analyses to identify outliers and potential confounding factors.
        - 03A_bayesian_modelling.R: R script dedicated to the Bayesian models run for the binary clinical groups contrasts.
        - 03_tissue_targeted/untargeted_univariate.Rmd: includes all tested univariate statistical approaches for both continuous clinical biomarkers and bayesian modelling results parsing.
    - ../modelling: directory with the Python scripts and Quarto notebooks for the ElasticNet regression modelling and Structural Equation Modelling (SEM) approaches. Script names are self-explanatory of what the script includes.
        - 000_dementia_ElasticNet_modelruns.py
        - 01_demenntiabiomarkers_elasticnet_discovery.qmd
        - 001_dementia_ElasticNEt_modelbootstrapping.py
        - 02_dementiabiomarkers_elasticnet_bootstrapping.qmd
        - 002_dementia_ElasticNet_validation.py
        - 03_dementiabiomarkers_elasticnet_validation.qmd
        - 003_dementia_SEM_mediation.py
        - 04_dementiabiomarkers_SEM_mediation.qmd
        - 004_SNPs_extraction.py
        - 004_SNPs_extraction.sh
        - 005_dementia_SEM_validation_genetics.py
        - 05_dementiabiomarkers_SEM_mediation_validation.qmd
    - metabolomics_validation_CruchagaCohort.R: script for the validation of both our metabolomics and proteomics findings in the Cruchaga (ADNI) cohort
    - polyols_focused_analysis.R: main script for the manuscript CSF polyols-centered analysis and figures generation.