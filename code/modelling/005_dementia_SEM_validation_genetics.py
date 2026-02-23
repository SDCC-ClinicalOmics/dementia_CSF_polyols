################################################################
####### Python script for Structural Equation Modelling ########
########### validation and genetics effects testing ############
################################################################

# Consideration notes
# We need to scale all numerical elements included in the model
# due to convergence failing. Since we are building latent variables,
# scaling everything to the same range will make sense.

# In SEM, the optimizer (Maximum Likelihood) works by minimizing the 
# difference between the observed covariance matrix and the model-implied 
# covariance matrix. If Age has a variance of 100 (Age≈75) and your
# CRMN-transformed metabolites have a variance of 0.01, 
# the gradients become "unbalanced".

# import libraries
import os as os
import numpy as np
import pandas as pd
from semopy import Model
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

import pickle

# ----------------------------------
# Setup
# ----------------------------------
OUTPUT_DIR = "../results/ADNI_genetics"
os.makedirs(OUTPUT_DIR, exist_ok=True)
np.random.seed(2548)

# ----------------------------------
# Dataset 1: Full metabolomics (for pure replication)
# ----------------------------------
df_full = pd.read_csv("../results/ADNI_genetics/Cruchaga_metabolomics_relevant.txt", sep = "\t", encoding = 'latin1')

# Encode diagnosis
df_full['DIAGNOSIS'] = df_full['DIAGNOSIS'].map({1: 'HC', 2: 'MCI', 3: 'AD'}) # we recode to 0/1 to clarify this is a binary variable

# Encode gender
df_full['Gender_encoded'] = df_full['PTGENDER'].map({1: 0, 2: 1}) # we recode to 0/1 to clarify this is a binary variable

# Transform biomarkers (same as discovery)
df_full['TAU'] = StandardScaler().fit_transform(
    np.log(df_full['TAU']).values.reshape(-1,1)
).ravel()

df_full['PTAU'] = StandardScaler().fit_transform(
    np.log(df_full['PTAU']).values.reshape(-1,1)
).ravel()

df_full['ABETA42'] = StandardScaler().fit_transform(
    df_full['ABETA42'].values.reshape(-1,1)
).ravel()

# ----------------------------------
# Dataset 2: Metabolomics + Genetics (for genetic adjustment)
# ----------------------------------
df_genetics = pd.read_csv("../results/ADNI_genetics/Cruchaga_metabolomics_relelevant_WITHGENETICS.txt", sep = "\t", encoding = 'latin1')

# Encode diagnosis
df_genetics['DIAGNOSIS'] = df_genetics['DIAGNOSIS'].map({1: 'HC', 2: 'MCI', 3: 'AD'}) # we recode to 0/1 to clarify this is a binary variable

# Encode gender
df_genetics['Gender_encoded'] = df_genetics['PTGENDER'].map({1: 0, 2: 1}) # we recode to 0/1 to clarify this is a binary variable

# Transform biomarkers
df_genetics['TAU'] = StandardScaler().fit_transform(
    np.log(df_genetics['TAU']).values.reshape(-1,1)
).ravel()

df_genetics['PTAU'] = StandardScaler().fit_transform(
    np.log(df_genetics['PTAU']).values.reshape(-1,1)
).ravel()

df_genetics['ABETA42'] = StandardScaler().fit_transform(
    df_genetics['ABETA42'].values.reshape(-1,1)
).ravel()

# ----------------------------------
# Define metabolite configurations
# ----------------------------------
# ADJUST these metabolite names to match your validation cohort columns
VALIDATION_CONFIG = {
    'TAU': [
        'arabitol.xylitol', 'threonate', 'ribonate', 'glucose', 'gluc.ery', 
        # Add more based on your actual column names
    ],
    'PTAU': [
        'arabitol.xylitol', 'threonate', 'ribonate', 'glucose'
    ],
    'ABETA42': ['threonate', 'mannitol.sorbitol'
    ]
}


# ----------------------------------
# Define necessary functions
# ----------------------------------
def run_mediation_validation(df_sub, mediator, mets, contrast_name, 
                             genetic_model='none', dataset_name='full',
                             n_iter=1000):
    """
    Run SEM mediation in validation cohort
    
    Parameters:
    -----------
    df_sub : DataFrame
        Subset of data for this contrast
    mediator : str
        Biomarker (TAU, PTAU, ABETA42)
    mets : list
        Metabolite names
    contrast_name : str
        e.g., "HC_vs_AD"
    genetic_model : str
        'none' - No genetic adjustment (uses full dataset)
        'adjusted' - Add GRS_unw as covariate (requires genetics dataset)
        'pathway' - GRS -> Metabolites (requires genetics dataset)
    dataset_name : str
        'full' or 'genetics' - for tracking which dataset was used
    n_iter : int
        Bootstrap iterations
    """
    
    print(f"\n{'='*70}")
    print(f"Contrast: {contrast_name} | Mediator: {mediator}")
    print(f"Genetic Model: {genetic_model} | Dataset: {dataset_name}")
    print(f"{'='*70}")
    
    # Validation checks
    if genetic_model != 'none' and 'GRS_unw' not in df_sub.columns:
        raise ValueError(f"Genetic model '{genetic_model}' requires GRS_unw column. Use genetics dataset.")
    
    print(f"Sample size: {len(df_sub)}")
    print(f"Metabolites: {len(mets)}")
    
    # Scale variables
    cols_to_scale = mets + [mediator, 'AGE']
    
    if genetic_model != 'none':
        cols_to_scale.append('GRS_unw')
    
    scaler = StandardScaler()
    df_sub[cols_to_scale] = scaler.fit_transform(df_sub[cols_to_scale])
    
    # Clean metabolite names
    met_map = {col: f"m{i}" for i, col in enumerate(mets)}
    inv_met_map = {v: k for k, v in met_map.items()}
    df_local = df_sub.copy().rename(columns=met_map)
    m_names = list(met_map.values())
    
    # Build model specification
    if genetic_model == 'none':
        model_spec = f"""
            SugarAlcohol =~ {" + ".join(m_names)}
            {mediator} ~ SugarAlcohol + AGE + Gender_encoded
            outcome ~ {mediator} + SugarAlcohol + AGE + Gender_encoded
        """
        print("Model: Replication (no genetic adjustment)")
    
    elif genetic_model == 'adjusted':
        model_spec = f"""
            SugarAlcohol =~ {" + ".join(m_names)}
            {mediator} ~ SugarAlcohol + AGE + Gender_encoded + GRS_unw
            outcome ~ {mediator} + SugarAlcohol + AGE + Gender_encoded + GRS_unw
        """
        print("Model: Genetic adjustment (GRS as covariate)")
    
    elif genetic_model == 'pathway':
        model_spec = f"""
            SugarAlcohol =~ {" + ".join(m_names)}
            SugarAlcohol ~ GRS_unw + AGE + Gender_encoded
            {mediator} ~ SugarAlcohol + AGE + Gender_encoded
            outcome ~ {mediator} + SugarAlcohol + AGE + Gender_encoded + GRS_unw
        """
        print("Model: Genetic pathway (GRS -> Metabolites)")
    
    # Bootstrap
    boots = {'ind': [], 'dir': [], 'tot': []}
    loadings_list = []
    success_count = 0
    
    print(f"\nBootstrapping {n_iter} iterations...")
    
    for i in range(n_iter):
        if i % 200 == 0 and i > 0:
            print(f"  Progress: {i}/{n_iter} ({success_count} converged)")
        
        sample = df_local.sample(n=len(df_local), replace=True)
        
        try:
            m = Model(model_spec)
            m.fit(sample)
            res = m.inspect()
            
            # Extract paths
            a = res.loc[(res['lval'] == mediator) & (res['rval'] == 'SugarAlcohol'), 'Estimate'].values[0]
            b = res.loc[(res['lval'] == 'outcome') & (res['rval'] == mediator), 'Estimate'].values[0]
            c_p = res.loc[(res['lval'] == 'outcome') & (res['rval'] == 'SugarAlcohol'), 'Estimate'].values[0]
            
            # Loadings
            it_loadings = res.loc[(res['op'] == '~') & (res['rval'] == 'SugarAlcohol') & (res['lval'].isin(m_names))]
            if len(it_loadings) > 0:
                loadings_list.append(dict(zip(it_loadings['lval'], it_loadings['Estimate'])))
            
            boots['ind'].append(a * b)
            boots['dir'].append(c_p)
            boots['tot'].append((a * b) + c_p)
            success_count += 1
            
        except:
            continue
    
    convergence_rate = (success_count / n_iter) * 100
    print(f"\n✓ Convergence: {convergence_rate:.1f}% ({success_count}/{n_iter})")
    
    if convergence_rate < 70:
        print(f"  ⚠ WARNING: Low convergence rate!")
    
    # Calculate statistics
    def get_stats(boot_list):
        if not boot_list or len(boot_list) == 0:
            return [np.nan] * 4
        median = np.percentile(boot_list, 50)
        lower = np.percentile(boot_list, 2.5)
        upper = np.percentile(boot_list, 97.5)
        p = 2 * min(np.mean(np.array(boot_list) > 0), np.mean(np.array(boot_list) < 0))
        return median, lower, upper, p
    
    ind_s = get_stats(boots['ind'])
    dir_s = get_stats(boots['dir'])
    tot_s = get_stats(boots['tot'])
    
    # Proportion mediated
    if np.sign(ind_s[0]) != np.sign(dir_s[0]):
        prop_label = "Ratio |Ind/Dir|"
        prop_val = abs(ind_s[0] / dir_s[0]) if dir_s[0] != 0 else np.nan
    else:
        prop_label = "Prop. Mediated"
        prop_val = (ind_s[0] / tot_s[0]) if tot_s[0] != 0 else np.nan
    
   
    fname = f"{contrast_name}_{mediator}_{genetic_model}_{dataset_name}"
    
    # Summary
    summary = pd.DataFrame({
        "Contrast": contrast_name,
        "Mediator": mediator,
        "Genetic_Model": genetic_model,
        "Dataset": dataset_name,
        "Sample_Size": len(df_sub),
        "N_Metabolites": len(mets),
        "Effect": ["Indirect", "Direct", "Total", prop_label],
        "Estimate": [ind_s[0], dir_s[0], tot_s[0], prop_val],
        "CI_Lower": [ind_s[1], dir_s[1], tot_s[1], np.nan],
        "CI_Upper": [ind_s[2], dir_s[2], tot_s[2], np.nan],
        "p_val": [ind_s[3], dir_s[3], tot_s[3], np.nan],
        "Convergence_Rate": [convergence_rate] * 4
    })
    
    summary.to_csv(f"{OUTPUT_DIR}/SEM_validation_{fname}.csv", index=False)
    
    # Bootstrap distributions
    with open(f"{OUTPUT_DIR}/SEM_validation_boot_{fname}.pkl", 'wb') as f:
        pickle.dump(boots, f)
    
    # Loadings
    print(f"Collected loadings from {len(loadings_list)} iterations")

    if loadings_list:
        loadings_df_raw = pd.DataFrame(loadings_list)

        loadings_summary = []
        for m_code in m_names:
            m_real_name = inv_met_map[m_code]
            m_data = loadings_df_raw[m_code].dropna()

            if len(m_data) == 0:
                continue

            med, low, upp, p = get_stats(m_data.tolist())

            loadings_summary.append({
                "Metabolite": m_real_name,
                "Median_Loading": med,
                "CI_Lower": low,
                "CI_Upper": upp,
                "p_val": p
            })

        load_df = pd.DataFrame(loadings_summary)

        if len(load_df) > 0:
            _, load_df['fdr_p'], _, _ = multipletests(
                load_df['p_val'], method='fdr_bh'
            )

            load_df.to_csv(
                f"{OUTPUT_DIR}/SEM_validation_loadings_{fname}.csv",
                index=False
            )
    
    # Print summary
    print(f"\n{'='*70}")
    print("RESULTS")
    print(f"{'='*70}")
    print(f"Indirect: {ind_s[0]:.4f} [{ind_s[1]:.4f}, {ind_s[2]:.4f}], p={ind_s[3]:.4f}")
    print(f"Direct:   {dir_s[0]:.4f} [{dir_s[1]:.4f}, {dir_s[2]:.4f}], p={dir_s[3]:.4f}")
    print(f"Total:    {tot_s[0]:.4f} [{tot_s[1]:.4f}, {tot_s[2]:.4f}], p={tot_s[3]:.4f}")
    print(f"{prop_label}: {prop_val:.4f}")
    
    return summary


def compare_datasets_and_genetics(contrast, mediator):
    """
    Compare results across:
    1. Full dataset (replication)
    2. Genetics subset (replication) 
    3. Genetics subset (with GRS adjustment)
    """

    print(f"\n{'='*70}")
    print(f"COMPREHENSIVE COMPARISON: {contrast} - {mediator}")
    print(f"{'='*70}")
    
    analyses = [
        ('none', 'full', 'Full Dataset Replication'),
        ('none', 'genetics', 'Genetics Subset Replication'),
        ('adjusted', 'genetics', 'Genetic Adjustment (GRS)'),
        ('pathway', 'genetics', 'Genetic Pathway Model')
    ]
    
    results = []
    
    for model, dataset, label in analyses:
        try:
            fname = f"../results/ADNI_genetics/SEM_validation_{contrast}_{mediator}_{model}_{dataset}.csv"
            df = pd.read_csv(fname)
            ind = df[df['Effect'] == 'Indirect'].iloc[0]
            
            results.append({
                'Analysis': label,
                'Dataset': dataset,
                'N': int(ind['Sample_Size']),
                'Indirect': ind['Estimate'],
                'CI_Lower': ind['CI_Lower'],
                'CI_Upper': ind['CI_Upper'],
                'p_value': ind['p_val'],
                'Sig': '***' if ind['p_val'] < 0.001 else '**' if ind['p_val'] < 0.01 else '*' if ind['p_val'] < 0.05 else 'ns',
                'Conv%': ind['Convergence_Rate']
            })
        except Exception as e:
            print(f"  ⚠ Missing: {label}")
    
    if len(results) >= 3:
        comp_df = pd.DataFrame(results)
        print("\n" + comp_df.to_string(index=False))
        
        # Compare full vs genetics subset
        print(f"\n{'='*70}")
        print("DATASET COMPARISON")
        print(f"{'='*70}")
        full_effect = results[0]['Indirect']
        subset_effect = results[1]['Indirect']
        print(f"Full dataset (n={results[0]['N']}):     {full_effect:.4f}")
        print(f"Genetics subset (n={results[1]['N']}): {subset_effect:.4f}")
        
        if abs(full_effect - subset_effect) / abs(full_effect) < 0.2:
            print("→ Genetics subset is REPRESENTATIVE of full cohort")
        else:
            print("→ ⚠ Genetics subset may differ from full cohort")
        
        # Genetic attenuation
        print(f"\n{'='*70}")
        print("GENETIC ATTENUATION")
        print(f"{'='*70}")
        base = results[1]['Indirect']  # Genetics subset, no adjustment
        adj = results[2]['Indirect']   # Genetics subset, with GRS
        
        if base != 0:
            atten = ((base - adj) / base) * 100
            print(f"Before GRS adjustment: {base:.4f}")
            print(f"After GRS adjustment:  {adj:.4f}")
            print(f"Attenuation:           {atten:.1f}%")
            
            if abs(atten) < 25:
                print("\n→ Effect PERSISTS after genetic adjustment")
                print("  Metabolic/environmental pathway independent of genetics")
            elif abs(atten) > 75:
                print("\n→ Effect LARGELY EXPLAINED by genetics")
                print("  Genetic confounding likely")
            else:
                print("\n→ PARTIAL genetic contribution")
        
        return comp_df
    
    return None

# ----------------------------------
# Define variables for run
# ----------------------------------
# Define contrasts
CONTRASTS = [
    ("HC", "MCI"),
    ("HC", "AD"),
    ("MCI", "AD")
]

# Storage
all_results = []

# ----------------------------------
# Main analysis loop
# ----------------------------------
for g1, g2 in CONTRASTS:
    contrast_name = f"{g1}_vs_{g2}"
    
    print(f"\n\n{'#'*70}")
    print(f"# CONTRAST: {contrast_name}")
    print(f"{'#'*70}")
    
    # Subset both datasets
    df_full_subset = df_full[df_full['DIAGNOSIS'].isin([g1, g2])].copy()
    df_full_subset['outcome'] = (df_full_subset['DIAGNOSIS'] == g2).astype(int)
    
    df_gen_subset = df_genetics[df_genetics['DIAGNOSIS'].isin([g1, g2])].copy()
    df_gen_subset['outcome'] = (df_gen_subset['DIAGNOSIS'] == g2).astype(int)
    
    print(f"\nFull dataset: {g1}={sum(df_full_subset['outcome']==0)}, {g2}={sum(df_full_subset['outcome']==1)}")
    print(f"Genetics dataset: {g1}={sum(df_gen_subset['outcome']==0)}, {g2}={sum(df_gen_subset['outcome']==1)}")
    
    for mediator, metabolites in VALIDATION_CONFIG.items():
        
        if len(metabolites) == 0:
            print(f"\n⚠ Skipping {mediator} - no metabolites")
            continue
        
        print(f"\n\n{'*'*70}")
        print(f"* MEDIATOR: {mediator}")
        print(f"* Metabolites: {metabolites}")
        print(f"{'*'*70}")
        
        # ---------------------------------------------------------------
        # Analysis 1: Replication with FULL dataset (maximum power)
        # ---------------------------------------------------------------
        print("\n" + "-"*70)
        print("ANALYSIS 1: REPLICATION (Full metabolomics dataset, no genetics)")
        print("-"*70)
        
        res1 = run_mediation_validation(
            df_full_subset, mediator, metabolites,
            contrast_name, genetic_model='none', dataset_name='full',
            n_iter=1000
        )
        all_results.append(res1)
        
        # ---------------------------------------------------------------
        # Analysis 2: Replication with GENETICS dataset (for comparison)
        # ---------------------------------------------------------------
        print("\n" + "-"*70)
        print("ANALYSIS 2: REPLICATION (Genetics subset, no genetics adjustment)")
        print("(To ensure genetic subset is representative)")
        print("-"*70)
        
        res2 = run_mediation_validation(
            df_gen_subset, mediator, metabolites,
            contrast_name, genetic_model='none', dataset_name='genetics',
            n_iter=1000
        )
        all_results.append(res2)
        
        # ---------------------------------------------------------------
        # Analysis 3: Genetic adjustment (GENETICS dataset only)
        # ---------------------------------------------------------------
        print("\n" + "-"*70)
        print("ANALYSIS 3: GENETIC ADJUSTMENT (Genetics dataset with GRS)")
        print("-"*70)
        
        res3 = run_mediation_validation(
            df_gen_subset, mediator, metabolites,
            contrast_name, genetic_model='adjusted', dataset_name='genetics',
            n_iter=1000
        )
        all_results.append(res3)
        
        # ---------------------------------------------------------------
        # Analysis 4: Genetic pathway (GENETICS dataset only)
        # ---------------------------------------------------------------
        print("\n" + "-"*70)
        print("ANALYSIS 4: GENETIC PATHWAY (Genetics dataset, GRS -> Metabolites)")
        print("-"*70)
        
        res4 = run_mediation_validation(
            df_gen_subset, mediator, metabolites,
            contrast_name, genetic_model='pathway', dataset_name='genetics',
            n_iter=1000
        )
        all_results.append(res4)

# Save all results
all_results_df = pd.concat(all_results, ignore_index=True)
all_results_df.to_csv('../results/ADNI_genetics/all_validation_results.csv', index=False)

print(f"\n\n{'='*70}")
print("✓ ALL ANALYSES COMPLETE")
print(f"{'='*70}")




# Run comparisons
compare_datasets_and_genetics("HC_vs_AD", "TAU")
compare_datasets_and_genetics("HC_vs_AD", "PTAU")
compare_datasets_and_genetics("HC_vs_AD", "ABETA42")