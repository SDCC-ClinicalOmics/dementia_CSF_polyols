################################################################
####### Python script for Structural Equation Modelling ########
################################################################

# Consideration notes
# We need to scale all numerical elements included in the model
# due to convergence failing. Since we are building latent variables,
# scaling everything to the same range will make sense.

# In SEM, the optimizer (Maximum Likelihood) works by minimizing the 
# difference between the observed covariance matrix and the model-implied 
# covariance matrix. If Age has a variance of 100 (Age≈75) and your
# CRMN-transformed metabolites have a variance of 0.01, 
# the gradients become "unbalanced.

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
OUTPUT_DIR = "../results/modelling"
os.makedirs(OUTPUT_DIR, exist_ok=True)
np.random.seed(2548)

# load data
df = pd.read_csv('../data/csf_ratios2.csv', sep=';', encoding='latin1')
## encode gender
df['Gender_encoded'] = (df['Gender'] == 'Male').astype(int)

## transform & scale tTau, pTay and AB42 as we did before for consistency
df['TAU'] =  StandardScaler().fit_transform(np.log(df['TAU']).values.reshape(-1,1)).ravel() 
df['PTAU'] = StandardScaler().fit_transform(np.log(df['PTAU']).values.reshape(-1,1)).ravel()   
df['ABETA42'] = StandardScaler().fit_transform(df['ABETA42'].values.reshape(-1,1)).ravel()

# Define metabolites to use in each biomarker (based on ElasticNetResults)
STUDY_CONFIG = {
    'TAU': ['L.....Arabitol', 'd.Galactose', 'D.Threitol', 'Ribonic.acid', 
    'd.Glucose', 'gluc.erythitol', 'L.Threonic.acid', 'X3.Deoxy.erythro.pentitol'],
    'PTAU': ['L.....Arabitol', 'd.Galactose', 'D.Threitol', 'Ribonic.acid', 
    'd.Glucose','L.Threonic.acid'],
    'ABETA42': ['L.Threonic.acid', 'D.....Ribofuranose', 'D.Threitol', 'D.....Xylose', 'Sorbitol']
}

# Define clinical contrasts to run
CONTRASTS = [
    ("HC", "MCI"),
    ("HC", "AD"), 
    ("HC", "VaD"),
    ("MCI", "AD"),
    ("MCI", "VaD"),
    ("AD", "VaD")
]

# -----------------------------------------
# Functions
# -----------------------------------------
def get_stats(boot_list):
    if not boot_list or len(boot_list) == 0:
        return [np.nan] * 4
    median = np.percentile(boot_list, 50)
    lower = np.percentile(boot_list, 2.5)
    upper = np.percentile(boot_list, 97.5)
    p = 2 * min(np.mean(np.array(boot_list) > 0), np.mean(np.array(boot_list) < 0))
    return median, lower, upper, p

def run_mediation(df_sub, mediator, mets, contrast_name, n_iter=1000):
    print(f"Running: {contrast_name} | Mediator: {mediator}")
    #
    # Scale everything that will be included in the latent variable
    cols_to_scale = mets + [mediator, 'AgeAtVisit']
    scaler = StandardScaler()
    df_sub[cols_to_scale] = scaler.fit_transform(df_sub[cols_to_scale])
    # Clean metabolite names for SEM syntax
    met_map = {col: f"m{i}" for i, col in enumerate(mets)}
    inv_met_map = {v: k for k, v in met_map.items()} # To map back to real names
    df_local = df_sub.copy().rename(columns=met_map)
    m_names = list(met_map.values())
    # try simplifying the model
    #sample['SugarAlcohol_Score'] = sample[m_names].mean(axis=1)
    #sample['SugarAlcohol_Score'] = StandardScaler().fit_transform(sample[['SugarAlcohol_Score']])
    model_spec = f"""
        SugarAlcohol =~ {" + ".join(m_names)}
        {mediator} ~ SugarAlcohol + AgeAtVisit + Gender_encoded
        outcome ~ {mediator} + SugarAlcohol + AgeAtVisit + Gender_encoded
    """
    # 2. Simplified Model (Path Analysis)
    # model_spec = f"""
    #    {mediator} ~ SugarAlcohol_Score + AgeAtVisit + Gender_encoded
    #    outcome ~ {mediator} + SugarAlcohol_Score + AgeAtVisit + Gender_encoded
    #"""
    #
    boots = {'ind': [], 'dir': [], 'tot': []}
    loadings_list = [] # list for loadinfgs store
    success_count = 0  # Track successful iterations
    #
    for _ in range(n_iter):
        sample = df_local.sample(n=len(df_local), replace=True)
        try:
            m = Model(model_spec)
            m.fit(sample)
            res = m.inspect()
            #       
            # Extract paths
            a = res.loc[(res['lval'] == mediator) & (res['rval'] == 'SugarAlcohol'), 'Estimate'].values[0]
            b = res.loc[(res['lval'] == 'outcome') & (res['rval'] == mediator), 'Estimate'].values[0]
            c_p = res.loc[(res['lval'] == 'outcome') & (res['rval'] == 'SugarAlcohol'), 'Estimate'].values[0]
            #
            # PINPOINT: Capture Loadings
            it_loadings = res.loc[(res['op'] == '~') & (res['rval'] == 'SugarAlcohol') & (res['lval'].isin(m_names))]
            loadings_list.append(dict(zip(it_loadings['lval'], it_loadings['Estimate'])))
            #
            boots['ind'].append(a * b)
            boots['dir'].append(c_p)
            boots['tot'].append((a * b) + c_p)
            success_count += 1
        except:
            continue
    #        
    print(f"   --> Convergence Rate: {(success_count/n_iter)*100:.1f}% ({success_count}/{n_iter})")    # Calculate Stats

    # back-transform names for loadings data
    if loadings_list:
        loadings_df_raw = pd.DataFrame(loadings_list)
        
        # Calculate median and CIs for each metabolite loading
        loadings_summary = []
        for m_code in m_names:
            m_real_name = inv_met_map[m_code]
            m_data = loadings_df_raw[m_code].dropna()
            
            median, lower, upper, p = get_stats(m_data.tolist())
            
            loadings_summary.append({
                "Metabolite": m_real_name,
                "Median_Loading": median,
                "CI_Lower": lower,
                "CI_Upper": upper,
                "p_val": p
            })
            
        load_df = pd.DataFrame(loadings_summary)
        # Apply FDR correction to the loadings
        _, load_df['fdr_p'], _, _ = multipletests(load_df['p_val'], method='fdr_bh')
        
        load_fname = f"SEM_mediation_loadings_{contrast_name}_{mediator}".replace(" ", "_")
        load_df.to_csv(f"{OUTPUT_DIR}/{load_fname}.csv", index=False)

    ind_s = get_stats(boots['ind'])
    dir_s = get_stats(boots['dir'])
    tot_s = get_stats(boots['tot'])
    #
    # Determine Ratio/Prop
    if np.sign(ind_s[0]) != np.sign(dir_s[0]):
        prop_label, prop_val = "Ratio |Ind/Dir|", abs(ind_s[0] / dir_s[0])
    else:
        prop_label, prop_val = "Prop. Mediated", (ind_s[0] / tot_s[0])
    #
    # Save CSV Summary
    summary = pd.DataFrame({
        "Contrast": contrast_name,
        "Mediator": mediator,
        "Effect": ["Indirect", "Direct", "Total", prop_label],
        "Estimate": [ind_s[0], dir_s[0], tot_s[0], prop_val],
        "CI_Lower": [ind_s[1], dir_s[1], tot_s[1], np.nan],
        "CI_Upper": [ind_s[2], dir_s[2], tot_s[2], np.nan],
        "p_val": [ind_s[3], dir_s[3], tot_s[3], np.nan]
    })
    #
    # Save Pickle for Quarto
    fname = f"{contrast_name}_{mediator}".replace(" ", "_")
    summary.to_csv(f"{OUTPUT_DIR}/SEM_mediation_table_{fname}.csv", index=False)
    with open(f"{OUTPUT_DIR}/SEM_mediation_boot_{fname}.pkl", 'wb') as f:
        pickle.dump(boots, f)

    # ROC Curves computation
    try:
        m_final = Model(model_spec)
        m_final.fit(df_local)
        scores = m_final.predict(df_local)
        
        # Check if semopy actually produced the latent score
        if 'SugarAlcohol' in scores.columns:
            df_local['Metabolic_Score'] = scores['SugarAlcohol']
        else:
            print(f"   --- Using weighted average fallback for {contrast_name}")
            # Get the median loadings we just calculated to create a manual score
            # This is mathematically very close to the latent factor
            weights = load_df.set_index('Metabolite')['Median_Loading'].to_dict()
            # Map weights back to 'm0', 'm1' names
            m_weights = {met_map[k]: v for k, v in weights.items()}
            
            df_local['Metabolic_Score'] = sum(df_local[m] * w for m, w in m_weights.items())

        models_to_run = {
            'Biomarker': [mediator, 'AgeAtVisit', 'Gender_encoded'],
            'Metabolites': ['Metabolic_Score', 'AgeAtVisit', 'Gender_encoded'],
            'Combined': [mediator, 'Metabolic_Score', 'AgeAtVisit', 'Gender_encoded']
        }

        roc_results = {}
        for name, features in models_to_run.items():
            X, y = df_local[features], df_local['outcome']
            lr = LogisticRegression(max_iter=1000).fit(X, y)
            probs = lr.predict_proba(X)[:, 1]
            fpr, tpr, _ = roc_curve(y, probs)
            roc_results[name] = {'fpr': fpr.tolist(), 'tpr': tpr.tolist(), 'auc': auc(fpr, tpr)}
        
        # Save ROC data
        with open(f"{OUTPUT_DIR}/SEM_mediation_roc_{contrast_name}_{mediator}.pkl", 'wb') as f:
            pickle.dump(roc_results, f)
    except Exception as e:
        print(f"   !!! ROC calculation failed: {e}")

    # 1. Clean names for SEM syntax inside this function too
    met_map = {col: f"m{i}" for i, col in enumerate(mets)}
    df_roc = df_local.copy().rename(columns=met_map)
    m_names = list(met_map.values())

    # 2. Rebuild the model_spec string
    model_spec = f"""
        SugarAlcohol =~ {" + ".join(m_names)}
        {mediator} ~ SugarAlcohol + AgeAtVisit + Gender_encoded
        outcome ~ {mediator} + SugarAlcohol + AgeAtVisit + Gender_encoded
    """

    # 3. Fit the final model to get Factor Scores
    m_final = Model(model_spec)
    m_final.fit(df_roc)
    scores = m_final.predict(df_roc)
    
    # Check if 'SugarAlcohol' exists in predicted scores
    if 'SugarAlcohol' in scores.columns:
        df_roc['Metabolic_Score'] = scores['SugarAlcohol']
    else:
        # Fallback if latent variable estimation is unstable
        df_roc['Metabolic_Score'] = df_roc[m_names].mean(axis=1)

    # 4. Define the three models
    models = {
        'Biomarker_Only': [mediator, 'AgeAtVisit', 'Gender_encoded'],
        'Metabolites_Only': ['Metabolic_Score', 'AgeAtVisit', 'Gender_encoded'],
        'Combined': [mediator, 'Metabolic_Score', 'AgeAtVisit', 'Gender_encoded']
    }

    roc_results = {}

    for name, features in models.items():
        X = df_local[features]
        y = df_local['outcome']
        
        lr = LogisticRegression()
        lr.fit(X, y)
        probs = lr.predict_proba(X)[:, 1]
        
        fpr, tpr, _ = roc_curve(y, probs)
        roc_results[name] = {'fpr': fpr.tolist(), 'tpr': tpr.tolist(), 'auc': auc(fpr, tpr)}

    # 5. Save as a pickle for the QMD plot
    fname = f"SEM_mediation_roc_data_{contrast_name}_{mediator}.pkl"
    with open(f"{OUTPUT_DIR}/{fname}", 'wb') as f:
        pickle.dump(roc_results, f)


# ----------------------------------------
# Run loop
# ----------------------------------------
for g1, g2 in CONTRASTS:
    contrast_name = f"{g1}_vs_{g2}"
    
    # Subset data for this contrast
    df_contrast = df[df['fullClass'].isin([g1, g2])].copy()
    df_contrast['outcome'] = (df_contrast['fullClass'] == g2).astype(int)
    
    for mediator, metabolites in STUDY_CONFIG.items():
        run_mediation(df_contrast, mediator, metabolites, contrast_name)

print("\nAll analyses complete. Files saved to results/modelling/")
