################################################################
###### Python script for ElasticNet models bootstrapping #######
################################################################

# import libraries
import os as os
import numpy as np
import pandas as pd

from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.utils import resample
from sklearn.metrics import r2_score, mean_squared_error

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

import pickle

OUTPUT_DIR = "../results/modelling"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# load data
df = pd.read_csv('../data/csf_ratios2.csv', sep=';', encoding='latin1')

# set random seed for reproducibility
np.random.seed(235)

# define functions and parameters
alphas = np.logspace(-3, 1, 50)

#### bootstrapping function
def run_bootstrap(X, y, metabolite_names, B):

    n, p = X.shape
    metabolite_idx = [X.columns.get_loc(m) for m in metabolite_names]
    
    selection_counts = np.zeros(len(metabolite_names))
    coef_matrix = np.zeros((B, len(metabolite_names)))
    r2_sample = [] # to be used for CV
    pred_matrix = np.zeros((B, len(y)))

    for b in range(B):

        # resample subjects
        Xb, yb = resample(X, y, replace=True)

        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('enet', ElasticNetCV(
                l1_ratio=[.1, .5, .9, .95, .99, 1],
                alphas=alphas,
                cv=5,
                n_jobs=-1))
        ])

        pipe.fit(Xb, yb)

        model = pipe.named_steps['enet']
        coefs = model.coef_

        # extract metabolite coefficients only
        meta_coefs = coefs[metabolite_idx]

        coef_matrix[b, :] = meta_coefs
        selection_counts += (meta_coefs != 0)

        # compute in-sample R2 for this bootstrap
        preds = pipe.predict(Xb)
        r2_sample.append(r2_score(yb, preds))

        # predicted scores
        pred_matrix[b,:] = pipe.predict(X)

    # selection frequency
    selection_freq = selection_counts / B

    # bootstrap prediction uncertainty
    prediction_sd = pred_matrix.std(axis=0)

    # bootstrap R2 scores (R2-CV)
    y_pred_cv = pred_matrix.mean(axis=0)
    r2_cv = r2_score(y, y_pred_cv)

    return selection_freq, coef_matrix, r2_sample, r2_cv, pred_matrix, prediction_sd

#### saving function
def save_bootstrap_results(filename, outcome_name, results, metabolites, feature_columns):

    obj = {
        "outcome": outcome_name,
        "metabolites": metabolites,
        "feature_columns": feature_columns,

        "selection_freq": results[0],
        "coef_matrix": results[1],
        "r2_sample": results[2],
        "r2_cv": results[3],
        "pred_matrix": results[4],
        "prediction_sd": results[5]
    }

    out_path = os.path.join(OUTPUT_DIR, filename)

    with open(out_path, "wb") as f:
        pickle.dump(obj, f)

    print(f"Saved bootstrap results to: {out_path}")

# -------------------------------------------------------
# Scale targets
# -------------------------------------------------------
df['TAU'] = StandardScaler().fit_transform(np.log(df['TAU']).values.reshape(-1,1)).ravel() 
df['PTAU'] = StandardScaler().fit_transform(np.log(df['PTAU']).values.reshape(-1,1)).ravel() 
df['ABETA42'] = StandardScaler().fit_transform(df['ABETA42'].values.reshape(-1,1)).ravel()

# -------------------------------------------------------
# Define the 15 metabolites and permitted confounders
# -------------------------------------------------------
metabolites = [
    'd.Glucose', 'gluc.erythitol', 'd.Galactose', 'gluc.sorbitol', 
    'D.....Xylose', 'D.....Ribofuranose', 'Sorbitol', 'L.....Arabitol', 
    'D.Threitol', 'meso.Erythritol', 'X3.Deoxy.erythro.pentitol', 
    'X1.5.Anhydrohexitol', 'Ribonic.acid', 'L.Threonic.acid', 'Glyceric.acid'
]

confounders = ["AgeAtVisit", "Gender"]

print(f"\n15 Metabolites: {metabolites}")
print(f"\n2 Confounders: {confounders}")

# -------------------------------------------------------
# Create a working dataframe with only the variables we need
# -------------------------------------------------------
working_cols = metabolites + confounders + ['TAU', 'PTAU', 'ABETA42']
df_working = df[working_cols].copy()

# Encode categorical variables (Gender)
df_working['Gender_encoded'] = (df_working['Gender'] == 'Male').astype(int) 

# -------------------------------------------------------
# Prepare feature matrix X and target y
# -------------------------------------------------------
feature_columns = metabolites + ['AgeAtVisit', 'Gender_encoded']

X = df_working[feature_columns]
y_ttau = df_working['TAU'].values
y_ptau = df_working['PTAU'].values
y_ab42 = df_working['ABETA42'].values

# -------------------------------------------------------
# Run bootstrapping
# -------------------------------------------------------
B = 1000

print("\nRunning bootstraps...\n")
# TAU
res_ttau = run_bootstrap(X, y_ttau, metabolites, B=B)
save_bootstrap_results(
    filename="bootstrap_TAU.pkl",
    outcome_name="TAU",
    results=res_ttau,
    metabolites=metabolites,
    feature_columns=feature_columns
)

# PTAU
res_ptau = run_bootstrap(X, y_ptau, metabolites, B=B)
save_bootstrap_results(
    filename="bootstrap_PTAU.pkl",
    outcome_name="PTAU",
    results=res_ptau,
    metabolites=metabolites,
    feature_columns=feature_columns
)

# ABETA42
res_ab42 = run_bootstrap(X, y_ab42, metabolites, B=B)
save_bootstrap_results(
    filename="bootstrap_ABETA42.pkl",
    outcome_name="ABETA42",
    results=res_ab42,
    metabolites=metabolites,
    feature_columns=feature_columns
)

print("\nAll bootstrap analyses completed and saved.")