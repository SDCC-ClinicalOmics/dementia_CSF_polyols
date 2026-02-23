################################################################
######## Python script for ElasticNet models validation ########
################################################################

# load libraries
import os
import numpy as np
import pandas as pd
import pickle
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold, cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.utils import resample
from sklearn.base import clone

OUTPUT_DIR = "../results/models_validation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# load data
df = pd.read_csv("../results/ADNI_genetics/Cruchaga_metabolomics_relevant.txt", sep = "\t", encoding = 'latin1')

# define modelling functions
def build_pipeline():
    return Pipeline([
        ("scaler", StandardScaler()),
        ("enet", ElasticNetCV(cv=5, 
                              random_state=42,
                              max_iter=1000,
                              n_jobs = -1))
    ])


def permutation_test(pipeline, X, y, n_permutations=1000):
    #
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    #
    # True score
    cv_scores = cross_val_score(
        pipeline, X, y,
        scoring="r2",
        cv=cv,
        n_jobs=-1
    )
    #
    true_score = cv_scores.mean()
    #
    permutation_scores = []
    #
    for i in range(n_permutations):
        #
        if i % 50 == 0:
            print(f"Permutation {i}/{n_permutations}")
        #
        y_perm = np.random.permutation(y)
        #
        fold_scores = []
        #
        for train_idx, test_idx in cv.split(X):
            model = clone(pipeline)
            model.fit(X[train_idx], y_perm[train_idx])
            #
            pred = model.predict(X[test_idx])
            fold_scores.append(r2_score(y_perm[test_idx], pred))
        #
        permutation_scores.append(np.mean(fold_scores))
    #
    permutation_scores = np.array(permutation_scores)
    #
    p_value = (np.sum(permutation_scores >= true_score) + 1) / (n_permutations + 1)
    #
    return true_score, permutation_scores, p_value

def run_configuration(X, y, subjectid, feature_set_name, outcome_name, feature_columns):
    #
    print("\n=============================================")
    print(f"Running: {feature_set_name} | Outcome: {outcome_name}")
    print("=============================================")
    #
    pipeline = build_pipeline()
    #
    # -------------------------------------------------
    # Fit final model on full data
    # -------------------------------------------------
    pipeline.fit(X, y)
    #
    # Optimal parameters
    alpha = pipeline.named_steps["enet"].alpha_
    l1_ratio = pipeline.named_steps["enet"].l1_ratio_
    #
    print(f"Optimal alpha: {alpha:.6f}")
    print(f"Optimal l1_ratio: {l1_ratio:.6f}")
    #
    # -------------------------------------------------
    # Predictions and performance metrics
    # -------------------------------------------------
    y_pred = pipeline.predict(X)
    #
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    mae = mean_absolute_error(y, y_pred)
    #
    print(f"R2: {r2:.4f} | RMSE: {rmse:.4f} | MAE: {mae:.4f}")
    #
    # -------------------------------------------------
    # Coefficients
    # -------------------------------------------------
    coefficients = pipeline.named_steps["enet"].coef_
    intercept = pipeline.named_steps["enet"].intercept_
    #
    non_zero = int((np.abs(coefficients) > 1e-6).sum())
    #
    # -------------------------------------------------
    # Scaled X for SHAP (for plotting in qmd)
    # -------------------------------------------------
    X_scaled = pipeline.named_steps["scaler"].transform(X)
    #
    # -------------------------------------------------
    # Permutation testing
    # -------------------------------------------------
    true_score, perm_scores, pval = permutation_test(
        pipeline, X, y,
        n_permutations=1000
    )
    #
    # -------------------------------------------------
    # Save results
    # -------------------------------------------------
    #
    result = {
        "feature_set": feature_set_name,
        "outcome": outcome_name,
        #
        # model info
        "fitted_model": pipeline,
        "alpha": alpha,
        "l1_ratio": l1_ratio,
        #
        # predictions & metrics
        "y_true": y,
        "y_pred": y_pred,
        "r2": r2,
        "rmse": rmse,
        "mae": mae,
        #
        # coefficients
        "coefficients": coefficients,
        "intercept": intercept,
        "non_zero_coefficients": non_zero,
        "feature_columns": feature_columns,
        #
        # data for SHAP
        "X_scaled": X_scaled,
        #
        # permutation test
        "true_r2_cv": true_score,
        "permutation_scores": perm_scores,
        "p_value": pval
    }
    #
    fname = os.path.join(OUTPUT_DIR, f"results_validation_{feature_set_name}_{outcome_name}.pkl")
    pickle.dump(result, open(fname, "wb"))
    #
    print(f"Saved results to: {fname}")
    #
    # Summary line for CSV
    return {
        "feature_set": feature_set_name,
        "outcome": outcome_name,
        "true_r2": true_score,
        "p_value": pval,
        "alpha": alpha,
        "l1_ratio": l1_ratio,
        "r2_full": r2,
        "rmse": rmse,
        "mae": mae,
        "non_zero_coefs": non_zero
    }

# -------------------------------------------------------
# Filter working DF for NAs
# -------------------------------------------------------
df = df.dropna(subset=['TAU', 'PTAU', 'ABETA42', 'RID', 'PTGENDER', 'AGE', 'myo.inositol', 'glucose', 
                       'glycerate', 'erythritol', 'ribonate', 'threonate', 'mannitol.sorbitol', 
                       'arabitol.xylitol', 'gluc.sorb', 'gluc.ery'])

# -------------------------------------------------------
# Scale targets
# -------------------------------------------------------
df['TAU'] = StandardScaler().fit_transform(np.log(df['TAU']).values.reshape(-1,1)).ravel() 
df['PTAU'] = StandardScaler().fit_transform(np.log(df['PTAU']).values.reshape(-1,1)).ravel() 
df['ABETA42'] = StandardScaler().fit_transform(df['ABETA42'].values.reshape(-1,1)).ravel()

# -------------------------------------------------------
# Prepare outcome dictionary
# -------------------------------------------------------
outcomes = {
    "TAU": df["TAU"].values,
    "PTAU": df["PTAU"].values,
    "ABETA42": df["ABETA42"].values
}

summary_rows = []

# =======================================================
# RUN 1 - ONLY METABOLITES
# =======================================================
print("\n--- RUN 1: Metabolites only ---")
feature_set_name = "metabolites_only"

metabolites = [
    'myo.inositol', 'glucose', 'glycerate', 'erythritol', 'ribonate',
    'threonate', 'mannitol.sorbitol', 'arabitol.xylitol', 'gluc.sorb',
    'gluc.ery'
]
confounders = []

print(f'\n10 Metabolites: {metabolites}')

X1 = df[metabolites].values
subjectid=df["RID"].values

for outcome_name, y in outcomes.items():
    res = run_configuration(X1, y, subjectid, feature_set_name, outcome_name, metabolites)
    summary_rows.append(res)

# =======================================================
# RUN 2 - METABOLITES + AGE + GENDER
# =======================================================
print("\n--- RUN 2: Metabolites + Age + Gender ---")
feature_set_name = "metabolites_plus_age_gender"

confounders = ['AGE', 'PTGENDER']
working_cols = metabolites + confounders + ['RID']
df_working = df[working_cols].copy()

## Encode categorical variables (Gender)
df_working['Gender_encoded'] = df_working['PTGENDER'].map({1: 0, 2: 1}) # we recode to 0/1 to clarify this is a binary variable

feature_columns = metabolites + ['AGE', 'Gender_encoded']

X2 = df_working[feature_columns].values
subjectid=df_working["RID"].values

for outcome_name, y in outcomes.items():
    res = run_configuration(X2, y, subjectid, feature_set_name, outcome_name, feature_columns)
    summary_rows.append(res)

# -------------------------------------------------------
# Save overall summary
# -------------------------------------------------------

summary = pd.DataFrame(summary_rows)
summary_file = os.path.join(OUTPUT_DIR, "all_model_summaries.csv")
summary.to_csv(summary_file, index=False)

print("\nAll analyses completed!")

# -------------------------------------------------------
# Bootstrap piece
# -------------------------------------------------------
print("\nNow running the bootstrap replications!")

# set random seed for reproducibility
np.random.seed(357)

# define functions and parameters
alphas = np.logspace(-3, 1, 50)


#### bootstrapping function
def run_bootstrap(X, y, metabolite_names, B):
    #
    n, p = X.shape
    metabolite_idx = [X.columns.get_loc(m) for m in metabolite_names]
    #
    selection_counts = np.zeros(len(metabolite_names))
    coef_matrix = np.zeros((B, len(metabolite_names)))
    r2_sample = [] # to be used for CV
    pred_matrix = np.zeros((B, len(y)))
    #
    for b in range(B):
        #
        # resample subjects
        Xb, yb = resample(X, y, replace=True)
        #
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('enet', ElasticNetCV(
                l1_ratio=[.1, .5, .9, .95, .99, 1],
                alphas=alphas,
                cv=5,
                n_jobs=-1))
        ])
        #
        pipe.fit(Xb, yb)
        #
        model = pipe.named_steps['enet']
        coefs = model.coef_
        #
        # extract metabolite coefficients only
        meta_coefs = coefs[metabolite_idx]
        #
        coef_matrix[b, :] = meta_coefs
        selection_counts += (meta_coefs != 0)
        #
        # compute in-sample R2 for this bootstrap
        preds = pipe.predict(Xb)
        r2_sample.append(r2_score(yb, preds))
        #
        # predicted scores
        pred_matrix[b,:] = pipe.predict(X)
        #
    # selection frequency
    selection_freq = selection_counts / B
    #
    # bootstrap prediction uncertainty
    prediction_sd = pred_matrix.std(axis=0)
    #
    # bootstrap R2 scores (R2-CV)
    y_pred_cv = pred_matrix.mean(axis=0)
    r2_cv = r2_score(y, y_pred_cv)
    #
    return selection_freq, coef_matrix, r2_sample, r2_cv, pred_matrix, prediction_sd

#### saving function
def save_bootstrap_results(filename, outcome_name, results, metabolites, feature_columns):
    #
    obj = {
        "outcome": outcome_name,
        "metabolites": metabolites,
        "feature_columns": feature_columns,
        #
        "selection_freq": results[0],
        "coef_matrix": results[1],
        "r2_sample": results[2],
        "r2_cv": results[3],
        "pred_matrix": results[4],
        "prediction_sd": results[5]
    }
    #
    out_path = os.path.join(OUTPUT_DIR, filename)
    #
    with open(out_path, "wb") as f:
        pickle.dump(obj, f)
    #
    print(f"Saved bootstrap results to: {out_path}")


# -------------------------------------------------------
# Create a working dataframe with only the variables we need
# -------------------------------------------------------
working_cols = metabolites + confounders + ['TAU', 'PTAU', 'ABETA42']
df_working = df[working_cols].copy()

# Encode categorical variables (Gender)
df_working['Gender_encoded'] = df_working['PTGENDER'].map({1: 0, 2: 1}) # we recode to 0/1 to clarify this is a binary variable

# -------------------------------------------------------
# Prepare feature matrix X and target y
# -------------------------------------------------------
feature_columns = metabolites + ['AGE', 'Gender_encoded']

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
    filename="validation_bootstrap_TAU.pkl",
    outcome_name="TAU",
    results=res_ttau,
    metabolites=metabolites,
    feature_columns=feature_columns
)

# PTAU
res_ptau = run_bootstrap(X, y_ptau, metabolites, B=B)
save_bootstrap_results(
    filename="validation_bootstrap_PTAU.pkl",
    outcome_name="PTAU",
    results=res_ptau,
    metabolites=metabolites,
    feature_columns=feature_columns
)

# ABETA42
res_ab42 = run_bootstrap(X, y_ab42, metabolites, B=B)
save_bootstrap_results(
    filename="validation_bootstrap_ABETA42.pkl",
    outcome_name="ABETA42",
    results=res_ab42,
    metabolites=metabolites,
    feature_columns=feature_columns
)

print("\nAll bootstrap analyses completed and saved.")