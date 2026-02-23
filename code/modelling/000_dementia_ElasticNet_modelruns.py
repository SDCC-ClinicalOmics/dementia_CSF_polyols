################################################################
####### Python script for ElasticNet models construction #######
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
from sklearn.base import clone

OUTPUT_DIR = "../results/modelling"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# load data
df = pd.read_csv('../data/csf_ratios2.csv', sep=';', encoding='latin1')

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

    cv = KFold(n_splits=5, shuffle=True, random_state=42)

    # True score
    cv_scores = cross_val_score(
        pipeline, X, y,
        scoring="r2",
        cv=cv,
        n_jobs=-1
    )

    true_score = cv_scores.mean()

    permutation_scores = []

    for i in range(n_permutations):

        if i % 50 == 0:
            print(f"Permutation {i}/{n_permutations}")

        y_perm = np.random.permutation(y)

        fold_scores = []

        for train_idx, test_idx in cv.split(X):
            model = clone(pipeline)
            model.fit(X[train_idx], y_perm[train_idx])

            pred = model.predict(X[test_idx])
            fold_scores.append(r2_score(y_perm[test_idx], pred))

        permutation_scores.append(np.mean(fold_scores))

    permutation_scores = np.array(permutation_scores)

    p_value = (np.sum(permutation_scores >= true_score) + 1) / (n_permutations + 1)

    return true_score, permutation_scores, p_value

def run_configuration(X, y, subjectid, feature_set_name, outcome_name, feature_columns):

    print("\n=============================================")
    print(f"Running: {feature_set_name} | Outcome: {outcome_name}")
    print("=============================================")

    pipeline = build_pipeline()

     # -------------------------------------------------
    # Fit final model on full data
    # -------------------------------------------------
    pipeline.fit(X, y)

    # -------------------------------------------------
    # Cross-validated predictions (Metabolite Risk Score)
    # -------------------------------------------------
    mrs_scores = cross_val_predict(
        pipeline, X, y,
        cv=10,
        n_jobs=-1
    )

    # Save MRS as independent file
    mrs_df = pd.DataFrame({
        "DDBBno": subjectid,
        "MRS": mrs_scores,
        "Outcome": y
    })

    mrs_fname = os.path.join(
        OUTPUT_DIR,
        f"MRS_{feature_set_name}_{outcome_name}.txt"
    )

    mrs_df.to_csv(mrs_fname, index=False, sep="\t")

    print(f"Saved MRS scores to: {mrs_fname}")

    # Optimal parameters
    alpha = pipeline.named_steps["enet"].alpha_
    l1_ratio = pipeline.named_steps["enet"].l1_ratio_

    print(f"Optimal alpha: {alpha:.6f}")
    print(f"Optimal l1_ratio: {l1_ratio:.6f}")

    # -------------------------------------------------
    # Predictions and performance metrics
    # -------------------------------------------------
    y_pred = pipeline.predict(X)

    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    mae = mean_absolute_error(y, y_pred)

    print(f"R2: {r2:.4f} | RMSE: {rmse:.4f} | MAE: {mae:.4f}")

    # -------------------------------------------------
    # Coefficients
    # -------------------------------------------------
    coefficients = pipeline.named_steps["enet"].coef_
    intercept = pipeline.named_steps["enet"].intercept_

    non_zero = int((np.abs(coefficients) > 1e-6).sum())

    # -------------------------------------------------
    # Scaled X for SHAP (for plotting in qmd)
    # -------------------------------------------------
    X_scaled = pipeline.named_steps["scaler"].transform(X)

    # -------------------------------------------------
    # Permutation testing
    # -------------------------------------------------
    true_score, perm_scores, pval = permutation_test(
        pipeline, X, y,
        n_permutations=1000
    )

    # -------------------------------------------------
    # Save results
    # -------------------------------------------------

    result = {
        "feature_set": feature_set_name,
        "outcome": outcome_name,

        # model info
        "fitted_model": pipeline,
        "alpha": alpha,
        "l1_ratio": l1_ratio,

        # predictions & metrics
        "y_true": y,
        "y_pred": y_pred,
        "r2": r2,
        "rmse": rmse,
        "mae": mae,

        # coefficients
        "coefficients": coefficients,
        "intercept": intercept,
        "non_zero_coefficients": non_zero,
        "feature_columns": feature_columns,

        # data for SHAP
        "X_scaled": X_scaled,

        # permutation test
        "true_r2_cv": true_score,
        "permutation_scores": perm_scores,
        "p_value": pval
    }

    fname = os.path.join(OUTPUT_DIR, f"results_{feature_set_name}_{outcome_name}.pkl")
    pickle.dump(result, open(fname, "wb"))

    print(f"Saved results to: {fname}")

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
    'd.Glucose', 'gluc.erythitol', 'd.Galactose', 'gluc.sorbitol', 
    'D.....Xylose', 'D.....Ribofuranose', 'Sorbitol', 'L.....Arabitol', 
    'D.Threitol', 'meso.Erythritol', 'X3.Deoxy.erythro.pentitol', 
    'X1.5.Anhydrohexitol', 'Ribonic.acid', 'L.Threonic.acid', 'Glyceric.acid'
]
confounders = []

print(f'\n15 Metabolites: {metabolites}')

X1 = df[metabolites].values
subjectid=df["DDBBno"].values

for outcome_name, y in outcomes.items():
    res = run_configuration(X1, y, subjectid, feature_set_name, outcome_name, metabolites)
    summary_rows.append(res)

# =======================================================
# RUN 2 - METABOLITES + AGE + GENDER
# =======================================================
print("\n--- RUN 2: Metabolites + Age + Gender ---")
feature_set_name = "metabolites_plus_age_gender"

confounders = ['AgeAtVisit', 'Gender']
working_cols = metabolites + confounders + ['DDBBno']
df_working = df[working_cols].copy()

## Encode categorical variables (Gender)
df_working['Gender_encoded'] = (df_working['Gender'] == 'Male').astype(int) 

feature_columns = metabolites + ['AgeAtVisit', 'Gender_encoded']

X2 = df_working[feature_columns].values
subjectid=df_working["DDBBno"].values

for outcome_name, y in outcomes.items():
    res = run_configuration(X2, y, subjectid, feature_set_name, outcome_name, feature_columns)
    summary_rows.append(res)


# =======================================================
# RUN 3 - METABOLITES + AGE + GENDER + BMI + STATINS
# =======================================================
print("\n--- RUN 3: Metabolites + Age + Gender + BMI + Statins ---")
feature_set_name = "metabolites_plus_full_covariates"

confounders = ['AgeAtVisit', 'Gender', 'BMI', 'Statins']
working_cols = metabolites + confounders + ['TAU', 'PTAU', 'ABETA42', 'DDBBno']
df_working = df[working_cols].copy()

df_complete = df_working.dropna()

# Encode Gender (Male=1, Female=0)
df_complete['Gender_encoded'] = (df_complete['Gender'] == 'Male').astype(int)

# Statins is already numeric (0.0, 1.0), just ensure it's integer
df_complete['Statins_encoded'] = df_complete['Statins'].astype(int)

feature_columns = metabolites + ['AgeAtVisit', 'Gender_encoded', 'BMI', 'Statins_encoded']

X3 = df_complete[feature_columns].values
subjectid=df_complete["DDBBno"].values

for outcome_name in outcomes.keys():
    y = df_complete[outcome_name].values
    res = run_configuration(X3, y, subjectid, feature_set_name, outcome_name, feature_columns)
    summary_rows.append(res)

# -------------------------------------------------------
# Save overall summary
# -------------------------------------------------------

summary = pd.DataFrame(summary_rows)
summary_file = os.path.join(OUTPUT_DIR, "all_model_summaries.csv")
summary.to_csv(summary_file, index=False)

print("\nAll analyses completed!")
print("Summary saved to: {summary_file}")