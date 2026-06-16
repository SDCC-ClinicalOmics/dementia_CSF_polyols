# =============================================================================
# Variance Partitioning: Disentangling Tau and Amyloid Pathways
# Reviewer response analysis — Clos-Garcia et al.
# =============================================================================
# Addresses reviewer concern: "tau associations are often (partially) explained
# by Abeta abnormalities. Show how much variance in tau is explained by Abeta
# and vice versa."
#
# Three complementary analyses:
#   1. Variance partitioning (vegan::varpart) — Venn decomposition
#   2. Hierarchical regression — ΔR² when adding polyols beyond Aβ42
#   3. Partial correlations — individual metabolite–tau associations
#      controlling for Aβ42
# =============================================================================

library(tidyverse)
library(vegan)       # varpart()
library(ppcor)       # pcor.test() for partial correlations
library(eulerr)      # Euler/Venn diagram
library(ggplot2)
library(patchwork)
library(broom)
library(here)

# =============================================================================
# 0. CONFIGURATION — ADAPT THESE TO YOUR DATA
# =============================================================================

# --- Variable names in your data frame ----------------------------------------
OUTCOME_VARS  <- c("tTau", "pTau")         # tau biomarkers (response)
AMYLOID_VAR   <- "Aβ42"                      # Aβ42 column name

COVARIATES    <- c("AgeAtVisit", "Gender", "BMI")      # standard covariate set

# Polyol metabolites selected by Elastic Net (from manuscript Table / Fig 5B)
# Edit to match exact column names in your data frame
POLYOL_VARS <- c(
  "L.....Arabitol",          # β = 0.374 — strongest contributor
  "Ribonic.acid",        # β = 0.123
  "D.Threitol",          # β = 0.122
  "Sorbitol",
  "meso.Erythritol",
  "gluc.erythitol" # ratio
)

# --- Data loading -------------------------------------------------------------
df <- read.table(here("data/csf_ratios2.csv"), header = T, sep = ";")
df$Gender <- ifelse(df$Gender == "Female", 1, 2)
df <- df |> 
  dplyr::mutate(fullClass = case_when(
    fullClass == "HC"  ~ "NCI",
    fullClass == "MCI" ~ "MCI-AD",
    fullClass == "AD"  ~ "Dementia-AD",
    fullClass == "VaD" ~ "MIX"
  )) |> 
  dplyr::rename("tTau" = TAU,
                "pTau" = PTAU,
                "Aβ42" = ABETA42)

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# Z-standardise all continuous variables (consistent with SEM models)
df_z <- df |>
  mutate(across(all_of(c(OUTCOME_VARS, AMYLOID_VAR, POLYOL_VARS,
                          "AgeAtVisit", "BMI")),
                ~ scale(.)[, 1]))

# Confirm no missing values in variables of interest
vars_needed <- c(OUTCOME_VARS, AMYLOID_VAR, POLYOL_VARS, COVARIATES)
df_complete <- df_z |> drop_na(all_of(vars_needed))
cat(sprintf("N complete cases: %d / %d\n", nrow(df_complete), nrow(df_z)))

# Build matrix objects for vegan
X_polyols    <- df_complete |> dplyr::select(dplyr::all_of(POLYOL_VARS)) |> as.matrix()
X_amyloid    <- df_complete |> dplyr::select(dplyr::all_of(AMYLOID_VAR)) |> as.matrix()
X_covariates <- df_complete |> dplyr::select(dplyr::all_of(COVARIATES))  |> as.matrix()

# =============================================================================
# 2. VARIANCE PARTITIONING  (vegan::varpart)
# =============================================================================
# Three-fraction partitioning: [Polyols | Aβ42 | Covariates]
# Fractions:
#   [a]  = unique to Polyols
#   [b]  = unique to Aβ42
#   [c]  = unique to Covariates
#   [ab] = shared Polyols ∩ Aβ42
#   etc.

run_varpart <- function(outcome, df_complete, X_polyols, X_amyloid, X_covariates) {
  Y <- df_complete[[outcome]]
  vp <- varpart(Y, X_polyols, X_amyloid, X_covariates)
  vp
}

vp_results <- map(OUTCOME_VARS, ~ run_varpart(.x, df_complete,
                                               X_polyols, X_amyloid,
                                               X_covariates))
names(vp_results) <- OUTCOME_VARS

# Print summaries
walk2(vp_results, OUTCOME_VARS, function(vp, name) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Variance Partitioning for", name, "\n")
  print(vp)
})

# Extract key fractions into a tidy table
extract_fractions <- function(vp, outcome) {
  tbl <- as.data.frame(vp$part$indfract)
  
  # Locate the Adj.R.squared column by name (case-insensitive) with
  # positional fallback (it is always the 3rd column in vegan's output).
  adj_col <- grep("Adj\\.R", colnames(tbl), value = TRUE, ignore.case = TRUE)
  r2_vals <- if (length(adj_col) > 0) {
    as.numeric(tbl[[adj_col[1]]])
  } else {
    as.numeric(tbl[[3]])   # positional fallback: Df | R.squared | Adj.R.squared | Testable
  }
  
  tibble(
    outcome     = outcome,
    fraction    = rownames(tbl),
    adj_R2      = r2_vals,
    description = c(
      "Unique to Polyols",
      "Unique to Aβ42",
      "Unique to Covariates",
      "Shared: Polyols ∩ Aβ42",
      "Shared: Polyols ∩ Covariates",
      "Shared: Aβ42 ∩ Covariates",
      "Shared: all three",
      "Residuals (unexplained)"
    )
  )
}

varpart_table <- map2_dfr(vp_results, OUTCOME_VARS, extract_fractions)
print(varpart_table)
write_csv(varpart_table, "results/varpart_fractions.csv")

# =============================================================================
# 3. HIERARCHICAL REGRESSION  — ΔR² approach
# =============================================================================
# Four nested models per tau outcome:
#   M0: tau ~ Covariates                             (baseline)
#   M1: tau ~ Covariates + Aβ42                      (amyloid step)
#   M2: tau ~ Covariates + Polyols                   (polyol step)
#   M3: tau ~ Covariates + Aβ42 + Polyols            (full model)
#
# Key comparisons:
#   M3 vs M1 → ΔR² of polyols over-and-above Aβ42
#   M3 vs M2 → ΔR² of Aβ42 over-and-above polyols

polyol_formula  <- paste(POLYOL_VARS,   collapse = " + ")
covariate_formula <- paste(COVARIATES,  collapse = " + ")

build_hierarchical_models <- function(outcome, df) {
  f_M0 <- as.formula(paste(outcome, "~", covariate_formula))
  f_M1 <- as.formula(paste(outcome, "~", covariate_formula, "+", AMYLOID_VAR))
  f_M2 <- as.formula(paste(outcome, "~", covariate_formula, "+", polyol_formula))
  f_M3 <- as.formula(paste(outcome, "~", covariate_formula, "+", AMYLOID_VAR,
                            "+", polyol_formula))

  M0 <- lm(f_M0, data = df)
  M1 <- lm(f_M1, data = df)
  M2 <- lm(f_M2, data = df)
  M3 <- lm(f_M3, data = df)

  r2 <- function(m) summary(m)$r.squared
  ar2 <- function(m) summary(m)$adj.r.squared

  # F-test for each increment
  f_test_M1_vs_M0 <- anova(M0, M1)
  f_test_M2_vs_M0 <- anova(M0, M2)
  f_test_M3_vs_M1 <- anova(M1, M3)   # polyols over-and-above Aβ42
  f_test_M3_vs_M2 <- anova(M2, M3)   # Aβ42 over-and-above polyols

  tibble(
    outcome       = outcome,
    model         = c("M0: Covariates only",
                      "M1: + Aβ42",
                      "M2: + Polyols",
                      "M3: + Aβ42 + Polyols (full)"),
    R2            = c(r2(M0), r2(M1), r2(M2), r2(M3)),
    adj_R2        = c(ar2(M0), ar2(M1), ar2(M2), ar2(M3)),
    delta_R2_vs_M0 = c(NA,
                        r2(M1) - r2(M0),
                        r2(M2) - r2(M0),
                        r2(M3) - r2(M0)),
    # Key increments addressing the reviewer
    F_stat        = c(NA,
                      f_test_M1_vs_M0$F[2],
                      f_test_M2_vs_M0$F[2],
                      NA),
    p_F           = c(NA,
                      f_test_M1_vs_M0$`Pr(>F)`[2],
                      f_test_M2_vs_M0$`Pr(>F)`[2],
                      NA)
  )
}

# Build dedicated comparison table for key increments
build_increment_table <- function(outcome, df) {
  f_M0 <- as.formula(paste(outcome, "~", covariate_formula))
  f_M1 <- as.formula(paste(outcome, "~", covariate_formula, "+", AMYLOID_VAR))
  f_M2 <- as.formula(paste(outcome, "~", covariate_formula, "+", polyol_formula))
  f_M3 <- as.formula(paste(outcome, "~", covariate_formula, "+", AMYLOID_VAR,
                            "+", polyol_formula))

  M0 <- lm(f_M0, data = df)
  M1 <- lm(f_M1, data = df)
  M2 <- lm(f_M2, data = df)
  M3 <- lm(f_M3, data = df)

  r2 <- function(m) summary(m)$r.squared

  anova_poly_over_ab  <- anova(M1, M3)  # polyols over-and-above Aβ42
  anova_ab_over_poly  <- anova(M2, M3)  # Aβ42 over-and-above polyols

  tibble(
    outcome              = outcome,
    comparison           = c(
      "Polyols unique (over-and-above Aβ42)",
      "Aβ42 unique (over-and-above Polyols)"
    ),
    delta_R2             = c(
      r2(M3) - r2(M1),
      r2(M3) - r2(M2)
    ),
    F_stat               = c(
      anova_poly_over_ab$F[2],
      anova_ab_over_poly$F[2]
    ),
    df_num               = c(
      anova_poly_over_ab$Df[2],
      anova_ab_over_poly$Df[2]
    ),
    df_den               = c(
      anova_poly_over_ab$Res.Df[2],
      anova_ab_over_poly$Res.Df[2]
    ),
    p_value              = c(
      anova_poly_over_ab$`Pr(>F)`[2],
      anova_ab_over_poly$`Pr(>F)`[2]
    )
  )
}

hier_models   <- map_dfr(OUTCOME_VARS, build_hierarchical_models, df = df_complete)
increment_tbl <- map_dfr(OUTCOME_VARS, build_increment_table,     df = df_complete)

cat("\n\n--- Hierarchical Regression Summary ---\n")
print(hier_models |> mutate(across(where(is.numeric), ~ round(.x, 3))))

cat("\n--- Key Increments (Reviewer Response) ---\n")
print(increment_tbl |> mutate(across(where(is.numeric), ~ round(.x, 4))))

write_csv(hier_models,   "results/hierarchical_regression.csv")
write_csv(increment_tbl, "results/increment_table.csv")

# =============================================================================
# 4. PARTIAL CORRELATIONS  — individual polyol × tau, controlling for Aβ42
# =============================================================================

run_partial_cor <- function(outcome, metabolite, df, control_var) {
  # partial correlation: metabolite ~ outcome | control_var + covariates
  vars_to_control <- c(control_var, COVARIATES)
  dat <- df |>
    dplyr::select(dplyr::all_of(c(outcome, metabolite, vars_to_control))) |>
    drop_na()

  res <- pcor.test(
    x = dat[[outcome]],
    y = dat[[metabolite]],
    z = dat[, vars_to_control]
  )

  tibble(
    outcome        = outcome,
    metabolite     = metabolite,
    partial_r      = res$estimate,
    p_value        = res$p.value,
    statistic      = res$statistic,
    n              = res$n,
    control_for    = paste(vars_to_control, collapse = " + ")
  )
}

# Grid of all combinations
grid <- expand_grid(outcome = OUTCOME_VARS, metabolite = POLYOL_VARS)

# Version 1: control for covariates only (unadjusted for Aβ42)
partial_cor_unadj <- pmap_dfr(grid, function(outcome, metabolite) {
  run_partial_cor(outcome, metabolite, df_complete, character(0)) |>
    mutate(adjustment = "Covariates only")
})

# Version 2: control for Aβ42 + covariates (the key reviewer-requested analysis)
partial_cor_ab42 <- pmap_dfr(grid, function(outcome, metabolite) {
  run_partial_cor(outcome, metabolite, df_complete, AMYLOID_VAR) |>
    mutate(adjustment = "Aβ42 + Covariates")
})

partial_cor_all <- bind_rows(partial_cor_unadj, partial_cor_ab42) |>
  group_by(adjustment) |>
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significant = p_adj < 0.05) |>
  ungroup()

cat("\n--- Partial Correlations (polyols × tau, controlling for Aβ42) ---\n")
print(partial_cor_all |>
  dplyr::select(outcome, metabolite, adjustment, partial_r, p_value, p_adj, significant) |>
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3))))

write_csv(partial_cor_all, "results/partial_correlations.csv")

# =============================================================================
# 5. VISUALISATION
# =============================================================================

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 5a. Euler diagram of variance partitioning ---
plot_euler <- function(vp, outcome) {
  fracs <- vp$part$indfract$Adj.R.square
  # fracs order from vegan: [a],[b],[c],[ab],[ac],bc],[abc]
  unique_poly <- max(fracs[1], 0)
  unique_ab42 <- max(fracs[2], 0)
  shared      <- max(fracs[4], 0)   # Polyols ∩ Aβ42 (marginalising covariates)

  euler_fit <- euler(c(
    "Polyols"   = unique_poly,
    "Aβ42"      = unique_ab42,
    "Polyols&Aβ42" = shared
  ))

  plot(euler_fit,
       fills  = list(fill = c("#4E79A7", "#F28E2B", "#76B7B2"), alpha = 0.6),
       labels = list(
         labels   = c(
           sprintf("Polyols\n%.1f%%", unique_poly * 100),
           sprintf("Aβ42\n%.1f%%",   unique_ab42 * 100),
           sprintf("Shared\n%.1f%%", shared * 100)
         ),
         fontsize = 12
       ),
       main = sprintf("Variance in %s (Adj. R²)", outcome))
}

# Save Euler plots
pdf("results/figures/euler_varpart_tTau.pdf", width = 5, height = 5)
plot_euler(vp_results[["tTau"]], "tTau")
dev.off()

pdf("results/figures/euler_varpart_pTau.pdf", width = 5, height = 5)
plot_euler(vp_results[["pTau"]], "pTau")
dev.off()

# --- 5b. Stacked bar chart of hierarchical R² ---
plot_data_hier <- hier_models |>
  dplyr::select(outcome, model, R2) |>
  pivot_wider(names_from = model, values_from = R2) |>
  pivot_longer(-outcome, names_to = "model", values_to = "R2") |>
  dplyr::mutate(model = factor(model, levels = c(
    "M0: Covariates only",
    "M1: + Aβ42",
    "M2: + Polyols",
    "M3: + Aβ42 + Polyols (full)"
  )))

p_hier <- ggplot(hier_models,
                 aes(x = model, y = R2, fill = outcome, group = outcome)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%.2f", R2)),
            position = position_dodge(0.7), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(tTau = "#4E79A7", pTau = "#F28E2B")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title    = "Hierarchical R² — Polyols vs. Aβ42",
    subtitle = "Unique contributions to tau variance",
    x        = NULL, y = "R² (proportion of variance explained)",
    fill     = "Tau biomarker"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave("results/figures/hierarchical_R2.pdf", p_hier, width = 8, height = 5)

# --- 5c. Partial correlation forest plot ---
p_pcor <- partial_cor_all |>
  dplyr::filter(outcome %in% OUTCOME_VARS) |>
  dplyr::mutate(
    metabolite  = factor(metabolite, levels = rev(POLYOL_VARS)),
    adjustment  = factor(adjustment,
                         levels = c("Covariates only", "Aβ42 + Covariates")),
    sig_label   = ifelse(significant, "*", "")
  ) |>
  ggplot(aes(x = partial_r, y = metabolite,
             colour = outcome, shape = adjustment, alpha = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sig_label),
            position = position_dodge(width = 0.5),
            vjust = -0.8, size = 5, show.legend = FALSE) +
  scale_colour_manual(values = c(tTau = "#4E79A7", pTau = "#F28E2B")) +
  scale_shape_manual(values = c("Covariates only" = 16, "Aβ42 + Covariates" = 17)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
  facet_wrap(~ outcome, ncol = 2) +
  labs(
    title    = "Partial correlations: polyols × tau",
    subtitle = "Before and after controlling for Aβ42  (* BH-adjusted p < 0.05)",
    x        = "Partial correlation (r)",
    y        = NULL,
    colour   = "Tau biomarker",
    shape    = "Adjustment"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#f0f0f0"))

ggsave("results/figures/partial_cor_forest.pdf", p_pcor, width = 9, height = 5)

# --- 5d. Summary table for supplement -----------------------------------------
summary_table <- increment_tbl |>
  dplyr::mutate(
    across(c(delta_R2, F_stat), ~ round(.x, 3)),
    p_label = ifelse(p_value < 0.001, "< 0.001",
              ifelse(p_value < 0.01,  "< 0.01",
                     sprintf("= %.3f", p_value)))
  ) |>
  dplyr::select(outcome, comparison, delta_R2, F_stat, df_num, df_den, p_label)

cat("\n--- Table for Reviewer Response / Supplement ---\n")
print(summary_table)
write_csv(summary_table, "results/increment_summary_for_paper.csv")

cat("\n\nAll outputs saved to results/ and results/figures/\n")
