# =============================================================================
# Stage-Stratified Polyol–Tau Analysis
# Reviewer 3, Comment 3 — Clos-Garcia et al.
# =============================================================================
# Addresses: "Separate analysis according to diagnostic stage — controls,
# MCI, AD dementia, mixed dementia — to ascertain if the findings are
# limited to any disease stage."
#
# Three outputs:
#   1. Within-group partial correlations (polyols × tau/pTau/Aβ42 per stage)
#   2. Stage-trend test: does polyol–tau effect size increase with severity?
#   3. Publication-ready figures: heatmap, dot-CI plot, trend line
# =============================================================================

library(here)
library(tidyverse)
library(ppcor)
library(broom)
library(ggplot2)
library(patchwork)

# =============================================================================
# 0. CONFIGURATION
# Match column names to the varpart script and the Python SEM script
# =============================================================================

# --- Diagnostic column -------------------------------------------------------
DIAGNOSIS_COL <- "fullClass"     # as in the Python SEM script
STAGE_ORDER   <- c("HC", "MCI", "AD", "VaD")   # ordered by severity
# If your R data uses "NCI" instead of "HC", change the first value:
# STAGE_ORDER <- c("NCI", "MCI", "AD", "VaD")

# Display labels for figures (matched to STAGE_ORDER)
STAGE_LABELS  <- c("HC", "MCI", "AD", "VaD")

# Numeric coding for trend test (0 = healthiest, 3 = most severe)
STAGE_NUMERIC <- setNames(0:3, STAGE_ORDER)

# --- Variable names (same as varpart script) ---------------------------------
POLYOL_VARS <- c(
  "L.....Arabitol",          # β = 0.374 — strongest contributor
  "Ribonic.acid",        # β = 0.123
  "D.Threitol",          # β = 0.122
  "Sorbitol",
  "meso.Erythritol",
  "gluc.erythitol" # ratio
)
TAU_VARS     <- c("TAU", "PTAU")
AMYLOID_VAR  <- "ABETA42"
ALL_OUTCOMES <- c(TAU_VARS, AMYLOID_VAR)
COVARIATES   <- c("AgeAtVisit", "Gender", "BMI")

# Clean metabolite labels for figures
MET_LABELS <- c(
  L.....Arabitol           = "L-Arabitol",
  Ribonic.acid         = "Ribonic acid",
  D.Threitol           = "D-Threitol",
  Sorbitol             = "Sorbitol",
  meso.Erythritol      = "meso-Erythritol",
  gluc.erythitol = "Gluc./Erythritol"
)


# --- Data loading -------------------------------------------------------------
df <- read.table(here("data/csf_ratios2.csv"), header = T, sep = ";")
df$Gender <- ifelse(df$Gender == "Female", 1, 2)

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# Z-standardise all continuous variables (consistent with SEM models)
df_z <- df |>
  mutate(across(all_of(c(TAU_VARS, AMYLOID_VAR, POLYOL_VARS,
                         "AgeAtVisit", "BMI")),
                ~ scale(.)[, 1]))

# Confirm no missing values in variables of interest
vars_needed <- c(TAU_VARS, AMYLOID_VAR, POLYOL_VARS, COVARIATES)
df_complete <- df_z |> drop_na(all_of(vars_needed))


# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

df_staged <- df_complete %>%
  dplyr::filter(.data[[DIAGNOSIS_COL]] %in% STAGE_ORDER) %>%
  dplyr::mutate(
    stage     = factor(.data[[DIAGNOSIS_COL]], levels = STAGE_ORDER),
    stage_num = STAGE_NUMERIC[as.character(.data[[DIAGNOSIS_COL]])]
  ) %>%
  # Z-score within the full dataset (consistent with main analysis)
  dplyr::mutate(across(all_of(c(POLYOL_VARS, ALL_OUTCOMES, "AgeAtVisit", "BMI")),
                ~ scale(.)[, 1]))

# Group sizes — critical for power caveats
group_ns <- df_staged %>%
  count(stage) %>%
  mutate(
    label     = paste0(as.character(stage), "\n(n=", n, ")"),
    power_note = case_when(
      n < 25 ~ "⚠ low power",
      n < 40 ~ "moderate power",
      TRUE   ~ "adequate"
    )
  )

cat("Group sizes:\n")
print(group_ns)
cat("\nNote: HC and MCI groups are underpowered for within-group inference.\n")
cat("Report directional consistency rather than individual p-values.\n\n")

# Named vector of n-labelled stage names for figure axes
stage_axis_labels <- setNames(group_ns$label, as.character(group_ns$stage))

# =============================================================================
# 2. WITHIN-GROUP PARTIAL CORRELATIONS
# =============================================================================
# For each group × metabolite × biomarker:
#   partial r (polyol ~ biomarker | Age + Sex + BMI)
# 95% CI via Fisher z-transformation (df = n - 2 - p, where p = # covariates)

pcor_with_ci <- function(x, y, z, n_cov) {
  # Returns partial r, Fisher-z 95% CI, and p-value
  res  <- pcor.test(x, y, z)
  r    <- res$estimate
  pval <- res$p.value
  df   <- length(x) - 2 - n_cov
  if (df < 3) return(tibble(r = r, ci_lo = NA, ci_hi = NA, p_value = pval))
  se   <- 1 / sqrt(df - 1)
  z_r  <- atanh(r)
  ci   <- tanh(c(z_r - 1.96 * se, z_r + 1.96 * se))
  tibble(r = r, ci_lo = ci[1], ci_hi = ci[2], p_value = pval)
}

run_staged_pcor <- function(stage_name, metabolite, biomarker, df) {
  dat <- df %>%
    dplyr::filter(stage == stage_name) %>%
    dplyr::select(all_of(c(metabolite, biomarker, COVARIATES))) %>%
    drop_na()

  n <- nrow(dat)
  if (n < 8) return(NULL)   # minimum viable df

  result <- tryCatch(
    pcor_with_ci(dat[[metabolite]], dat[[biomarker]],
                 dat[, COVARIATES], length(COVARIATES)),
    error = function(e) NULL
  )
  if (is.null(result)) return(NULL)

  result %>% mutate(stage = stage_name, metabolite = metabolite,
                    biomarker = biomarker, n = n)
}

# All combinations
grid <- expand_grid(
  stage_name = STAGE_ORDER,
  metabolite = POLYOL_VARS,
  biomarker  = ALL_OUTCOMES
)

pcor_staged <- pmap_dfr(grid, ~ run_staged_pcor(..1, ..2, ..3, df_staged)) %>%
  mutate(stage = factor(stage, levels = STAGE_ORDER)) %>%
  group_by(stage, biomarker) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = p_adj < 0.05
  ) %>%
  ungroup()

# Print tau results
cat("=== Within-group partial correlations (polyols × tau) ===\n")
pcor_staged %>%
  filter(biomarker %in% TAU_VARS) %>%
  mutate(across(c(r, ci_lo, ci_hi, p_adj), ~ round(.x, 3)),
         sig = ifelse(sig, "*", "")) %>%
  arrange(biomarker, metabolite, stage) %>%
  print(n = Inf)

write_csv(pcor_staged, "results/staged_partial_correlations.csv")

# =============================================================================
# 3. STAGE-TREND TEST
# =============================================================================
# Tests whether polyol–tau association STRENGTH changes with disease stage.
# Model: tau ~ polyol * stage_num + covariates
#
# Positive, significant polyol:stage_num interaction = effect strengthens
# as disease progresses. Run separately for the AD spectrum (HC→MCI→AD)
# and the full cohort (including VaD, which has a different aetiology).

run_trend_test <- function(metabolite, biomarker, df, label = "full") {
  dat <- df %>%
    dplyr::select(all_of(c(metabolite, biomarker, "stage_num", COVARIATES))) %>%
    drop_na()

  f <- as.formula(paste(
    biomarker, "~", metabolite, "* stage_num +",
    paste(COVARIATES, collapse = " + ")
  ))
  m <- lm(f, data = dat)

  tidy(m, conf.int = TRUE) %>%
    dplyr::filter(term == paste0(metabolite, ":stage_num")) %>%
    dplyr::mutate(metabolite = metabolite, biomarker = biomarker, cohort = label) %>%
    dplyr::select(cohort, biomarker, metabolite, estimate,
           conf.low, conf.high, p.value)
}

trend_grid <- expand_grid(metabolite = POLYOL_VARS, biomarker = ALL_OUTCOMES)

# Full cohort (HC → MCI → AD → VaD)
trend_full <- pmap_dfr(trend_grid, ~ run_trend_test(..1, ..2, df_staged, "All stages"))

# AD spectrum only (HC → MCI → AD): removes VaD which follows different path
trend_ad_spectrum <- pmap_dfr(trend_grid, ~ run_trend_test(
  ..1, ..2,
  df_staged %>% filter(stage %in% c("HC", "MCI", "AD")),
  "AD spectrum (HC→MCI→AD)"
))

trend_results <- bind_rows(trend_full, trend_ad_spectrum) %>%
  group_by(cohort, biomarker) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         sig   = p_adj < 0.05) %>%
  ungroup()

cat("\n=== Stage-trend test (polyol × disease severity interaction) ===\n")
cat("Positive estimate = polyol–tau association STRENGTHENS with disease stage\n\n")
trend_results %>%
  filter(biomarker %in% TAU_VARS) %>%
  mutate(across(c(estimate, conf.low, conf.high, p_adj), ~ round(.x, 3)),
         sig = ifelse(sig, "*", "")) %>%
  arrange(biomarker, cohort, p_adj) %>%
  print(n = Inf)

write_csv(trend_results, "results/stage_trend_test.csv")

# =============================================================================
# 4. VISUALISATION
# =============================================================================

# Shared colour palette: blue (positive) / red (negative) for partial r
r_palette <- scale_fill_gradient2(
  low = "#B2182B", mid = "white", high = "#2166AC",
  midpoint = 0, limits = c(-0.8, 0.8), name = "Partial r",
  oob = scales::squish
)

# --- 4a. Heatmap: stage × metabolite, coloured by partial r ------------------
# Three panels side by side: tTau | pTau | Aβ42
# Most useful for the supplement — shows full picture at a glance

make_heatmap_panel <- function(bm, pcor_data) {
  dat <- pcor_data %>%
    filter(biomarker == bm) %>%
    mutate(
      met_label   = factor(MET_LABELS[metabolite],
                           levels = rev(MET_LABELS)),
      stage_label = factor(stage_axis_labels[as.character(stage)],
                           levels = stage_axis_labels[STAGE_ORDER]),
      tile_label  = sprintf("%.2f%s", r, ifelse(sig, "*", ""))
    )

  ggplot(dat, aes(x = stage_label, y = met_label, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = tile_label), size = 3.2,
              colour = ifelse(abs(dat$r) > 0.45, "white", "grey20")) +
    r_palette +
    labs(title = bm, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid  = element_blank(),
          axis.text.x = element_text(hjust = 0.5, size = 9),
          plot.title  = element_text(face = "bold", hjust = 0.5))
}

p_heat <- (make_heatmap_panel("TAU", pcor_staged) +
           make_heatmap_panel("PTAU", pcor_staged) +
           make_heatmap_panel("ABETA42", pcor_staged)) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title    = "Stage-stratified polyol–biomarker partial correlations",
    subtitle = "Partial r (Age, Sex, BMI adjusted) | * BH-adj p < 0.05",
    theme    = theme(plot.title    = element_text(face = "bold"),
                     plot.subtitle = element_text(colour = "grey40"))
  ) &
  theme(legend.position = "right")

ggsave("results/figures/staged_heatmap.pdf", p_heat, width = 12, height = 5)
ggsave("results/figures/staged_heatmap.png", p_heat, width = 12, height = 5,
       dpi = 300)

# --- 4b. Dot-CI plot: tau only, all stages, with 95% CIs --------------------
# Main figure for the reviewer response / manuscript supplement.
# Shows uncertainty clearly — important given small NCI and MCI groups.

stage_colours <- c(HC  = "#4E79A7", MCI = "#F28E2B",
                   AD  = "#E15759", VaD = "#76B7B2")

p_dotci <- pcor_staged %>%
  filter(biomarker %in% TAU_VARS) %>%
  mutate(
    met_label   = factor(MET_LABELS[metabolite], levels = rev(MET_LABELS)),
    stage_label = factor(stage_axis_labels[as.character(stage)],
                         levels = stage_axis_labels[STAGE_ORDER]),
    biomarker   = factor(biomarker, levels = TAU_VARS)
  ) %>%
  ggplot(aes(x = r, y = met_label,
             colour = stage, group = stage, alpha = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi),
    height = 0, linewidth = 0.65,
    position = position_dodge(width = 0.65)
  ) +
  geom_point(
    size = 2.6,
    position = position_dodge(width = 0.65)
  ) +
  scale_colour_manual(values = stage_colours, name = "Stage",
                      labels = stage_axis_labels[STAGE_ORDER]) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.28), guide = "none") +
  facet_wrap(~ biomarker, ncol = 2) +
  labs(
    title    = "Within-stage polyol–tau partial correlations",
    subtitle = "95% CI (Fisher z)  |  faded = BH-adj p ≥ 0.05  |  adjusted for Age, Sex, BMI",
    x        = "Partial correlation (r)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "#f2f2f2"),
    strip.text       = element_text(face = "bold")
  )

ggsave("results/figures/staged_dotci.pdf",  p_dotci, width = 9,  height = 6)
ggsave("results/figures/staged_dotci.png",  p_dotci, width = 9,  height = 6,
       dpi = 300)

# --- 4c. Trend line: mean polyol–tau r across stages ------------------------
# Summarises the directional trend across stages using the mean partial r
# across all polyol metabolites. Error bars = SE across metabolites.

trend_viz <- pcor_staged %>%
  filter(biomarker %in% TAU_VARS) %>%
  group_by(stage, biomarker) %>%
  summarise(
    mean_r  = mean(r,  na.rm = TRUE),
    se_r    = sd(r, na.rm = TRUE) / sqrt(sum(!is.na(r))),
    n_group = first(n),
    .groups = "drop"
  ) %>%
  mutate(
    stage_label = factor(stage_axis_labels[as.character(stage)],
                         levels = stage_axis_labels[STAGE_ORDER]),
    biomarker   = factor(biomarker, levels = TAU_VARS)
  )

p_trend <- ggplot(trend_viz,
                  aes(x = stage_label, y = mean_r,
                      colour = biomarker, group = biomarker)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin = mean_r - se_r, ymax = mean_r + se_r),
                width = 0.12, linewidth = 0.7,
                position = position_dodge(0.15)) +
  geom_line(linewidth = 1.0, position = position_dodge(0.15)) +
  geom_point(size = 3.5, position = position_dodge(0.15)) +
  scale_colour_manual(values = c(TAU = "#4E79A7", PTAU = "#F28E2B"),
                      name = "Tau biomarker") +
  scale_y_continuous(
    breaks = seq(-0.2, 0.5, 0.1),
    limits = c(-0.25, 0.55),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title    = "Polyol–tau association across disease stages",
    subtitle = "Mean partial r ± SE across 6 polyol metabolites  |  Age, Sex, BMI adjusted",
    x        = "Diagnostic stage",
    y        = "Mean partial correlation (r)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.subtitle    = element_text(colour = "grey40", size = 9))



ggsave("results/figures/staged_trend.pdf", p_trend, width = 6, height = 5)
ggsave("results/figures/staged_trend.png", p_trend, width = 6, height = 5,
       dpi = 300)

# =============================================================================
# 5. SUPPLEMENTARY TABLES
# =============================================================================

# Table S_staged_pcor: full within-group partial correlation table
supp_pcor <- pcor_staged %>%
  dplyr::mutate(
    ci_str  = sprintf("[%.3f, %.3f]", ci_lo, ci_hi),
    p_label = case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01  ~ "< 0.01",
      TRUE            ~ sprintf("= %.3f", p_value)
    ),
    sig_str = ifelse(sig, "*", ""),
    r       = round(r, 3)
  ) %>%
  dplyr::arrange(biomarker, metabolite, stage) %>%
  dplyr::select(biomarker, stage, n, metabolite, r, ci_str, p_label, sig_str)

write_csv(supp_pcor, "results/staged_supplement_pcor.csv")

# Table S_staged_trend: trend test results
supp_trend <- trend_results %>%
  dplyr::mutate(
    ci_str  = sprintf("[%.3f, %.3f]", conf.low, conf.high),
    p_label = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01  ~ "< 0.01",
      TRUE            ~ sprintf("= %.3f", p.value)
    ),
    sig_str  = ifelse(sig, "*", ""),
    estimate = round(estimate, 3)
  ) %>%
  dplyr::arrange(biomarker, cohort, p_adj) %>%
  dplyr::select(cohort, biomarker, metabolite, estimate, ci_str, p_label, sig_str)

write_csv(supp_trend, "results/staged_supplement_trend.csv")

# Console summary
cat("\n=== Summary for reviewer response text ===\n")
cat("Fill in the bracketed values from your output:\n\n")

pcor_staged %>%
  filter(biomarker %in% TAU_VARS, sig == TRUE) %>%
  count(stage, biomarker, name = "n_sig_metabolites") %>%
  mutate(total = length(POLYOL_VARS),
         summary = sprintf("%s / %d polyols significant in %s (%s)",
                           n_sig_metabolites, total,
                           as.character(stage), biomarker)) %>%
  pull(summary) %>%
  cat(sep = "\n")

cat("\n\nAll files saved to results/ and results/figures/\n")
cat("Figures: staged_heatmap, staged_dotci, staged_trend (PDF + PNG)\n")
cat("Tables:  staged_partial_correlations.csv, stage_trend_test.csv,\n")
cat("         staged_supplement_pcor.csv, staged_supplement_trend.csv\n")
