# ------------------------------------------------------
#  Bayesian differential abundance analysis (metabolomics)
#  Requirements:
#  - your_data: samples x metabolites (numeric)
#  - class_vector: binary vector (0 = class1, 1 = class2)
# ------------------------------------------------------

# Install if needed
if (!require(brms)) install.packages("brms", dependencies = TRUE)
if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE)
if (!require(furrr)) install.packages("furrr", dependencies = TRUE)

library(brms)
library(tidyverse)
library(furrr)
library(here)
library(dplyr)

# ------------------------------------------------------
# 1. Load and prepare data
# ------------------------------------------------------
plasma <- read.table(file.path(here(), "Plasma_CRMNnormalized_ratios_UNtargeted.txt"), header =T)
row.names(plasma) <- plasma$sample
comparisons <- list("HC_AD" = c("HC", "AD"), 
                    "HC_MCI" = c("HC", "MCI"),
                    "HC_MCI_AD" = c("HC", "MCI", "AD"))
metadata <- plasma[, 1:20]
metabolites <- plasma[, 21:ncol(plasma)]

csf <- read.table(file.path(here(), "CSF_CRMNnormalized_ratios_UNtargeted.txt"), header =T)
row.names(csf) <- csf$sample
comparison <- c("HC", "VaD")
metadata <- csf[, 1:20]
metabolites <- csf[, 21:ncol(csf)]

csf_targeted <- read.table(here("CSF_targeted_logtransformed.txt"), header = T)
csf_targeted$fullClass <- gsub("No cognitive impairment", "HC", csf_targeted$fullClass)
row.names(csf_targeted) <- csf_targeted$sample
comparisons <- list("HC_AD" = c("HC", "AD"), 
                    "HC_MCI" = c("HC", "MCI"),
                    "HC_MCI_AD" = c("HC", "MCI", "AD"),
                    "HC_VaD" = c("HC", "VaD"))
metadata <- csf_targeted[, 1:23]
metabolites <- csf_targeted[, 24:ncol(csf_targeted)]

plasma_targeted <- read.table(here("Plasma_targeted_logtransformed.txt"), header = T)
row.names(plasma_targeted) <- plasma_targeted$sample
comparisons <- list("HC_AD" = c("HC", "AD"), 
                    "HC_MCI" = c("HC", "MCI"),
                    "HC_MCI_AD" = c("HC", "MCI", "AD"),
                    "HC_VaD" = c("HC", "VaD"))
metadata <- plasma_targeted[, 1:23]
metabolites <- plasma_targeted[, 24:ncol(plasma_targeted)]


for(i in 1:length(comparisons)){
  comparison <- comparisons[[i]]
  mtdt <- metadata %>% mutate(fullClass = as.character(fullClass)) %>% filter(fullClass %in% comparison, )
  mets <- metabolites[row.names(metabolites) %in% mtdt$sample, ]
  
class_vector <- factor(ifelse(mtdt$fullClass == "HC", 0, 1))
# Ensure data is numeric and scaled

# ------------------------------------------------------
# 2. Define a function to fit one metabolite
# ------------------------------------------------------
## run our own metabolites
fit_metabolite_bayes <- function(response, class_vector) {
  df <- data.frame(y = response, class = as.factor(class_vector), age = mtdt$AgeAtVisit, gender = mtdt$Gender, statins = mtdt$Statins)
  
  # Bayesian model for a single metabolite
  model <- brm(
    formula = y ~ class + age + gender + statins,
    data = df,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "b"),           # effect size (β)
      prior(normal(0, 5), class = "Intercept"),   # baseline mean (μ)
      prior(student_t(3, 0, 1), class = "sigma")  # variance (σ)
    ),
    iter = 2000, chains = 2, cores = 2, silent = TRUE, refresh = 0
  )
  
  post <- posterior_samples(model)
  eff <- post$b_class1  # coefficient for class effect
  
  summary <- tibble(
    mean_beta = mean(eff),
    lower_95 = quantile(eff, 0.025),
    upper_95 = quantile(eff, 0.975),
    p_pos = mean(eff > 0),
    p_neg = mean(eff < 0)
  )
  
  return(summary)
}

# ------------------------------------------------------
# 3. Run Bayesian models for all metabolites
# ------------------------------------------------------
results <- map_dfr(
  mets %>% as.data.frame() %>% as.list(),
  ~ fit_metabolite_bayes(.x, class_vector),
  .id = "Metabolite"
)

# ------------------------------------------------------
# 4. Interpret results
# ------------------------------------------------------
results <- results %>%
  mutate(
    significant = if_else(p_pos > 0.95 | p_neg > 0.95, TRUE, FALSE),
    direction = case_when(
      p_pos > 0.95 ~ "Higher in Class 2",
      p_neg > 0.95 ~ "Lower in Class 2",
      TRUE ~ "No difference"
    )
  )


# ------------------------------------------------------
# 6. Save / visualize results
# ------------------------------------------------------
#write.csv(results, "HCvsVaD_bayesian_results.csv", row.names = FALSE)
#write.csv(results, paste0("Plasma_", paste(comparison, collapse="_"), "_bayesian_results.csv"), row.names = FALSE)
#write.csv(results, paste0("CSF_TARGETED_", paste(comparison, collapse="_"), "_bayesian_results.csv"), row.names = FALSE)
write.csv(results, paste0("Plasma_TARGETED_", paste(comparison, collapse="_"), "_bayesian_results.csv"), row.names = FALSE)

}

