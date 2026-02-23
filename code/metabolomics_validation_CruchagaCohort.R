################################################################################
############################# ADNI METABOLOMICS ################################
################################################################################

# libraries
library(here)
library(tidyverse)
library(ggpubr)
library(cowplot)

# load clinical data
dx <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/DXSUM_24Dec2025.csv"), sep = ",", header = T)
metadata <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/ADNI_SelectedVars.csv"), sep = ",", header = T)
metadata <- metadata |> 
  dplyr::filter(!is.na(TAU), !is.na(PTAU), !is.na(ABETA42))

dx <- dx[, c("PTID", "RID")]
metadata <- metadata |> 
  dplyr::left_join(dx, by = c("subject_id" = "PTID")) |> 
  dplyr::select(RID, subject_id, VISCODE2 = visit, AGE = entry_age, PTGENDER, DIAGNOSIS, GENOTYPE, ABETA42, PTAU, TAU)

# METABOLOMICS ------------------------------------------------------------------
# load metabolomics 
mets <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/CruchagaLab_CSF_metabolomic_matrix_20230620.csv"), sep = ",", header = T)
mets <- mets |> 
  dplyr::left_join(metadata, by = c("RID", "VISCODE2")) |> 
  dplyr::select(RID, subject_id, EXAMDATE, VISCODE2, AGE, PTGENDER, DIAGNOSIS, GENOTYPE, ABETA42, PTAU, TAU, everything()) |> 
  dplyr::select(-GUSPECID, -Metabolon_Barcode) |> 
  distinct() 

# annotate metabolites
mets_annot <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/ADNI_Cruchaga_lab_ADNI_CSF_metabolomic_info_20_06_2023.csv"), sep = ",", header = T)
mets_annot <- mets_annot |> 
  dplyr::filter(PLOT_NAME %in% c( "glucose", "glucose 6-phosphate", "UDP-glucose", "erythritol", "UDP-galactose",
                                  "galactose 1-phosphate", "mannitol/sorbitol", "xylose", "ribose", "arabitol/xylitol", "myo-inositol",
                                  "ribonate", "gluconate", "threonate", "glycerate")) 
mets_annot <- mets_annot |> 
  dplyr::select(ChemID, PLOT_NAME)

# Filter and check normality for metabolomics
mets <- mets |> 
  dplyr::select(RID, subject_id, EXAMDATE, VISCODE2, AGE, PTGENDER, DIAGNOSIS, GENOTYPE, ABETA42, PTAU, TAU, all_of(mets_annot$ChemID))

qqplots <- sapply(mets[, 12:20], ggqqplot)
cowplot::plot_grid(plotlist = qqplots)

qqplots_log <- sapply(mets[, 12:20], function(x) ggqqplot(log(x)))
cowplot::plot_grid(plotlist = qqplots_log)

mets[, 12:20] <- sapply(mets[, 12:20], log)

### compute ratios glucose/sorbitol and glucose/erythritol
mets$gluc.sorb <- mets$X572 - mets$X100001740
mets$gluc.ery <- mets$X572 - mets$X100000846

# transform tau and AB
mets$ABETA42 <- scale(mets$ABETA42)
mets$TAU <- scale(log(mets$TAU))
mets$PTAU <- scale(log(mets$PTAU))

mets$DIAGNOSIS <- factor(mets$DIAGNOSIS)

# run regressions ----
ptau_assoc <- lapply(mets[, 12:22], function(x){
  mod <- glm(PTAU ~ x + PTGENDER + AGE, data = mets)
  mod.s <- summary(mod)
  mod.s <- as.data.frame(coefficients(mod.s))
  
  out <- c("r2" = mod.s$Estimate[2],
           "std" = mod.s$`Std. Error`[2],
           "p" = mod.s$`Pr(>|t|)`[2])
})

ptau_assoc <- bind_rows(ptau_assoc, .id = "metabolite") |> 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) |> 
  dplyr::left_join(mets_annot, by = c("metabolite" = "ChemID"))

ptau_assoc$PLOT_NAME[11:12] <- c("Glucose / Sorbitol", "Glucose / Erythritol")

ttau_assoc <- lapply(mets[, 12:22], function(x){
  mod <- glm(TAU ~ x + PTGENDER + AGE, data = mets)
  mod.s <- summary(mod)
  mod.s <- as.data.frame(coefficients(mod.s))
  
  out <- c("r2" = mod.s$Estimate[2],
           "std" = mod.s$`Std. Error`[2],
           "p" = mod.s$`Pr(>|t|)`[2])
})

ttau_assoc <- bind_rows(ttau_assoc, .id = "metabolite") |> 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) |> 
  dplyr::left_join(mets_annot, by = c("metabolite" = "ChemID"))

ttau_assoc$PLOT_NAME[11:12] <- c("Glucose / Sorbitol", "Glucose / Erythritol")

AB_assoc <- lapply(mets[, 12:22], function(x){
  mod <- glm(ABETA42 ~ x + PTGENDER + AGE, data = mets)
  mod.s <- summary(mod)
  mod.s <- as.data.frame(coefficients(mod.s))
  
  out <- c("r2" = mod.s$Estimate[2],
           "std" = mod.s$`Std. Error`[2],
           "p" = mod.s$`Pr(>|t|)`[2])
})

AB_assoc <- bind_rows(AB_assoc, .id = "metabolite") |> 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) |> 
  dplyr::left_join(mets_annot, by = c("metabolite" = "ChemID"))

AB_assoc$PLOT_NAME[11:12] <- c("Glucose / Sorbitol", "Glucose / Erythritol")

# results table
bind_rows(ttau_assoc,
          ptau_assoc,
          AB_assoc) |> 
  View()

# plot metabolites
mets$DIAGNOSIS <- case_when(
  mets$DIAGNOSIS == 1 ~ "HC",
  mets$DIAGNOSIS == 2 ~ "MCI",
  mets$DIAGNOSIS == 3 ~ "AD"
) 
mets$DIAGNOSIS <- factor(mets$DIAGNOSIS, levels = c("HC", "MCI", "AD"))

ptau_plots <- list()
for(m in names(mets)[12:22]){
  stats <- ptau_assoc |> 
    dplyr::filter(metabolite == m) |> 
    dplyr::mutate(plot = paste0("r2 = ", round(r2, digits = 2), ", FDR = ", round(fdr, digits = 2)))
  df <- mets |> dplyr::select(PTAU, DIAGNOSIS, met = all_of(m))
  ptau_plots[[m]] <- ggplot(df, aes(PTAU, met)) +
    geom_point(aes(color = DIAGNOSIS)) +
    geom_smooth(method="lm", se = F, color = "black") +
    annotate("text", label = stats$plot, size = 5,
             x = -Inf, y = Inf,
             hjust = -0.05,
             vjust = 1.1) +
    theme_minimal() +
    labs(y = stats$PLOT_NAME, 
         x = "pTau") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = .8))
}

plot_grid(plotlist = ptau_plots)

ab_plots <- list()
for(m in names(mets)[12:22]){
  stats <- AB_assoc |> 
    dplyr::filter(metabolite == m) |> 
    dplyr::mutate(plot = paste0("r2 = ", round(r2, digits = 2), ", FDR = ", round(fdr, digits = 2)))
  df <- mets |> dplyr::select(ABETA42, DIAGNOSIS, met = all_of(m))
  ab_plots[[m]] <- ggplot(df, aes(ABETA42, met)) +
    geom_point(aes(color = DIAGNOSIS)) +
    geom_smooth(method="lm", se = F, color = "black") +
    annotate("text", label = stats$plot, size = 5,
             x = -Inf, y = Inf,
             hjust = -0.05,
             vjust = 1.1) +
    theme_minimal() +
    labs(y = stats$PLOT_NAME, 
         x = "pTau") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = .8))
}

plot_grid(plotlist = ab_plots)











#### PROTEOMICS ----------------------------------------------------------------
# load & filter proteomics
prots <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"), sep = ",", header = T)
prots_anno <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv"), sep = ",", header = T)

polyol <- read.table(here("data/ExternalCohorts/CruchagaLab_CSF_metabolomics/POLYOLS_NETWORK_EXPANDED.txt"), header = T, sep  = "\t")
prots_anno <- prots_anno |> 
  dplyr::filter(EntrezGeneSymbol %in% unique(polyol$gene_symbol))

prots <- prots |> 
  dplyr::select(RID, EXAMDATE, VISCODE2, PlateId, all_of(prots_anno$Analytes)) |> 
  dplyr::left_join(metadata, by = c("RID", "VISCODE2")) |> 
  distinct() |> 
  dplyr::select(RID, EXAMDATE, VISCODE2, AGE, PTGENDER, DIAGNOSIS, GENOTYPE, ABETA42, PTAU, TAU, everything()) 

# Filter and check normality for metabolomics
qqplots <- sapply(prots[, 12:51], ggqqplot)
cowplot::plot_grid(plotlist = qqplots)

qqplots_log <- sapply(prots[, 12:51], function(x) ggqqplot(log(x)))
cowplot::plot_grid(plotlist = qqplots_log)

prots[, 12:51] <- sapply(prots[, 12:51], log)

# transform tau and AB
prots$ABETA42 <- scale(prots$ABETA42)
prots$TAU <- scale(log(prots$TAU))
prots$PTAU <- scale(log(prots$PTAU))

prots$DIAGNOSIS <- factor(prots$DIAGNOSIS)


# run regressions ----
ptau_assoc_px <- lapply(prots[, 12:51], function(x){
  mod <- glm(PTAU ~ x + PTGENDER + AGE, data = prots)
  mod.s <- summary(mod)
  mod.s <- as.data.frame(coefficients(mod.s))
  
  out <- c("r2" = mod.s$Estimate[2],
           "std" = mod.s$`Std. Error`[2],
           "p" = mod.s$`Pr(>|t|)`[2])
})

ptau_assoc_px <- bind_rows(ptau_assoc_px, .id = "Analytes") |> 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) |> 
  dplyr::left_join(prots_anno, by = "Analytes") |> 
  dplyr::select(Analytes, Protein = EntrezGeneSymbol, r2, std, p, fdr) |> 
  arrange(fdr)

ttau_assoc_px <- lapply(prots[, 12:51], function(x){
  mod <- glm(TAU ~ x + PTGENDER + AGE, data = prots)
  mod.s <- summary(mod)
  mod.s <- as.data.frame(coefficients(mod.s))
  
  out <- c("r2" = mod.s$Estimate[2],
           "std" = mod.s$`Std. Error`[2],
           "p" = mod.s$`Pr(>|t|)`[2])
})

ttau_assoc_px <- bind_rows(ttau_assoc_px, .id = "Analytes") |> 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) |> 
  dplyr::left_join(prots_anno, by = "Analytes") |> 
  dplyr::select(Analytes, Protein = EntrezGeneSymbol, r2, std, p, fdr) |> 
  arrange(fdr)

AB_assoc_px <- lapply(prots[, 12:51], function(x){
  mod <- glm(ABETA42 ~ x + PTGENDER + AGE, data = prots)
  mod.s <- summary(mod)
  mod.s <- as.data.frame(coefficients(mod.s))
  
  out <- c("r2" = mod.s$Estimate[2],
           "std" = mod.s$`Std. Error`[2],
           "p" = mod.s$`Pr(>|t|)`[2])
})

AB_assoc_px <- bind_rows(AB_assoc_px, .id = "Analytes") |> 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) |> 
  dplyr::left_join(prots_anno, by = "Analytes") |> 
  dplyr::select(Analytes, Protein = EntrezGeneSymbol, r2, std, p, fdr) |> 
  arrange(fdr)


# results table
bind_rows(ttau_assoc_px,
          ptau_assoc_px,
          AB_assoc_px) |> 
  View()

# plot metabolites
prots$DIAGNOSIS <- case_when(
  prots$DIAGNOSIS == 1 ~ "HC",
  prots$DIAGNOSIS == 2 ~ "MCI",
  prots$DIAGNOSIS == 3 ~ "AD"
) 
prots$DIAGNOSIS <- factor(prots$DIAGNOSIS, levels = c("HC", "MCI", "AD"))

ptau_px_plots <- list()
for(m in names(prots)[12:51]){
  stats <- ptau_assoc_px |> 
    dplyr::filter(Analytes == m) |> 
    dplyr::mutate(plot = paste0("r2 = ", round(r2, digits = 2), ", FDR = ", round(fdr, digits = 2)))
  df <- prots |> dplyr::select(PTAU, DIAGNOSIS, met = all_of(m))
  ptau_px_plots[[m]] <- ggplot(df, aes(PTAU, met)) +
    geom_point(aes(color = DIAGNOSIS)) +
    geom_smooth(method="lm", se = F, color = "black") +
    annotate("text", label = stats$plot, size = 5,
             x = -Inf, y = Inf,
             hjust = -0.05,
             vjust = 1.1) +
    theme_minimal() +
    labs(y = stats$Protein, 
         x = "pTau") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = .8))
}

plot_grid(plotlist = ptau_px_plots)

ab_px_plots <- list()
for(m in names(prots)[11:50]){
  stats <- AB_assoc_px |> 
    dplyr::filter(Analytes == m) |> 
    dplyr::mutate(plot = paste0("r2 = ", round(r2, digits = 2), ", FDR = ", round(fdr, digits = 2)))
  df <- prots |> dplyr::select(ABETA42, DIAGNOSIS, met = all_of(m))
  ab_px_plots[[m]] <- ggplot(df, aes(ABETA42, met)) +
    geom_point(aes(color = DIAGNOSIS)) +
    geom_smooth(method="lm", se = F, color = "black") +
    annotate("text", label = stats$plot, size = 5,
             x = -Inf, y = Inf,
             hjust = -0.05,
             vjust = 1.1) +
    theme_minimal() +
    labs(y = stats$Protein, 
         x = "pTau") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = .8))
}

plot_grid(plotlist = ab_px_plots)

#### plot polyol pathway proteins
ptau_assoc_px[grepl("AKR1B1", ptau_assoc_px$Protein), ]
ptau_assoc_px[grepl("SORD", ptau_assoc_px$Protein), ]

df <- prots |> 
  dplyr::select(TAU, PTAU, ABETA42, DIAGNOSIS, X16606.85, X15447.45) |> 
  dplyr::filter(complete.cases(across(everything()))) |> 
  dplyr::rename(AKR1B1 = X16606.85, SORD = X15447.45) |> 
  dplyr::mutate(SORD = log(SORD)) |> 
  dplyr::mutate(AKR1B1 = log(AKR1B1))

png("C:/Users/MGAR0052/Desktop/Papers/Fig3.png", res = 300, width = 3500, height = 1350)
cowplot::plot_grid(ggpubr::ggarrange(ggplot(df, aes(TAU, AKR1B1)) + geom_point() + geom_smooth(method = "lm", se = F) +
                                       theme_bw() + labs(x = "tTau (scaled, log-transformed)",
                                                         y = "AKR1B1 (log-transformed)") + 
                                       stat_cor(geom = "label", fill = "white", alpha = 0.8),
                                     ggplot(df, aes(PTAU, AKR1B1)) + geom_point() + geom_smooth(method = "lm", se = F) +
                                       theme_bw() + labs(x = "pTau (scaled, log-transformed)",
                                                         y = "AKR1B1 (log-transformed)") + 
                                       stat_cor(geom = "label", fill = "white", alpha = 0.8),
                                     ggplot(df, aes(ABETA42, AKR1B1)) + geom_point() + geom_smooth(method = "lm", se = F) +
                                       theme_bw() + labs(x = "Aβ42 (scaled)",
                                                         y = "AKR1B1 (log-transformed)") + 
                                       stat_cor(geom = "label", fill = "white", alpha = 0.8), 
                                     ncol = 3),
                   ggpubr::ggarrange(ggplot(df, aes(TAU, SORD)) + geom_point() + geom_smooth(method = "lm", se = F) +
                                       theme_bw() + labs(x = "tTau (scaled, log-transformed)",
                                                         y = "SORD (log-transformed)") + 
                                       stat_cor(geom = "label", fill = "white", alpha = 0.8),
                                     ggplot(df, aes(PTAU, SORD)) + geom_point() + geom_smooth(method = "lm", se = F) +
                                       theme_bw() + labs(x = "pTau (scaled, log-transformed)",
                                                         y = "SORD (log-transformed)") + 
                                       stat_cor(geom = "label", fill = "white", alpha = 0.8),
                                     ggplot(df, aes(ABETA42, SORD)) + geom_point() + geom_smooth(method = "lm", se = F) +
                                       theme_bw() + labs(x = "Aβ42 (scaled)",
                                                         y = "SORD (log-transformed)") + 
                                       stat_cor(geom = "label", fill = "white", alpha = 0.8), 
                                     ncol = 3),
                   ncol = 1)
dev.off()

#### Test overlap of differentially regulated proteins in SomaScan with polyol network
library(clusterProfiler)

background <- read.table("ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv", sep = ",", header = T)
background <- unique(background$EntrezGeneSymbol)

netw_proteins <- unique(polyol$gene_symbol)

ptau_sig <- ptau_assoc_px |> dplyr::filter(fdr <= 0.05) |> pull(Protein)
ttau_sig <- ttau_assoc_px |> dplyr::filter(fdr <= 0.05) |> pull(Protein)
ab_sig <- AB_assoc_px |> dplyr::filter(fdr <= 0.05) |> pull(Protein)


run_enrichment_cp <- function(hits, network, background) {
  
  # Ensure everything is inside the universe
  hits       <- intersect(hits, background)
  network    <- intersect(network, background)
  background <- unique(background)
  
  # TERM2GENE: one pathway = your network
  term2gene <- tibble(
    term = "Polyol_Network",
    gene = network
  )
  
  # Run enrichment
  enr <- enricher(
    gene          = hits,
    TERM2GENE     = term2gene,
    universe      = background,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  
  # If no enrichment detected
  if (is.null(enr) || nrow(enr@result) == 0) {
    return(
      tibble(
        overlap = 0,
        network_size = length(network),
        hit_size = length(hits),
        odds_ratio = NA,
        conf_low = NA,
        conf_high = NA,
        p_value = 1
      )
    )
  }
  
  res <- enr@result[1, ]
  return(res)
  
}

prots_validation <- list(
  ttau = ttau_sig,
  ptau = ptau_sig, 
  ab = ab_sig
)

prots_validation <- lapply(prots_validation, function(x){
  run_enrichment_cp(hits = x,
                    network = netw_proteins,
                    background = background)
})

bind_rows(prots_validation, .id = "biomarker")

run_enrichment <- function(hits, network, background) { 
  # Ensure all sets are within background 
  hits <- intersect(hits, background)
  network <- intersect(network, background)
  a <- length(intersect(hits, network)) # network & significant 
  b <- length(setdiff(network, hits)) # network & not significant 
  c <- length(setdiff(hits, network)) # non-network & significant 
  d <- length(setdiff(background, union(hits, network))) # non-network & not significant 
  
  contingency <- matrix(c(a, b, c, d), 
                        nrow = 2, 
                        dimnames = list( 
                          Network = c("Yes", "No"), 
                          Significant = c("Yes", "No") 
                        )) 
  
  test <- fisher.test(contingency) 
  tibble( 
    overlap = a,
    network_size = length(network),
    hit_size = length(hits), 
    odds_ratio = unname(test$estimate), 
    conf_low = test$conf.int[1], 
    conf_high = test$conf.int[2], 
    p_value = test$p.value 
  ) 
}


prots_validation <- list(
  ttau = ttau_sig,
  ptau = ptau_sig, 
  ab = ab_sig
)

prots_validation <- lapply(prots_validation, function(x){
  run_enrichment(hits = x,
                    network = netw_proteins,
                    background = background)
})

bind_rows(prots_validation, .id = "biomarker")

