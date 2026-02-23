################################################################################
#################### Polyol and Tau focused analysis ###########################
################################################################################

# set up -------------------
## libraries
library(here)
library(tidyverse)
library(ggnewscale)
library(patchwork)
library(clusterProfiler)
library(visNetwork)
library(stringdist)
library(dplyr)
library(rlang)
library(RSQLite)
library(STRINGdb)
library(igraph)
library(org.Hs.eg.db)
library(aPEAR)

## directories
dir_data <- here("data")
if (!dir.exists(dir_data)) dir.create(dir_data)
dir_res <- here("results")
if (!dir.exists(dir_res)) dir.create(dir_res)
dir_code <- here("code")
if (!dir.exists(dir_code)) dir.create(dir_code)

## functions --------------
## Best matches between two columns (HMDB annotation)
get_best_match <- function(df, query_col, ref1_col, ref2_col) {
  query_sym <- sym(query_col)

  df |>
    rowwise() |>
    mutate(
      dist_ref1 = if_else(get(query_col) == get(ref1_col), 0,
        stringdist::stringdist(get(query_col), get(ref1_col), method = "lv")
      ),
      dist_ref2 = if_else(get(query_col) == get(ref2_col), 0,
        stringdist::stringdist(get(query_col), get(ref2_col), method = "lv")
      ),
      combined_score = dist_ref1 + dist_ref2
    ) |>
    ungroup() |>
    arrange(!!query_sym, combined_score) |> # lowest combined score wins
    distinct(metabolite, .keep_all = TRUE) # keep only the first entry, as its the lowest score
}

############# function for plotting half-tiles
make_half_triangles <- function(df, level1, level2) {
  df |>
    mutate(
      x = as.numeric(factor(comparison)),
      y = as.numeric(factor(metabolite)),
      triangle = case_when(
        tissue == level1 ~ "lower",
        tissue == level2 ~ "upper"
      )
    ) |>
    rowwise() |>
    do({
      x <- .$x
      y <- .$y

      if (.$triangle == "lower") {
        data.frame(
          x = c(x - 0.5, x + 0.5, x - 0.5),
          y = c(y - 0.5, y - 0.5, y + 0.5),
          tissue = .$tissue,
          beta = .$beta,
          comparison = .$comparison,
          metabolite = .$metabolite
        )
      } else {
        data.frame(
          x = c(x + 0.5, x + 0.5, x - 0.5),
          y = c(y + 0.5, y - 0.5, y + 0.5),
          tissue = .$tissue,
          beta = .$beta,
          comparison = .$comparison,
          metabolite = .$metabolite
        )
      }
    })
}


########### extract EnrichmentNetworks function
extract_network <- function(enrichment, enrichment_networks) {
  # Extract data layers from the enrichment_network
  nodes <- enrichment_networks$plot$layers$geom_point$data

  # Merge enrichment DF to nodes DF using the 'process' as a column name
  nodes_df <- dplyr::left_join(
    nodes,
    enrichment |> dplyr::select(Description), # Select 'Description' and the column from process
    by = c("ID" = "Description")
  )

  edges_df <- enrichment_networks$plot$layers$geom_link0$data # Corrected here
  labels_df <- enrichment_networks$plot$layers$geom_text$data # Corrected here


  ## Return the result list containing modified nodes, edges and label data for each process
  result <- list(
    nodes = nodes_df,
    edges = edges_df,
    labels = labels_df
  )
  return(result)
}


## load data ------
# diff abundance
met.res <- read.table(file.path(dir_res, "AllRelevantHits_metabolomicsFINAL.txt"), header = T, sep = "\t")

# annotation table
hmdb.annot <- read.delim(file.path(dir_res, "HMDB_ANNOTATIONS_ALL.txt"), header = T)

# HMDB pathways
hmdb_ontology <- read.table(here("all_metabolite_pathways.csv"), sep = ",", header = T)
hmdb_ontology <- hmdb_ontology |> dplyr::filter(pathway_category == "kegg")

## filter results for biomarkers
met.bmks <- met.res |>
  dplyr::filter(Comparison %in% c("betaAmyloid", "pTau", "tTau"))

# we work with sig.mets.pattern
## add HMDB codes to metabolites ----
metabolites.to.map <- unique(met.bmks$Metabolite.name)
#### split ratios
metabolites.to.map <- unique(unlist(sapply(metabolites.to.map, function(x) {
  if (grepl(" / ", x)) {
    strsplit(x, " / ")[[1]]
  } else {
    x
  }
})))

## map to HMDB
hmdb_out <- list()
for (met in metabolites.to.map) {
  ## match to query name
  hmdb.out <- hmdb.annot[which(hmdb.annot$query_term == met), ]
  if (nrow(hmdb.out) >= 1) {
    hmdb_out[[met]] <- hmdb.out
  } else {
    hmdb.out <- hmdb.annot[grepl(met, hmdb.annot$query_term, ignore.case = TRUE) |
      grepl(met, hmdb.annot$name, ignore.case = TRUE), ]
    hmdb_out[[met]] <- hmdb.out
  }
}
hmdb_out <- do.call(rbind, hmdb_out) |>
  tibble::rownames_to_column("metabolite") |>
  tidyr::separate(metabolite, sep = "\\.", into = "metabolite")

## clean the annotation DF a bit


hmdb_annot <- get_best_match(hmdb_out, "metabolite", "query_term", "name") |> dplyr::select(metabolite, query_term, hmdb_id, name)
hmdb_annot[which(hmdb_annot$query_term == "Oxalic acid"), ]$name <- "Oxalic acid"


## metabolite-protein interaction ----
con <- OmnipathR::metalinksdb_sqlite()

query <- "
SELECT
  ? AS hmdb_id,
  e.uniprot,
  p.gene_symbol,
  e.type,
  e.mor
FROM edges e
JOIN proteins p ON e.uniprot = p.uniprot
WHERE e.hmdb = ?;
"


met_prot <- lapply(hmdb_annot$hmdb_id, function(x) {
  out <- dbGetQuery(con, query, params = list(x, x))
})
names(met_prot) <- hmdb_annot$query_term

met_prot <- do.call(rbind, met_prot) %>%
  distinct() %>%
  tibble::rownames_to_column("metabolite") %>%
  tidyr::separate(metabolite, sep = "\\.", into = "metabolite") %>%
  dplyr::select(hmdb_id, uniprot, gene_symbol, metabolite, type, mor)


cat("Number of proteins per metabolite: \n\n")
table(met_prot$metabolite)

## create original network
met_prot_backup <- met_prot
met_prot$uniprot <- NULL

g <- graph_from_data_frame(met_prot)

## extend with STRING or IntAct ----
met_prot <- met_prot_backup

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  network_type = "full"
)

met_prot_string <- string_db$map(
  met_prot,
  "uniprot"
)

met_prot_string <- met_prot_string[complete.cases(met_prot_string), ]
ppi_edges <- string_db$get_interactions(met_prot_string$STRING_id)

# Reverse map STRING -> UniProt
string_to_uniprot <- met_prot_string[, c("STRING_id", "gene_symbol")]
ppi_edges_merged <- merge(ppi_edges,
  string_to_uniprot,
  by.x = "from", by.y = "STRING_id", all.x = TRUE
)
ppi_edges_merged <- merge(ppi_edges_merged,
  string_to_uniprot,
  by.x = "to", by.y = "STRING_id", all.x = TRUE
)
ppi_edges_merged <- ppi_edges_merged %>% distinct()

# Keep only mapped pairs
ppi_edges_final <- na.omit(ppi_edges_merged[, c("gene_symbol.x", "gene_symbol.y", "combined_score")])
ppi_edges_final <- ppi_edges_final %>% dplyr::rename(confidence = combined_score)
ppi_edges_final$confidence <- as.numeric(ppi_edges_final$confidence)

# Build igraph of new PPIs
ppi_graph <- graph_from_data_frame(ppi_edges_final, directed = FALSE)

# Combine with your existing network
ppi_graph <- as_directed(ppi_graph, mode = "mutual")
g_expanded <- igraph::union(g, ppi_graph)

## create edges_df dataframe
edges_df <- data.frame(
  from = ends(g_expanded, E(g_expanded))[, 1],
  to   = ends(g_expanded, E(g_expanded))[, 2]
)

# create nodes data.frame
nodes_df <- data.frame(id = V(g_expanded)$name, label = V(g_expanded)$name)
nodes_df <- nodes_df %>%
  dplyr::left_join(hmdb_annot, by = c("label" = "hmdb_id")) %>%
  dplyr::select(id, label, name) %>%
  dplyr::mutate(label = ifelse(grepl("HMDB", label), name, label)) %>%
  dplyr::select(id, label)

nodes_df$type <- ifelse(grepl("HMDB", nodes_df$id), "Metabolite", "MetaLinksDB") ## identify node type based on its origin

## enrichments - metabolites =========
mets_enr <- met.bmks |> 
  dplyr::left_join(hmdb_annot, by = c("Metabolite.name" = "query_term"))

df <- mets_enr |> 
  dplyr::filter(Comparison == "pTau", Tissue == "CSF") |> 
  dplyr::mutate(enr_score = -log10(pvalue + 0.00001) * sign(Mean.beta)) |>
  dplyr::arrange(-enr_score) |> 
  pull(hmdb_id) |> unique() |> 
  clusterProfiler::enricher(TERM2GENE = hmdb_ontology[, c("pathway_name", "sourceId")],
                              minGSSize = 5, pvalueCutoff = 0.5)
df <- df@result


## enrichments - whole network ==========
go_bp <- enrichGO(unique(met_prot$gene_symbol),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)
enrichmentNetwork(go_bp@result)

### POLYOL(s) NETWORK -----
# Polyol pathway–associated metabolite selection
# Polyol pathway–associated metabolites were defined to include (i) aldose substrates of aldose reductase, (ii) their corresponding polyol (sugar alcohol) products, and (iii) immediate downstream oxidation products arising from polyol–sugar interconversion. Aldose substrates included hexoses and pentoses known to enter the polyol pathway (e.g. glucose, galactose, xylose, ribose, and rhamnose). Polyols comprised sorbitol, arabitol, erythritol, threitol, and related pentitols and hexitols. In addition, selected sugar acids (e.g. ribonic acid, gluconic acid, threonic acid, and glyceric acid) were retained as downstream oxidation products of polyol-related sugars, reflecting redox-linked flux through the pathway. Metabolites not directly connected to aldose reduction or polyol oxidation (e.g. disaccharides, glycolytic intermediates, or TCA cycle metabolites) were excluded to maintain pathway specificity.
# Optional shorter version (if space is tight)
# Polyol pathway–associated metabolites were defined as aldose substrates, their reduced polyol products, and immediate oxidation products reflecting redox flux through the polyol pathway. Only metabolites directly linked to aldose reductase or polyol-associated oxidation were included.

##### select metabolites ----
polyols <- list(
  Upstream_substrates = c(
    "d-Glucose",
    "d-Glucose / Erythritol",
    "d-Galactose",
    "d-Glucose / Sorbitol",
    "D-(-)-Xylose",
    "D-(+)-Xylose",
    "D-(-)-Ribofuranose",
    "D-(-)-Rhamnose"
  ),
  Polyols = c(
    "Sorbitol",
    "L-(-)-Arabitol",
    "D-Threitol",
    "meso-Erythritol",
    "3-Deoxy-erythro-pentitol",
    "1,5-Anhydrohexitol",
    "Myo-inositol"
  ),
  Downstream_alternative_substrates = c(
    "Ribonic acid",
    "D-Gluconic acid",
    "L-Threonic acid",
    "Glyceric acid"
  )
)

lookup <- purrr::imap_dfr(polyols, ~ tibble(
  metabolite = .x,
  group = .y
))


met.polyols <- met.bmks |>
  dplyr::filter(Metabolite.name %in% unlist(polyols))

## add HMDB codes to metabolites ======
metabolites.to.map <- unique(met.polyols$Metabolite.name)
#### split ratios
metabolites.to.map <- unique(unlist(sapply(metabolites.to.map, function(x) {
  if (grepl(" / ", x)) {
    strsplit(x, " / ")[[1]]
  } else {
    x
  }
})))

## map to HMDB
hmdb_out <- list()
for (met in metabolites.to.map) {
  ## match to query name
  hmdb.out <- hmdb.annot[which(hmdb.annot$query_term == met), ]
  if (nrow(hmdb.out) >= 1) {
    hmdb_out[[met]] <- hmdb.out
  } else {
    hmdb.out <- hmdb.annot[grepl(met, hmdb.annot$query_term, ignore.case = TRUE) |
      grepl(met, hmdb.annot$name, ignore.case = TRUE), ]
    hmdb_out[[met]] <- hmdb.out
  }
}
hmdb_out <- do.call(rbind, hmdb_out) |>
  tibble::rownames_to_column("metabolite") |>
  tidyr::separate(metabolite, sep = "\\.", into = "metabolite")

## clean the annotation DF a bit

hmdb_annot <- get_best_match(hmdb_out, "metabolite", "query_term", "name") |> dplyr::select(metabolite, query_term, hmdb_id, name)
hmdb_annot[which(hmdb_annot$query_term == "Oxalic acid"), ]$name <- "Oxalic acid"


## metabolite-protein interaction =======
con <- OmnipathR::metalinksdb_sqlite()

met_polyols <- lapply(hmdb_annot$hmdb_id, function(x) {
  out <- dbGetQuery(con, query, params = list(x, x))
})
names(met_polyols) <- hmdb_annot$query_term

met_polyols <- do.call(rbind, met_polyols) %>%
  distinct() %>%
  tibble::rownames_to_column("metabolite") %>%
  tidyr::separate(metabolite, sep = "\\.", into = "metabolite") %>%
  dplyr::select(hmdb_id, uniprot, gene_symbol, metabolite, type, mor)

met_polyols_string <- string_db$map(
  met_polyols,
  "uniprot"
)
met_polyols_string <- met_polyols_string[!is.na(met_polyols_string$STRING_id), ]


met_polyols$uniprot <- NULL

## plot network ======
g <- graph_from_data_frame(met_polyols)

# Convert igraph → visNetwork
## create edges_df dataframe
edges_df <- data.frame(
  from = ends(g, E(g))[, 1],
  to   = ends(g, E(g))[, 2]
)

# create nodes data.frame
nodes_df <- data.frame(id = V(g)$name, label = V(g)$name)
nodes_df <- nodes_df %>%
  dplyr::left_join(hmdb_annot, by = c("label" = "hmdb_id")) %>%
  dplyr::select(id, label, name) %>%
  dplyr::mutate(label = ifelse(grepl("HMDB", label), name, label)) %>%
  dplyr::select(id, label)

## enrichment - metabolites =====
clusterProfiler::enricher(met_polyols |> pull(hmdb_id) |> unique(),
                          TERM2GENE = hmdb_ontology[, c("pathway_name", "sourceId")])

as.data.frame(clusterProfiler::enrichKEGG(AnnotationDbi::select(org.Hs.eg.db, 
                                                  keys = met_polyols |> pull(gene_symbol) |> unique(), 
                                                  columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL") |> 
                              pull(ENTREZID) |> unique(), 
                            organism = "hsa", keyType = "ncbi-geneid")) |> write.table(file.path(dir_res, "polyol_network_KEGG_enrichments_genes.txt"), sep = "\t", quote = F, row.names = F)

## enrichment - BP =====
go_bp <- enrichGO(unique(met_polyols$gene_symbol),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

go_bp <- go_bp@result
write.table(go_bp, file.path(dir_res, "polyol_network_GOBP_enrichments_genes.txt"), sep = "\t", quote = F, row.names = F)
go_bp |>
  filter(grepl("polyol", Description)) |>
  ggplot(aes(FoldEnrichment, reorder(Description, FoldEnrichment), fill = -log10(p.adjust))) +
  geom_col() +
  theme_minimal() +
  labs(x = "Fold Enrichment",
       y = "Pathways")
  

aPEAR::enrichmentNetwork(go_bp)

## Extend to STRING ====
ppi_edges <- string_db$get_interactions(met_polyols_string$STRING_id) |> dplyr::filter(combined_score >= 700)
# Reverse map STRING -> UniProt
string_to_uniprot <- met_polyols_string[, c("STRING_id", "gene_symbol")]
ppi_edges_merged <- merge(ppi_edges,
  string_to_uniprot,
  by.x = "from", by.y = "STRING_id", all.x = TRUE
)
ppi_edges_merged <- merge(ppi_edges_merged,
  string_to_uniprot,
  by.x = "to", by.y = "STRING_id", all.x = TRUE
)
ppi_edges_merged <- ppi_edges_merged %>% distinct()

# Keep only mapped pairs
ppi_edges_final <- na.omit(ppi_edges_merged[, c("gene_symbol.x", "gene_symbol.y", "combined_score")])
ppi_edges_final <- ppi_edges_final %>% dplyr::rename(confidence = combined_score)
ppi_edges_final$confidence <- as.numeric(ppi_edges_final$confidence)

# Build igraph of new PPIs
ppi_graph <- graph_from_data_frame(ppi_edges_final, directed = FALSE)

# Combine with your existing network
ppi_graph <- as_directed(ppi_graph, mode = "mutual")
g_expanded <- igraph::union(g, ppi_graph)

## create edges_df dataframe
edges_df <- data.frame(
  from = ends(g_expanded, E(g_expanded))[, 1],
  to   = ends(g_expanded, E(g_expanded))[, 2]
)
# create nodes data.frame
nodes_df <- data.frame(id = V(g_expanded)$name, label = V(g_expanded)$name)
nodes_df <- nodes_df %>%
  dplyr::left_join(hmdb_annot, by = c("label" = "hmdb_id")) %>%
  dplyr::select(id, label, name) %>%
  dplyr::mutate(label = ifelse(grepl("HMDB", label), name, label)) %>%
  dplyr::select(id, label)

nodes_df$type <- ifelse(grepl("HMDB", nodes_df$id), "Metabolite", "MetaLinksDB") ## identify node type based on its origin

write.table(met_polyols_string, file.path(dir_res, "POLYOLS_NETWORK_EXPANDED.txt"), sep = "\t", quote = F, row.names = F)

## enrichments - whole network ==========
net_polyols <- as_data_frame(g_expanded)
net_polyols <- net_polyols |>
  dplyr::left_join(met.bmks, by = c("metabolite" = "Metabolite.name")) |>
  dplyr::select(from, to, metabolite, type, mor, beta = Mean.beta, pvalue, comparison = Comparison, tissue = Tissue) |>
  dplyr::filter(!is.na(beta)) |>
  dplyr::mutate(interaction = case_when(
    beta > 0 & mor == 1 ~ "Increased_Activator",
    beta > 0 & mor == -1 ~ "Increased_Inhibitor",
    beta < 0 & mor == 1 ~ "Decreased_Activator",
    beta < 0 & mor == -1 ~ "Decreased_Inhibitor",
  ))

net_polyols |>
  dplyr::filter(tissue == "CSF", comparison == "pTau") |>
  split(~interaction) |>
  purrr::map(\(net) as.data.frame(enrichGO(unique(net$to),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )))


## metabolite-wise ----
metabolite_BP <- net_polyols |>
  dplyr::filter(tissue == "CSF") |>
  split(~metabolite) |>
  purrr::map(\(net) as.data.frame(enrichGO(unique(net$to),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )))


metabolite_BP_nets <- lapply(metabolite_BP, function(BP) {
  enrichmentNetwork(BP, plotOnly = F)
})
bp_enr_nets <- bp_enr |>
  purrr::map(\(net) (enrichmentNetwork(net, plotOnly = FALSE, minClusterSize = 10)))

metabolite_BP_nets_xtr <- lapply(names(metabolite_BP), function(met){
  enrichment <- metabolite_BP[[met]]
  enrichment_networks <- metabolite_BP_nets[[met]]
  clusters <- extract_network(enrichment, enrichment_networks)
  return(clusters)
}
)
names(metabolite_BP_nets_xtr) <- names(metabolite_BP)
metabolite_BP_clusters <- lapply(metabolite_BP_nets_xtr, function(net) {
  net[["nodes"]] |>
    dplyr::select(Cluster, "Cluster size", ID) |>
    distinct() |>
    arrange(Cluster) |>
    as.data.frame()
})

metabolite_BP_clusters <- metabolite_BP_clusters |> dplyr::bind_rows()


## with only pathways
metabolite_BP |> 
  dplyr::bind_rows(.id = "metabolite") |> 
  dplyr::left_join(metabolite_BP_clusters, by = c("Description" = "ID")) |> 
  ggplot(aes(metabolite, Description, fill = RichFactor)) +
    geom_tile() +
    facet_grid(rows = vars(Cluster), scales = "free", space = "free")

## HEATMAP -----
pol.h <- met.res |>
  dplyr::filter(Metabolite.name %in% unlist(polyols)) |>
  dplyr::filter(pvalue <= 0.05) |>
  dplyr::mutate(contrast = ifelse(grepl("HC", Comparison), "Clinical_Grouping", "Biomarkers")) |>
  dplyr::left_join(lookup, by = c("Metabolite.name" = "metabolite")) |>
  dplyr::select(-pvalue) |>
  dplyr::rename(beta = Mean.beta, metabolite = Metabolite.name, comparison = Comparison, tissue = Tissue) |>
  dplyr::mutate(
    tissue = factor(tissue, levels = c("CSF", "Plasma")),
    metabolite = factor(metabolite),
    comparison = factor(comparison, levels = c("HCvsMCI", "HCvsAD", "HCvsVaD", "betaAmyloid", "pTau", "tTau"))
  ) |>
  arrange(metabolite)

met_order <- unlist(polyols)

pol.h <- pol.h |>
  tidyr::complete(metabolite = met_order, comparison, tissue)
pol.h$metabolite <- factor(pol.h$metabolite, levels = met_order)
pol.h <- make_half_triangles(pol.h, "CSF", "Plasma")
pol.h <- pol.h |>
  dplyr::left_join(lookup, by = "metabolite") |>
  dplyr::mutate(contrast = ifelse(grepl("HC", comparison), "Clinical\nGrouping", "Biomarkers")) |>
  dplyr::mutate(contrast = factor(contrast, levels = c("Clinical\nGrouping", "Biomarkers"))) |>
  dplyr::mutate(group = gsub("_", "\n", group))


pol.h$comparison <- gsub("HC", "NCI", pol.h$comparison)
comparisons_order <- c("NCIvsMCI", "NCIvsAD", "NCIvsVaD", "betaAmyloid", "pTau", "tTau")

p1 <- ggplot() +
  ## Plasma layer
  geom_polygon(
    data = pol.h |> filter(tissue == "Plasma"),
    aes(x, y,
      group = interaction(metabolite, comparison, tissue),
      fill = beta
    ), color = "black", na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = "darkblue", high = "darkorange",
    mid = "gray75", midpoint = 0, na.value = NA,
    name = "Plasma β"
  ) +
  new_scale_fill() +

  ## CSF layer
  geom_polygon(
    data = pol.h |> filter(tissue == "CSF"),
    aes(x, y,
      group = interaction(metabolite, comparison, tissue),
      fill = beta
    ), color = "black", na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = "darkolivegreen3", high = "darkred",
    mid = "gray75", midpoint = 0, na.value = NA,
    name = "CSF β"
  ) +

  ## axes
  scale_x_continuous(
    breaks = seq_along(comparisons_order),
    labels = comparisons_order
  ) +
  scale_y_continuous(
    breaks = seq_along(met_order),
    labels = met_order
  ) +
  # coord_fixed() +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(
    rows = vars(group),
    cols = vars(contrast),
    scales = "free",
    space = "free",
    labeller = labeller(group = label_wrap_gen(width = 10)) # wrap left strips
  )

# ggsave(
#   filename = file.path(dir_res, "POLYOLS_ClinicalGroups_Biomarkers_BothTissues.pdf"),
#   plot = p1, width = 12, height = 10, dpi = 300
# )
# png(file.path(dir_res, "POLYOLS_ClinicalGroups_Biomarkers_BothTissues.png"), width = 30, height = 13, units = "cm", res = 300)
# p1
# dev.off()


###### whole set of metabolites HEATMAP ===== 
mets.all <- list(
  Sugars_related = c(
    "2,3-Dihydroxypropyl phosphate",
    "Ribonic acid"
  ),
  Gut_derived = c(
    "2-Hydroxypyridine",
    "Hippuric acid / p-Cresol",
    "Benzoic acid"
  ),
  Nitrogen_waste_and_purines = c(
    "Uric acid",
    "Creatinine / Urea",
    "Creatinine"
  ),
  Lipids_related = c(
    "Pantothenic acid",
    "1-Monopalmitin",
    "Cholesterol"
  ),
  BCAA_and_organic.keto.hydroxy_acids = c(
    "Glycolic acid",
    "Glyceric acid / Ribonic acid",
    "Glyceric acid / Quinic acid",
    "Glyceric acid",
    "2,4-Dihydroxybutanoic acid",
    "2,3-Dihydroxybutanoic acid",
    "2-Hydroxybutyric acid",
    "3-Hydroxyisovaleric acid",
    "2-Keto-3-methylvaleric acid / 2-Ketoisocaproic acid",
    "2-Keto-3-methylvaleric acid"
  ),
  Amino_acids = c(
    "Valine",
    "Leucine",
    "Isoleucine",
    "Lysine",
    "Glycine",
    "Glutamic acid",
    "L-Tryptophan",
    "5-Aminopentanoic acid",
    "R-(-)-1-Amino-2-propanol",
    "Aminoethanolamine / Ethanolamine",
    "L-5-Oxoproline"
  ),
  TCA_cycle = c(
    "Citric acid",
    "Oxalic acid",
    "L-(-)-Tartaric acid",
    "L-(+)-Tartaric acid",
    "L-Threonic acid"
  ),
  Central_energy_metabolism_and_Carbohydrates = c(
    "d-Galactose",
    "d-Glucose / Erythritol",
    "D-(-)-Ribofuranose",
    "D-(-)-Xylose",
    "Cellobiose",
    "Sorbitol",
    "meso-Erythritol",
    "L-(-)-Arabitol",
    "D-Threitol",
    "3-Deoxy-erythro-pentitol"
  )
)
lookup2 <- purrr::imap_dfr(mets.all, ~ tibble(
  metabolite = .x,
  group = .y
))

met_order2 <- unlist(mets.all)

all.h <- met.res |>
  mutate(Metabolite.name = factor(Metabolite.name, levels = met_order2))

all.h <- all.h |>
  dplyr::filter(pvalue <= 0.05) |>
  dplyr::mutate(contrast = ifelse(grepl("HC", Comparison), "Clinical_Grouping", "Biomarkers")) |>
  dplyr::left_join(lookup2, by = c("Metabolite.name" = "metabolite")) |>
  dplyr::select(-pvalue) |>
  dplyr::rename(beta = Mean.beta, metabolite = Metabolite.name, comparison = Comparison, tissue = Tissue) |>
  dplyr::mutate(
    tissue = factor(tissue, levels = c("CSF", "Plasma")),
    metabolite = factor(metabolite),
    comparison = factor(comparison, levels = c("HCvsMCI", "HCvsAD", "HCvsVaD", "betaAmyloid", "pTau", "tTau"))
  ) |>
  arrange(metabolite)


all.h <- all.h |>
  tidyr::complete(metabolite = met_order2, comparison, tissue)
all.h$metabolite <- factor(all.h$metabolite, levels = met_order2)
all.h <- make_half_triangles(all.h, "CSF", "Plasma")
all.h <- all.h |>
  dplyr::left_join(lookup2, by = "metabolite") |>
  dplyr::mutate(contrast = ifelse(grepl("HC", comparison), "Clinical\nGrouping", "Biomarkers")) |>
  dplyr::mutate(contrast = factor(contrast, levels = c("Clinical\nGrouping", "Biomarkers"))) |>
  dplyr::mutate(group = gsub("_", "\n", group)) 

comparisons_order <- c("NCIvsMCI", "NCIvsAD", "NCIvsVaD", "betaAmyloid", "pTau", "tTau")

p2 <- ggplot() +
  ## Plasma layer
  geom_polygon(
    data = all.h |> filter(tissue == "Plasma", !is.na(group)),
    aes(x, y,
        group = interaction(metabolite, comparison, tissue),
        fill = beta
    ), color = "black", na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = "darkblue", high = "darkorange",
    mid = "gray75", midpoint = 0, na.value = NA,
    name = "Plasma β"
  ) +
  new_scale_fill() +
  
  ## CSF layer
  geom_polygon(
    data = all.h |> filter(tissue == "CSF", !is.na(group)),
    aes(x, y,
        group = interaction(metabolite, comparison, tissue),
        fill = beta
    ), color = "black", na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = "darkolivegreen3", high = "darkred",
    mid = "gray75", midpoint = 0, na.value = NA,
    name = "CSF β"
  ) +
  
  ## axes
  scale_x_continuous(
    breaks = seq_along(comparisons_order),
    labels = comparisons_order
  ) +

  scale_y_continuous(
    breaks = seq_along(met_order2),
    labels = met_order2
  ) +
  coord_fixed() +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(
    rows = vars(group),
    cols = vars(contrast),
    scales = "free",
    space = "free",
    labeller = labeller(group = label_wrap_gen(width = 10)) # wrap left strips
  )


# ggsave(
#   filename = file.path(dir_res, "ALLMETS_ClinicalGroups_Biomarkers_BothTissues.pdf"),
#   plot = p2, width = 12, height = 20, dpi = 300
# )
# png(file.path(dir_res, "ALLMETS_ClinicalGroups_Biomarkers_BothTissues.png"), width = 30, height = 36, units = "cm", res = 300)
# p2
# dev.off()



## metabolites abundance plots ----
library(ggpubr)
crmn_csf <- read.table(file.path(dir_data, "CSF/CSF_CRMNnormalized_ratios_UNtargeted.txt"), header = T)
crmn_csf_red <- crmn_csf |>
  dplyr::select(DDBBno, Gender, AgeAtVisit, CognitiveSyndrome, fullClass,
    betaAmyloid = Cerebrospinal.fluids_Amyloid.beta.42.pg.ml,
    pTau = Cerebrospinal.fluids_Phospho.Tau.pg.ml,
    tTau = Cerebrospinal.fluids_Total.Tau.pg.ml,
    BMI, CSF_glucose, Statins,
    d.Glucose, gluc.erythitol, d.Galactose, gluc.sorbitol,
    D.....Xylose, D.....Ribofuranose, Sorbitol,
    L.....Arabitol, D.Threitol, meso.Erythritol,
    X3.Deoxy.erythro.pentitol, X1.5.Anhydrohexitol, Ribonic.acid,
    L.Threonic.acid, Glyceric.acid) |> 
  dplyr::filter(fullClass %in% c("HC", "MCI", "AD", "VaD"))

write.table(crmn_csf_red, file.path(dir_res, "polyols_for_predictive_CSF.txt"), sep = "\t", quote = F, row.names = F)

crmn_csf_red$pTau[crmn_csf_red$pTau < 1] <- NA ## remove low outliers
crmn_csf_red$tTau[crmn_csf_red$tTau < 1] <- NA
## create target ratio variables
crmn_csf_red <- crmn_csf_red |> 
  dplyr::filter(!is.na(betaAmyloid) & !is.na(pTau) & !is.na(tTau)) |>
  dplyr::mutate(ptau.ab42 = pTau / betaAmyloid) |> 
  dplyr::mutate(ttau.ab42 = tTau / betaAmyloid)

crmn_csf_red$dementia <- ifelse(crmn_csf_red$fullClass == "HC", 0, 1)
crmn_csf_red <- crmn_csf_red |> 
  dplyr::rename(ABETA42 = betaAmyloid, PTAU= pTau, TAU = tTau)

write.table(crmn_csf_red, file.path(dir_data, "csf_ratios2.csv"), sep = ";", quote = F, row.names = F)

## cholesterol check
csf_targ <- read.table(file.path(dir_data, "CSF/CSF_targeted_logtransformed.txt"), sep = "\t", header = T)
csf_targ <- csf_targ |> 
  dplyr::rename(AB = Cerebrospinal.fluids_Amyloid.beta.42.pg.ml,
                pTAU = Cerebrospinal.fluids_Phospho.Tau.pg.ml,
                tTAU = Cerebrospinal.fluids_Total.Tau.pg.ml) |> 
  dplyr::mutate(across(.cols = c("AB", "pTAU", "tTAU"), log)) |> 
  dplyr::mutate(across(.cols = c("AB", "pTAU", "tTAU"), ~ifelse(.x < 0, NA_real_, .x))) |> 
  as.data.frame()


csf_targ <- split(csf_targ, csf_targ$Statins)
names(csf_targ) <- c("No_statin", "Statin")

vars <- c("AB", "pTAU", "tTAU")
statins_cholesterol_ev <- lapply(csf_targ, function(df){
  output_dataframe <- data.frame()
  for(var in vars){
    tmp_results <- glm(df[[var]] ~ 1 + Cholesterol + AgeAtVisit + Gender, data = df, family = "gaussian")
    tmp_results_s <- summary(tmp_results)
    null_deviance <- tmp_results$null.deviance
    residual_deviance <- tmp_results$deviance
    # Calculate R-squared
    r_squared <- 1 - (residual_deviance / null_deviance)
    CIs <- confint(tmp_results)
    
    output_dataframe[var, "r2"] <- r_squared
    output_dataframe[var, "estimate"] <- tmp_results$coefficients[2]
    output_dataframe[var, "se"] <- tmp_results_s$coefficients[2, 2]
    output_dataframe[var, "CI_upper"] <- CIs[2, 2]
    output_dataframe[var, "CI_lower"] <- CIs[2, 1]
    output_dataframe[var, "p_value"] <- tmp_results_s$coefficients[2, 4] 
  }
  return(output_dataframe)
})

bind_rows(statins_cholesterol_ev, .id = "Group") |> 
  tibble::rownames_to_column("Biomarker") |> 
  tidyr::separate(Biomarker, sep = "\\.", into = "Biomarker") |> 
  dplyr::mutate(FDR = p.adjust(p_value, method = "fdr"))
