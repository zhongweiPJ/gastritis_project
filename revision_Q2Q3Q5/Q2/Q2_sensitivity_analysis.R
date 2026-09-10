# ============================================================
# Q2 Sensitivity Analysis
# Nominal-P-value DEG list vs FDR-controlled DEG list
#
# Primary goal:
#   Compare pathway-level ORA results obtained from:
#   1) nominal DEG: P < 0.05 & |log2FC| > 0.5
#   2) FDR DEG: BH-adjusted P < 0.05 & |log2FC| > 0.5
#
# The analysis is performed separately for genes higher in:
#   - subtype B
#   - subtype A
#
# Ranked-list GSEA should remain the PRIMARY pathway analysis.
# ORA here is a sensitivity analysis / complementary analysis.
# ============================================================

# ---------- 0. Packages ----------
# If needed:
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
# install.packages(c("dplyr", "readr", "ggplot2", "pheatmap"))
# install.packages("msigdbr")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
})

# ---------- 1. Paths ----------
# Change this to the folder containing the supplied CSV files
indir <- "Q2_sensitivity_analysis"
outdir <- file.path(indir, "ORA_results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Read gene lists ----------
nom_B <- read_csv(file.path(indir, "nominal_B_higher.csv"), show_col_types = FALSE)
nom_A <- read_csv(file.path(indir, "nominal_A_higher.csv"), show_col_types = FALSE)
fdr_B <- read_csv(file.path(indir, "FDR_B_higher.csv"), show_col_types = FALSE)
fdr_A <- read_csv(file.path(indir, "FDR_A_higher.csv"), show_col_types = FALSE)

# ---------- 3. Convert SYMBOL -> ENTREZ ----------
symbol_to_entrez <- function(genes) {
  x <- bitr(
    unique(genes),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  x <- x[!duplicated(x$SYMBOL), ]
  x
}

# ---------- 4. ORA for GO Biological Process ----------
run_GO <- function(genes, label) {
  map <- symbol_to_entrez(genes)
  ego <- enrichGO(
    gene          = map$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  res <- as.data.frame(ego)
  if (nrow(res) > 0) {
    res <- res %>% arrange(p.adjust)
    write_csv(res, file.path(outdir, paste0(label, "_GO_BP.csv")))
  }
  ego
}

# ---------- 5. ORA for KEGG ----------
run_KEGG <- function(genes, label) {
  map <- symbol_to_entrez(genes)
  ek <- enrichKEGG(
    gene          = map$ENTREZID,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  res <- as.data.frame(ek)
  if (nrow(res) > 0) {
    res <- res %>% arrange(p.adjust)
    write_csv(res, file.path(outdir, paste0(label, "_KEGG.csv")))
  }
  ek
}

# ---------- 6. ORA using MSigDB Hallmark ----------
# msigdbr API differs slightly between package versions.
# The following is compatible with current versions; if your
# installed version uses collection=, replace category/subcategory
# accordingly.
run_Hallmark <- function(genes, label) {

  msig <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ) %>%
    select(gs_name, gene_symbol)

  # Create TERM2GENE and run enricher
  term2gene <- msig %>% select(gs_name, gene_symbol)

  enr <- enricher(
    gene          = unique(genes),
    TERM2GENE     = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )

  res <- as.data.frame(enr)
  if (nrow(res) > 0) {
    res <- res %>% arrange(p.adjust)
    write_csv(res, file.path(outdir, paste0(label, "_Hallmark.csv")))
  }
  enr
}

# ---------- 7. Run all four analyses ----------
ora_objects <- list(
  nominal_B = list(genes = nom_B$gene, label = "nominal_B"),
  nominal_A = list(genes = nom_A$gene, label = "nominal_A"),
  FDR_B     = list(genes = fdr_B$gene, label = "FDR_B"),
  FDR_A     = list(genes = fdr_A$gene, label = "FDR_A")
)

for (x in ora_objects) {
  run_GO(x$genes, x$label)
  run_KEGG(x$genes, x$label)
  run_Hallmark(x$genes, x$label)
}

# ============================================================
# 8. Sensitivity comparison
# ============================================================

# Helper to read enrichment output and add metadata
read_enrich <- function(file, method, direction, database) {
  if (!file.exists(file)) return(NULL)
  x <- read_csv(file, show_col_types = FALSE)
  if (nrow(x) == 0) return(NULL)
  x %>%
    mutate(
      method = method,
      direction = direction,
      database = database
    )
}

files <- list(
  read_enrich(file.path(outdir, "nominal_B_GO_BP.csv"),
              "Nominal", "B_higher", "GO_BP"),
  read_enrich(file.path(outdir, "FDR_B_GO_BP.csv"),
              "FDR", "B_higher", "GO_BP"),
  read_enrich(file.path(outdir, "nominal_A_GO_BP.csv"),
              "Nominal", "A_higher", "GO_BP"),
  read_enrich(file.path(outdir, "FDR_A_GO_BP.csv"),
              "FDR", "A_higher", "GO_BP"),

  read_enrich(file.path(outdir, "nominal_B_KEGG.csv"),
              "Nominal", "B_higher", "KEGG"),
  read_enrich(file.path(outdir, "FDR_B_KEGG.csv"),
              "FDR", "B_higher", "KEGG"),
  read_enrich(file.path(outdir, "nominal_A_KEGG.csv"),
              "Nominal", "A_higher", "KEGG"),
  read_enrich(file.path(outdir, "FDR_A_KEGG.csv"),
              "FDR", "A_higher", "KEGG"),

  read_enrich(file.path(outdir, "nominal_B_Hallmark.csv"),
              "Nominal", "B_higher", "Hallmark"),
  read_enrich(file.path(outdir, "FDR_B_Hallmark.csv"),
              "FDR", "B_higher", "Hallmark"),
  read_enrich(file.path(outdir, "nominal_A_Hallmark.csv"),
              "Nominal", "A_higher", "Hallmark"),
  read_enrich(file.path(outdir, "FDR_A_Hallmark.csv"),
              "FDR", "A_higher", "Hallmark")
)

comparison <- bind_rows(files)

if (nrow(comparison) > 0) {
  write_csv(
    comparison,
    file.path(outdir, "sensitivity_all_enrichment_results.csv")
  )
}

# ---------- 9. Define significant pathways ----------
# For ORA, use BH-adjusted P < 0.05.
sig <- comparison %>%
  filter(!is.na(p.adjust), p.adjust < 0.05) %>%
  mutate(pathway = Description)

write_csv(
  sig,
  file.path(outdir, "sensitivity_significant_pathways.csv")
)

# ---------- 10. Identify pathway overlap ----------
# Compare pathway names within each direction/database.
overlap_summary <- sig %>%
  group_by(direction, database, pathway) %>%
  summarise(
    nominal_sig = any(method == "Nominal"),
    FDR_sig     = any(method == "FDR"),
    .groups = "drop"
  ) %>%
  mutate(
    classification = case_when(
      nominal_sig & FDR_sig ~ "Concordant",
      nominal_sig & !FDR_sig ~ "Nominal-only",
      !nominal_sig & FDR_sig ~ "FDR-only",
      TRUE ~ "Other"
    )
  )

write_csv(
  overlap_summary,
  file.path(outdir, "sensitivity_pathway_overlap.csv")
)

# ---------- 11. Summary counts ----------
summary_counts <- overlap_summary %>%
  group_by(direction, database, classification) %>%
  summarise(n = n(), .groups = "drop")

write_csv(
  summary_counts,
  file.path(outdir, "sensitivity_pathway_summary.csv")
)

# ---------- 12. Plot: concordance ----------
plot_df <- summary_counts %>%
  filter(classification != "Other")

p <- ggplot(
  plot_df,
  aes(x = classification, y = n)
) +
  geom_col() +
  facet_grid(direction ~ database, scales = "free_y") +
  labs(
    x = NULL,
    y = "Number of pathways",
    title = "Sensitivity analysis: nominal vs FDR-controlled enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "Figure_Sensitivity_pathway_concordance.pdf"),
  p,
  width = 10,
  height = 6
)

# ---------- 13. Recommended biological interpretation ----------
# Use the pathway-level overlap, not the raw number of nominal DEGs,
# to judge robustness.
#
# Recommended reporting:
#   - Concordant pathways = significant under both thresholds
#   - Nominal-only pathways = exploratory findings
#   - FDR-only pathways = should be noted but are usually uncommon
#
# Do NOT describe nominal-only pathways as statistically significant
# differential-expression-derived findings.
# ============================================================
