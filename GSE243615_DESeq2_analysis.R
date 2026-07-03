# ============================================
# RNA-seq DE Analysis — GSE243615
# Author: Charlie Stewart-Lewis
# Date: July 2026
# Description: DESeq2 analysis of MZ1/PDD/GSK
#              combination drug treatments
# ============================================
#
# KEY REFERENCES
# --------------
# Dataset:
#   Mori Y, Akizuki Y, Honda R, Takao M, et al. & Ohtake F (2024).
#   Intrinsic signaling pathways modulate targeted protein degradation.
#   Nat Commun. https://doi.org/10.1038/s41467-024-49519-z
#   GEO accession: GSE243615
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243615
#
# ashr:
#   Stephens M, et al. (2023). ashr: Methods for Adaptive Shrinkage,
#   using Empirical Bayes
#   https://doi.org/10.32614/CRAN.package.ashr
#   https://CRAN.R-project.org/package=ashr
#
# clusterProfiler:
#   S. Xu, et al. (2024). Using clusterProfiler to characterize
#   multiomics data
#   Nature Protocols, 19(11):3292-3320
#   https://doi.org/10.1038/s41596-024-01020-z
#
# DESeq2:
#   Love MI, Huber W, Anders S (2014). Moderated estimation of
#   fold change and dispersion for RNA-seq data with DESeq2.
#   Genome Biology, 15(12), 550.
#   https://doi.org/10.1186/s13059-014-0550-8
#
# enrichplot:
#   Yu G (2026). enrichplot: Visualization of Functional Enrichment Result.
#   Bioconductor. https://doi.org/10.18129/B9.bioc.enrichplot
#
# GEOquery:
#   Davis S, Meltzer P (2007). GEOquery: a bridge between the
#   Gene Expression Omnibus (GEO) and BioConductor.
#   Bioinformatics, 14, 1846-1847.
#   https://doi.org/10.1093/bioinformatics/btm254
#
# knitr:
#   Xie Y (2025).knitr: A General-Purpose Package for
#   Dynamic Report Generation in R.
#   https://yihui.org/knitr
#
# org.Hs.eg.db:
#   Carlson M (2026). org.Hs.eg.db: Genome wide annotation for Human.
#   Bioconductor. https://doi.org/10.18129/B9.bioc.org.Hs.eg.db
#
# pheatmap:
#   Kolde R (2025). pheatmap: Pretty Heatmaps.
#   https://doi.org/10.32614/CRAN.package.pheatmap
#   https://CRAN.R-project.org/package=pheatmap
#
# R:
#   R Core Team (2024). R: A language and environment for
#   statistical computing. R Foundation for Statistical Computing.
#   https://www.R-project.org/
# ============================================
#
# === 1. LOAD LIBRARIES ===
# Run install lines once only, then comment out:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "tidyverse", "pheatmap",
#                        "clusterProfiler", "org.Hs.eg.db", "GEOquery",
#                        "ashr", "enrichplot"))
# install.packages("knitr")

library(clusterProfiler)
library(DESeq2)
library(GEOquery)
library(knitr)
library(org.Hs.eg.db)
library(pheatmap)
library(tidyverse)

# === 2. LOAD DATA ===
gse <- getGEO("GSE243615", GSEMatrix = TRUE)
getGEOSuppFiles("GSE243615", makeDirectory = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# Note: raw counts file must be downloaded manually from GEO
# URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243615
counts <- read.table("GSE243615/GSE243615_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                     header       = TRUE,
                     row.names    = 1,
                     sep          = "\t",
                     check.names  = FALSE,
                     comment.char = "")

# === 3. BUILD METADATA ===
## replace omitted geo_accession: "GSM7791896"
col_data <- pData(gse[[1]])
col_data <- col_data[col_data$geo_accession != "GSM_MISSING", ]

col_data[, c("geo_accession", "title", "characteristics_ch1.3")]
colnames(counts)

col_data$condition <- gsub("treatment: ", "", col_data$characteristics_ch1.3)
col_data$condition <- gsub(" and ", "_", col_data$condition)

col_data_ordered <- col_data[colnames(counts), ]

metadata <- data.frame(
  row.names = colnames(counts),
  condition = col_data_ordered$condition
)

metadata$condition <- factor(metadata$condition,
                             levels = c("vehicle", "PDD", "MZ1",
                                        "GSK", "MZ1_PDD", "MZ1_GSK"))

# === 4. BUILD DESeqDataSet & PRE-FILTER  ===
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# === 5. RUN DESeq2 ===
dds <- DESeq(dds)

# === 6. QC ===
## PCA (Principle Component Analysis) plot
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA of all samples")

png("outputs/pca_plot.png", width = 8, height = 6, units = "in", res = 300)
print(plotPCA(vsd, intgroup = "condition") +
        theme_minimal() +
        ggtitle("PCA of all samples"))
dev.off()

## sample distance heat map
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample distances")

png("outputs/sample_distances.png", width = 8, height = 7, units = "in", res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample distances")
dev.off()

# === 7. DE TABLE & LFC SHRINKAGE  ===
comparisons <- list(
  "PDD vs vehicle"     = c("PDD",     "vehicle"),
  "MZ1 vs vehicle"     = c("MZ1",     "vehicle"),
  "GSK vs vehicle"     = c("GSK",     "vehicle"),
  "MZ1+PDD vs vehicle" = c("MZ1_PDD", "vehicle"),
  "MZ1+GSK vs vehicle" = c("MZ1_GSK", "vehicle"),
  "MZ1+PDD vs MZ1"     = c("MZ1_PDD", "MZ1"),
  "MZ1+PDD vs PDD"     = c("MZ1_PDD", "PDD"),
  "MZ1+GSK vs MZ1"     = c("MZ1_GSK", "MZ1"),
  "MZ1+GSK vs GSK"     = c("MZ1_GSK", "GSK")
)

results_list       <- list()
results_list_shrunk <- list()

for (name in names(comparisons)) {
  contrast <- c("condition", comparisons[[name]])
  results_list[[name]]        <- results(dds, contrast = contrast, alpha = 0.05)
  results_list_shrunk[[name]] <- lfcShrink(dds, contrast = contrast, type = "ashr")
}

de_summary <- data.frame(
  comparison  = names(results_list),
  total_genes = sapply(results_list, function(r) sum(!is.na(r$padj))),
  up          = sapply(results_list, function(r) sum(r$padj < 0.05 & r$log2FoldChange >  1, na.rm = TRUE)),
  down        = sapply(results_list, function(r) sum(r$padj < 0.05 & r$log2FoldChange < -1, na.rm = TRUE)),
  total_sig   = sapply(results_list, function(r) sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE))
)

kable(de_summary, format = "markdown")
write.csv(de_summary, "outputs/de_summary.csv", row.names = FALSE)

# === 8. VISUALISATION ===
## volcano plots
make_volcano <- function(res, title) {
  df <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1,
           direction   = case_when(
             significant & log2FoldChange > 0 ~ "up",
             significant & log2FoldChange < 0 ~ "down",
             TRUE ~ "ns"
           ))
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), colour = direction)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_colour_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1),     linetype = "dashed") +
    theme_minimal() +
    ggtitle(title)
}
make_volcano(results_list_shrunk[["PDD vs vehicle"]],     "PDD vs Vehicle")
make_volcano(results_list_shrunk[["MZ1 vs vehicle"]],     "MZ1 vs Vehicle")
make_volcano(results_list_shrunk[["GSK vs vehicle"]],     "GSK vs Vehicle")
make_volcano(results_list_shrunk[["MZ1+PDD vs vehicle"]], "MZ1+PDD vs Vehicle")
make_volcano(results_list_shrunk[["MZ1+GSK vs vehicle"]], "MZ1+GSK vs Vehicle")
make_volcano(results_list_shrunk[["MZ1+PDD vs MZ1"]],    "MZ1+PDD vs MZ1")
make_volcano(results_list_shrunk[["MZ1+PDD vs PDD"]],    "MZ1+PDD vs PDD")
make_volcano(results_list_shrunk[["MZ1+GSK vs MZ1"]],    "MZ1+GSK vs MZ1")
make_volcano(results_list_shrunk[["MZ1+GSK vs GSK"]],    "MZ1+GSK vs GSK")

dir.create("outputs", showWarnings = FALSE)

pdf("outputs/volcano_plots.pdf", width = 8, height = 6)
for (name in names(results_list_shrunk)) {
  print(make_volcano(results_list_shrunk[[name]], name))
}
dev.off()

## heatmap of top DE genes
annotation_col <- data.frame(
  condition = metadata$condition,
  row.names = colnames(counts)
)

ann_colours <- list(
  condition = c(
    vehicle  = "grey80",
    PDD      = "steelblue",
    MZ1      = "darkorange",
    GSK      = "forestgreen",
    MZ1_PDD  = "purple",
    MZ1_GSK  = "firebrick"
  )
)
make_heatmap <- function(comparison_name, n_genes = 50) {
  top_genes <- results_list[[comparison_name]] %>%
    as.data.frame() %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj) %>%
    head(n_genes) %>%
    rownames()
  
  if (length(top_genes) < 2) {
    message("Too few significant genes for: ", comparison_name)
    return(NULL)
  }
  
  mat <- assay(vsd)[top_genes, ]
  mat_scaled <- t(scale(t(mat)))
  
  pheatmap(mat_scaled,
           annotation_col    = annotation_col,
           annotation_colors = ann_colours,
           show_rownames     = TRUE,
           show_colnames     = FALSE,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           scale             = "none",
           color             = colorRampPalette(c("navy", "white", "firebrick"))(100),
           main              = paste("Top", n_genes, "DE genes:", comparison_name),
           fontsize_row      = 6)
}
make_heatmap("PDD vs vehicle")
make_heatmap("MZ1 vs vehicle")
make_heatmap("GSK vs vehicle")
make_heatmap("MZ1+PDD vs vehicle")
make_heatmap("MZ1+GSK vs vehicle")
make_heatmap("MZ1+PDD vs MZ1")
make_heatmap("MZ1+PDD vs PDD")
make_heatmap("MZ1+GSK vs MZ1")
make_heatmap("MZ1+GSK vs GSK")

pdf("outputs/heatmaps.pdf", width = 10, height = 12)
for (name in names(comparisons)) {
  tryCatch(
    make_heatmap(name),
    error = function(e) message("Error in heatmap for ", name, ": ", e$message)
  )
}
dev.off()

# === 9. PATHWAY ENRICHMENT ===
## collecting significant genes
run_enrichment <- function(res_obj, title) {
  sig <- res_obj %>%
    as.data.frame() %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    rownames()
  
  if (length(sig) < 10) {
    message("Too few significant genes for ", title, " — skipping")
    return(NULL)
  }
## GO enrichment
  go <- enrichGO(gene          = sig,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = "ENTREZID",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)
  
  print(dotplot(go, showCategory = 20))
  return(go)
}

go_results <- list()
go_results[["PDD vs vehicle"]]     <- run_enrichment(results_list[["PDD vs vehicle"]],     "PDD vs Vehicle")
go_results[["MZ1 vs vehicle"]]     <- run_enrichment(results_list[["MZ1 vs vehicle"]],     "MZ1 vs Vehicle")
go_results[["GSK vs vehicle"]]     <- run_enrichment(results_list[["GSK vs vehicle"]],     "GSK vs Vehicle")
go_results[["MZ1+PDD vs vehicle"]] <- run_enrichment(results_list[["MZ1+PDD vs vehicle"]], "MZ1+PDD vs Vehicle")
go_results[["MZ1+GSK vs vehicle"]] <- run_enrichment(results_list[["MZ1+GSK vs vehicle"]], "MZ1+GSK vs Vehicle")
go_results[["MZ1+PDD vs MZ1"]]    <- run_enrichment(results_list[["MZ1+PDD vs MZ1"]],    "MZ1+PDD vs MZ1")
go_results[["MZ1+PDD vs PDD"]]    <- run_enrichment(results_list[["MZ1+PDD vs PDD"]],    "MZ1+PDD vs PDD")
go_results[["MZ1+GSK vs MZ1"]]    <- run_enrichment(results_list[["MZ1+GSK vs MZ1"]],    "MZ1+GSK vs MZ1")
go_results[["MZ1+GSK vs GSK"]]    <- run_enrichment(results_list[["MZ1+GSK vs GSK"]],    "MZ1+GSK vs GSK")

pdf("outputs/go_enrichment_dotplots.pdf", width = 10, height = 8)
for (name in names(go_results)) {
  if (!is.null(go_results[[name]])) {
    print(dotplot(go_results[[name]], 
                  showCategory = 20, 
                  title        = name))
  }
}
dev.off()

# === 10. SESSION INFO ===
sessionInfo()
