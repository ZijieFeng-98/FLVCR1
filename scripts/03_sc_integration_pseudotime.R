# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 03: Cross-dataset Integration & Pseudotime Analysis
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

cat("=== FLVCR Pipeline: Integration & Pseudotime ===\n\n")

# ── Configuration ──
dir_in <- "data_intermediate"
dir_out <- "data_results"
dir_figs <- "figures"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# ── Load annotated objects ──
cat("Loading annotated scRNA-seq objects...\n")

load_if_exists <- function(file_path) {
  if (file.exists(file_path)) {
    cat("  ✓ Loading:", basename(file_path), "\n")
    return(readRDS(file_path))
  } else {
    cat("  ✗ Not found:", basename(file_path), "\n")
    return(NULL)
  }
}

obj_120508 <- load_if_exists(file.path(dir_in, "GSE120508_qc_annotated.rds"))
obj_109037 <- load_if_exists(file.path(dir_in, "GSE109037_qc_annotated.rds"))
obj_107644 <- load_if_exists(file.path(dir_in, "GSE107644_qc_annotated.rds"))

# Filter to available objects
objs_human <- Filter(Negate(is.null), list(obj_120508, obj_109037))

if (length(objs_human) == 0) {
  stop("No human datasets found. Please run 02_scRNA_QC_annotate.R first.")
}

# ── Integrate human datasets ──
cat("\n=== Integrating Human Datasets ===\n")

if (length(objs_human) > 1) {
  cat("Integrating", length(objs_human), "human datasets...\n")
  
  # Select integration features
  features <- SelectIntegrationFeatures(objs_human, nfeatures = 3000)
  
  # Prepare for integration
  objs_human <- lapply(objs_human, function(x) {
    DefaultAssay(x) <- "SCT"
    return(x)
  })
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(
    object.list = objs_human,
    normalization.method = "SCT",
    anchor.features = features,
    verbose = TRUE
  )
  
  # Integrate
  int_human <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    verbose = TRUE
  )
  
  # Standard workflow on integrated data
  DefaultAssay(int_human) <- "integrated"
  int_human <- RunPCA(int_human, npcs = 50, verbose = FALSE)
  int_human <- RunUMAP(int_human, dims = 1:30, verbose = FALSE)
  int_human <- FindNeighbors(int_human, dims = 1:30, verbose = FALSE)
  int_human <- FindClusters(int_human, resolution = 0.6, verbose = FALSE)
  
} else {
  cat("Single human dataset - skipping integration\n")
  int_human <- objs_human[[1]]
  DefaultAssay(int_human) <- "SCT"
}

# Visualize integrated data
p1 <- DimPlot(int_human, reduction = "umap", group.by = "celltype", label = TRUE) +
  ggtitle("Integrated Human Testis - Cell Types")

p2 <- DimPlot(int_human, reduction = "umap", group.by = "dataset") +
  ggtitle("Batch Effect Check")

p_combined <- p1 + p2
ggsave(file.path(dir_figs, "UMAP_integrated_human.png"), p_combined, width = 14, height = 6)

# ── Pseudotime analysis using Slingshot ──
cat("\n=== Pseudotime Analysis ===\n")

# Filter to germ cells for trajectory
germ_cells <- c("Spermatogonia", "Spermatocytes", "RoundSpermatids", "ElongatingSpermatids")
int_germ <- subset(int_human, subset = celltype %in% germ_cells)

cat("Germ cells for trajectory:", ncol(int_germ), "\n")

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(int_germ, assay = "SCT")

# Run Slingshot
cat("Running Slingshot trajectory inference...\n")
sce <- slingshot(
  sce,
  clusterLabels = int_germ$celltype,
  reducedDim = "UMAP",
  start.clus = "Spermatogonia",
  end.clus = "ElongatingSpermatids"
)

# Extract pseudotime
pseudotime <- slingPseudotime(sce)[, 1]  # First lineage
int_germ$pseudotime <- pseudotime

# Visualize pseudotime
p_pseudo <- ggplot(data.frame(UMAP_1 = reducedDim(sce, "UMAP")[, 1],
                               UMAP_2 = reducedDim(sce, "UMAP")[, 2],
                               Pseudotime = pseudotime,
                               CellType = int_germ$celltype),
                   aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Spermatogenesis Pseudotime")

ggsave(file.path(dir_figs, "UMAP_pseudotime_human.png"), p_pseudo, width = 8, height = 6)

# ── FLVCR1/2 expression along pseudotime ──
cat("\n=== FLVCR1/2 Expression Dynamics ===\n")

genes_of_interest <- c("FLVCR1", "FLVCR2")

# Check if genes are present
genes_present <- intersect(genes_of_interest, rownames(int_germ))
if (length(genes_present) == 0) {
  warning("FLVCR1/2 not found in data. Check gene names.")
} else {
  cat("Analyzing genes:", paste(genes_present, collapse = ", "), "\n")
  
  # Extract expression
  DefaultAssay(int_germ) <- "SCT"
  expr_data <- FetchData(int_germ, vars = c("pseudotime", "celltype", genes_present))
  
  # Plot expression vs pseudotime
  library(tidyr)
  expr_long <- expr_data %>%
    pivot_longer(cols = all_of(genes_present), names_to = "Gene", values_to = "Expression")
  
  p_expr <- ggplot(expr_long, aes(x = pseudotime, y = Expression, color = Gene)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 10), se = TRUE, linewidth = 1.5) +
    scale_color_manual(values = c("FLVCR1" = "#E64B35", "FLVCR2" = "#4DBBD5")) +
    theme_minimal(base_size = 12) +
    labs(x = "Pseudotime", y = "Expression (SCT)", 
         title = "FLVCR1/2 Expression Along Spermatogenesis") +
    theme(legend.position = "top")
  
  ggsave(file.path(dir_figs, "FLVCR1_FLVCR2_pseudotime_human.png"), p_expr, 
         width = 8, height = 6, dpi = 300)
  
  # Stage-wise boxplots
  expr_data$celltype <- factor(expr_data$celltype, 
                                levels = c("Spermatogonia", "Spermatocytes", 
                                           "RoundSpermatids", "ElongatingSpermatids"))
  
  p_stage <- ggplot(expr_long, aes(x = celltype, y = Expression, fill = Gene)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("FLVCR1" = "#E64B35", "FLVCR2" = "#4DBBD5")) +
    theme_minimal(base_size = 12) +
    labs(x = "Cell Type", y = "Expression (SCT)", 
         title = "FLVCR1/2 Expression by Stage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  
  ggsave(file.path(dir_figs, "FLVCR1_FLVCR2_by_stage_human.png"), p_stage, 
         width = 8, height = 6, dpi = 300)
  
  # Statistical test: tradeSeq for differential expression along pseudotime
  cat("\nFitting tradeSeq GAM models...\n")
  
  tryCatch({
    # Prepare count matrix for tradeSeq
    counts_matrix <- GetAssayData(int_germ, assay = "RNA", slot = "counts")[genes_present, ]
    
    # Fit GAM
    gam_fit <- fitGAM(
      counts = as.matrix(counts_matrix),
      pseudotime = pseudotime,
      cellWeights = matrix(1, ncol = 1, nrow = ncol(counts_matrix)),
      nknots = 6,
      verbose = TRUE
    )
    
    # Test for association with pseudotime
    assoc_test <- associationTest(gam_fit)
    print(assoc_test)
    
    # Save results
    write.csv(assoc_test, file.path(dir_out, "FLVCR_pseudotime_association_test.csv"))
    
  }, error = function(e) {
    cat("tradeSeq fitting failed:", e$message, "\n")
    cat("Proceeding with descriptive statistics.\n")
  })
  
  # Calculate fold-changes by stage
  stage_summary <- expr_long %>%
    group_by(celltype, Gene) %>%
    summarize(
      Mean = mean(Expression),
      Median = median(Expression),
      SD = sd(Expression),
      .groups = "drop"
    ) %>%
    arrange(Gene, celltype)
  
  write.csv(stage_summary, file.path(dir_out, "FLVCR_stage_summary_human.csv"), row.names = FALSE)
  cat("\n✓ Stage summary saved\n")
  print(stage_summary)
}

# ── Save integrated object ──
saveRDS(int_human, file.path(dir_in, "integrated_human.rds"))
saveRDS(int_germ, file.path(dir_in, "integrated_human_germ_pseudotime.rds"))

cat("\n=== Integration & Pseudotime Complete ===\n")
cat("Integrated objects saved to:", dir_in, "\n")
cat("Results saved to:", dir_out, "\n")
cat("Figures saved to:", dir_figs, "\n")
cat("\nNext step: Run 04_stage_profiles_FLVCR.R\n")

