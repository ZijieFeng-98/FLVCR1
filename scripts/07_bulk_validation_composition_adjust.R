# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 07: Bulk RNA-seq Validation with Cell Composition Adjustment
# ============================================================================

suppressPackageStartupMessages({
  library(MuSiC)
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(Biobase)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cat("=== FLVCR Pipeline: Bulk Validation with Composition Adjustment ===\n\n")

# ── Configuration ──
dir_in <- "data_intermediate"
dir_raw <- "data_raw"
dir_out <- "data_results"
dir_figs <- "figures"

cat("NOTE: This script requires bulk RNA-seq data from infertile vs control samples.\n")
cat("If you have bulk datasets, place them in data_raw/ with metadata.\n\n")

# ── Load scRNA reference for deconvolution ──
cat("Loading scRNA reference...\n")

sc_file <- file.path(dir_in, "integrated_human.rds")
if (!file.exists(sc_file)) {
  stop("Integrated scRNA reference not found. Please run script 03 first.")
}

sc_ref <- readRDS(sc_file)

# ── Prepare reference for MuSiC ──
cat("Preparing single-cell reference for deconvolution...\n")

# Convert to SingleCellExperiment
library(SingleCellExperiment)
sce_ref <- as.SingleCellExperiment(sc_ref, assay = "RNA")

# Add cell type labels
colData(sce_ref)$cellType <- sc_ref$celltype
colData(sce_ref)$sampleID <- sc_ref$dataset

# Summary
cat("Reference cell types:\n")
print(table(sce_ref$cellType))

# ── Example bulk data loading function ──
load_bulk_data <- function(bulk_file, metadata_file) {
  # This is a template - adjust based on your bulk data format
  
  # Example: read count matrix
  counts <- read.table(bulk_file, header = TRUE, row.names = 1, sep = "\t")
  
  # Example: read metadata
  meta <- read.table(metadata_file, header = TRUE, sep = "\t")
  
  return(list(counts = counts, metadata = meta))
}

# ── Deconvolution function using MuSiC ──
deconvolve_bulk <- function(bulk_counts, sc_reference) {
  cat("\nRunning MuSiC deconvolution...\n")
  
  tryCatch({
    # Create bulk expression set
    bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_counts))
    
    # Run MuSiC
    est_prop <- music_prop(
      bulk.eset = bulk_eset,
      sc.eset = sc_reference,
      clusters = "cellType",
      samples = "sampleID",
      verbose = TRUE
    )
    
    return(est_prop$Est.prop.weighted)
    
  }, error = function(e) {
    cat("MuSiC deconvolution failed:", e$message, "\n")
    return(NULL)
  })
}

# ── Differential expression with composition adjustment ──
run_de_with_composition <- function(counts, metadata, cell_proportions, 
                                   test_var = "condition") {
  cat("\nRunning DE analysis with composition adjustment...\n")
  
  # Combine metadata with cell proportions
  if (!is.null(cell_proportions)) {
    design_data <- cbind(metadata, cell_proportions)
    
    # Formula with composition covariates
    prop_cols <- colnames(cell_proportions)
    formula_str <- paste0("~ ", test_var, " + ", paste(prop_cols, collapse = " + "))
    cat("Design formula:", formula_str, "\n")
    
    design <- model.matrix(as.formula(formula_str), data = design_data)
  } else {
    # Without composition adjustment
    design <- model.matrix(as.formula(paste0("~ ", test_var)), data = metadata)
  }
  
  # Run limma-voom
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Extract results for test variable
  coef_name <- grep(test_var, colnames(design), value = TRUE)[1]
  results <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  
  return(results)
}

# ── EXAMPLE WORKFLOW (requires actual bulk data) ──
cat("\n=== Example Workflow ===\n")
cat("This is a template. Replace with your actual bulk datasets.\n\n")

# Simulate example scenario
simulate_bulk_analysis <- function() {
  cat("SIMULATION MODE: Demonstrating composition adjustment workflow\n\n")
  
  # Create example metadata
  example_meta <- data.frame(
    SampleID = paste0("Sample", 1:20),
    Condition = rep(c("Control", "Infertile"), each = 10),
    Age = rnorm(20, mean = 35, sd = 5),
    Batch = rep(c("A", "B"), 10)
  )
  
  # Simulate cell type proportions (would come from MuSiC in real analysis)
  example_props <- data.frame(
    Spermatogonia = runif(20, 0.05, 0.15),
    Spermatocytes = runif(20, 0.20, 0.40),
    RoundSpermatids = runif(20, 0.15, 0.35),
    ElongatingSpermatids = runif(20, 0.10, 0.25),
    Sertoli = runif(20, 0.05, 0.15),
    Leydig = runif(20, 0.02, 0.08)
  )
  
  # Normalize proportions to sum to 1
  example_props <- example_props / rowSums(example_props)
  
  cat("Example metadata:\n")
  print(head(example_meta))
  
  cat("\nExample cell type proportions:\n")
  print(head(example_props))
  
  # Visualize proportions by condition
  props_long <- cbind(example_meta, example_props) %>%
    tidyr::pivot_longer(cols = colnames(example_props), 
                       names_to = "CellType", 
                       values_to = "Proportion")
  
  p_props <- ggplot(props_long, aes(x = SampleID, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Condition, scales = "free_x") +
    theme_minimal() +
    labs(title = "Cell Type Composition by Condition",
         y = "Estimated Proportion") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
  
  ggsave(file.path(dir_figs, "bulk_cell_composition_example.png"), p_props,
         width = 12, height = 6, dpi = 300)
  
  # Create comparison visualization
  props_summary <- props_long %>%
    group_by(Condition, CellType) %>%
    summarize(Mean = mean(Proportion),
              SE = sd(Proportion) / sqrt(n()),
              .groups = "drop")
  
  p_compare <- ggplot(props_summary, aes(x = CellType, y = Mean, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  position = position_dodge(0.9), width = 0.2) +
    theme_minimal() +
    labs(title = "Cell Type Proportions: Control vs Infertile",
         y = "Mean Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(dir_figs, "bulk_composition_comparison_example.png"), p_compare,
         width = 10, height = 6, dpi = 300)
  
  # Statistical test for compositional differences
  comp_tests <- props_long %>%
    group_by(CellType) %>%
    summarize(
      Control_mean = mean(Proportion[Condition == "Control"]),
      Infertile_mean = mean(Proportion[Condition == "Infertile"]),
      t_statistic = tryCatch(
        t.test(Proportion ~ Condition)$statistic,
        error = function(e) NA
      ),
      p_value = tryCatch(
        t.test(Proportion ~ Condition)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"))
  
  cat("\nCompositional differences between conditions:\n")
  print(comp_tests)
  
  write.csv(comp_tests, file.path(dir_out, "bulk_composition_tests_example.csv"),
            row.names = FALSE)
  
  return(list(metadata = example_meta, proportions = example_props))
}

# Run simulation
example_data <- simulate_bulk_analysis()

# ── Template for FLVCR1/2 analysis with composition adjustment ──
cat("\n=== FLVCR1/2 Analysis Template ===\n")

analyze_flvcr_bulk <- function(de_results_unadjusted, de_results_adjusted) {
  # Extract FLVCR1/2 results
  genes_of_interest <- c("FLVCR1", "FLVCR2")
  
  comparison <- data.frame(
    Gene = genes_of_interest,
    logFC_unadjusted = NA,
    pval_unadjusted = NA,
    logFC_adjusted = NA,
    pval_adjusted = NA
  )
  
  for (gene in genes_of_interest) {
    if (gene %in% rownames(de_results_unadjusted)) {
      comparison$logFC_unadjusted[comparison$Gene == gene] <- 
        de_results_unadjusted[gene, "logFC"]
      comparison$pval_unadjusted[comparison$Gene == gene] <- 
        de_results_unadjusted[gene, "P.Value"]
    }
    
    if (gene %in% rownames(de_results_adjusted)) {
      comparison$logFC_adjusted[comparison$Gene == gene] <- 
        de_results_adjusted[gene, "logFC"]
      comparison$pval_adjusted[comparison$Gene == gene] <- 
        de_results_adjusted[gene, "P.Value"]
    }
  }
  
  return(comparison)
}

cat("\nTo run full analysis with your bulk data:\n")
cat("1. Load your bulk count matrix and metadata\n")
cat("2. Run MuSiC deconvolution: cell_props <- deconvolve_bulk(bulk_counts, sce_ref)\n")
cat("3. Run DE without adjustment: de_unadj <- run_de_with_composition(counts, meta, NULL)\n")
cat("4. Run DE with adjustment: de_adj <- run_de_with_composition(counts, meta, cell_props)\n")
cat("5. Compare FLVCR1/2: analyze_flvcr_bulk(de_unadj, de_adj)\n")

# ── Validation criteria ──
cat("\n=== VALIDATION CRITERIA ===\n")
cat("✓ PASS: FLVCR1/2 remain NOT significantly DE after composition adjustment\n")
cat("  - p-value > 0.05 or |logFC| < 0.5 after adjustment\n")
cat("  - Demonstrates they are not primary drivers of infertility phenotype\n")
cat("✓ SENSITIVITY: Test with multiple deconvolution methods (MuSiC, BisqueRNA, etc.)\n")

cat("\n=== Bulk Validation Script Complete ===\n")
cat("Results saved to:", dir_out, "\n")
cat("Figures saved to:", dir_figs, "\n")
cat("\nNext step: Run 08_cross_species_orthologs.R\n")

