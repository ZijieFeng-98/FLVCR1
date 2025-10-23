# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 04: Stage-Resolved Expression Profiles
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(patchwork)
})

cat("=== FLVCR Pipeline: Stage-Resolved Expression Profiles ===\n\n")

# ── Configuration ──
dir_in <- "data_intermediate"
dir_out <- "data_results"
dir_figs <- "figures"

# Genes of interest
genes_flvcr <- c("FLVCR1", "FLVCR2")
genes_context <- c("PRM1", "PRM2", "TNP1", "TNP2", "CRISP2", "DNAJC5B", 
                   "ACR", "SPACA7", "AKAP4", "ODF2")

# ── Helper function: Calculate stage averages with bootstrap CI ──
calc_stage_avg <- function(obj, genes, n_boot = 100) {
  cat("Calculating stage averages for", unique(obj$dataset), "...\n")
  
  # Filter to genes present
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    warning("No genes found in dataset")
    return(NULL)
  }
  
  # Set identity to cell type
  Idents(obj) <- "celltype"
  
  # Average expression by stage
  avg_expr <- AverageExpression(obj, features = genes, assays = "SCT")$SCT
  
  # Bootstrap confidence intervals
  stages <- levels(Idents(obj))
  results_list <- list()
  
  for (stage in stages) {
    cells <- WhichCells(obj, idents = stage)
    if (length(cells) == 0) next
    
    expr_mat <- GetAssayData(obj, assay = "SCT", slot = "data")[genes, cells, drop = FALSE]
    
    # Bootstrap
    boot_means <- replicate(n_boot, {
      sampled_cells <- sample(cells, size = min(length(cells), 100), replace = TRUE)
      rowMeans(expr_mat[, sampled_cells, drop = FALSE])
    })
    
    # Calculate statistics
    for (i in seq_along(genes)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        Dataset = unique(obj$dataset),
        Stage = stage,
        Gene = genes[i],
        Mean = mean(boot_means[i, ]),
        SE = sd(boot_means[i, ]),
        CI_lower = quantile(boot_means[i, ], 0.025),
        CI_upper = quantile(boot_means[i, ], 0.975),
        N_cells = length(cells)
      )
    }
  }
  
  return(bind_rows(results_list))
}

# ── Load and process datasets ──
cat("\nLoading datasets...\n")

datasets <- list(
  GSE120508 = file.path(dir_in, "GSE120508_qc_annotated.rds"),
  GSE109037 = file.path(dir_in, "GSE109037_qc_annotated.rds"),
  GSE107644 = file.path(dir_in, "GSE107644_qc_annotated.rds")
)

results_all <- list()

for (name in names(datasets)) {
  if (file.exists(datasets[[name]])) {
    obj <- readRDS(datasets[[name]])
    
    # Determine genes based on species
    if (grepl("107644", name)) {  # Mouse dataset
      genes_test <- c("Flvcr1", "Flvcr2", "Prm1", "Prm2", "Tnp1", "Tnp2", 
                      "Crisp2", "Dnajc5b", "Acr", "Akap4", "Odf2")
    } else {
      genes_test <- c(genes_flvcr, genes_context)
    }
    
    result <- calc_stage_avg(obj, genes_test, n_boot = 100)
    if (!is.null(result)) {
      results_all[[name]] <- result
    }
  } else {
    cat("  ✗ Not found:", name, "\n")
  }
}

# Combine results
if (length(results_all) > 0) {
  combined_results <- bind_rows(results_all)
  write.csv(combined_results, file.path(dir_out, "FLVCR_stage_profiles_all.csv"), 
            row.names = FALSE)
  cat("\n✓ Combined results saved\n")
} else {
  stop("No datasets processed. Please run script 02 first.")
}

# ── Visualization: FLVCR1 vs FLVCR2 across stages ──
cat("\nGenerating visualizations...\n")

# Filter to FLVCR genes (human datasets)
flvcr_data <- combined_results %>%
  filter(Gene %in% c("FLVCR1", "FLVCR2")) %>%
  mutate(Stage = factor(Stage, 
                        levels = c("Spermatogonia", "Spermatocytes", 
                                   "RoundSpermatids", "ElongatingSpermatids")))

# Line plot with error bars
p1 <- ggplot(flvcr_data, aes(x = Stage, y = Mean, color = Gene, group = Gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, linewidth = 0.8) +
  facet_wrap(~Dataset, ncol = 1) +
  scale_color_manual(values = c("FLVCR1" = "#E64B35", "FLVCR2" = "#4DBBD5")) +
  theme_minimal(base_size = 12) +
  labs(x = "Spermatogenic Stage", y = "Mean Expression (SCT)",
       title = "FLVCR1/2 Expression Across Spermatogenesis",
       subtitle = "Error bars: 95% bootstrap CI") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave(file.path(dir_figs, "FLVCR1_FLVCR2_stage_profiles_line.png"), p1, 
       width = 8, height = 10, dpi = 300)

# Bar plot comparison
p2 <- ggplot(flvcr_data, aes(x = Stage, y = Mean, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                position = position_dodge(0.9), width = 0.3) +
  facet_wrap(~Dataset, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("FLVCR1" = "#E64B35", "FLVCR2" = "#4DBBD5")) +
  theme_minimal(base_size = 12) +
  labs(x = "Spermatogenic Stage", y = "Mean Expression (SCT)",
       title = "FLVCR1 vs FLVCR2: Stage-Specific Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave(file.path(dir_figs, "FLVCR1_FLVCR2_stage_profiles_bar.png"), p2, 
       width = 10, height = 10, dpi = 300)

# ── Calculate fold-changes ──
cat("\nCalculating fold-changes...\n")

fc_results <- flvcr_data %>%
  group_by(Dataset, Gene) %>%
  arrange(Stage) %>%
  mutate(
    FC_vs_earliest = Mean / first(Mean),
    Log2FC_vs_earliest = log2(Mean / first(Mean))
  ) %>%
  ungroup()

# Spermatids vs Spermatogonia fold-change
fc_summary <- flvcr_data %>%
  filter(Stage %in% c("Spermatogonia", "RoundSpermatids", "ElongatingSpermatids")) %>%
  group_by(Dataset, Gene) %>%
  summarize(
    SG_mean = Mean[Stage == "Spermatogonia"][1],
    RS_mean = Mean[Stage == "RoundSpermatids"][1],
    ES_mean = Mean[Stage == "ElongatingSpermatids"][1],
    RS_vs_SG_FC = RS_mean / SG_mean,
    ES_vs_SG_FC = ES_mean / SG_mean,
    Max_spermatid_FC = max(c(RS_vs_SG_FC, ES_vs_SG_FC), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(fc_summary, file.path(dir_out, "FLVCR_fold_changes.csv"), row.names = FALSE)

cat("\nFold-change summary (Spermatids vs Spermatogonia):\n")
print(fc_summary)

# ── Validation criteria check ──
cat("\n=== VALIDATION CRITERIA CHECK ===\n")

flvcr2_fc <- fc_summary %>% filter(Gene == "FLVCR2")

cat("\n✓ PASS CRITERIA: FLVCR2 in spermatids >10× vs spermatogonia\n")
for (i in 1:nrow(flvcr2_fc)) {
  dataset <- flvcr2_fc$Dataset[i]
  fc <- flvcr2_fc$Max_spermatid_FC[i]
  status <- ifelse(fc > 10, "✓ PASS", "✗ FAIL")
  cat(sprintf("  %s: %.1f× %s\n", dataset, fc, status))
}

# ── Heatmap of all genes ──
if (nrow(combined_results) > 0) {
  # Prepare matrix for heatmap
  heatmap_data <- combined_results %>%
    filter(!is.na(Mean)) %>%
    select(Dataset, Stage, Gene, Mean) %>%
    pivot_wider(names_from = Stage, values_from = Mean, values_fill = 0)
  
  # Only plot if we have data
  if (nrow(heatmap_data) > 1) {
    library(pheatmap)
    
    mat <- as.matrix(heatmap_data[, -c(1:2)])
    rownames(mat) <- paste(heatmap_data$Dataset, heatmap_data$Gene, sep = "_")
    
    # Z-score normalize
    mat_scaled <- t(scale(t(mat)))
    
    png(file.path(dir_figs, "FLVCR_context_genes_heatmap.png"), 
        width = 8, height = 10, units = "in", res = 300)
    pheatmap(mat_scaled,
             cluster_cols = FALSE,
             cluster_rows = TRUE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = "Stage Expression Profiles (Z-scored)",
             fontsize = 10)
    dev.off()
    
    cat("\n✓ Heatmap saved\n")
  }
}

# ── Summary statistics ──
summary_stats <- combined_results %>%
  group_by(Gene, Stage) %>%
  summarize(
    N_datasets = n(),
    Mean_expr = mean(Mean),
    SD_expr = sd(Mean),
    .groups = "drop"
  )

write.csv(summary_stats, file.path(dir_out, "FLVCR_stage_summary_stats.csv"), 
          row.names = FALSE)

cat("\n=== Stage Profiles Complete ===\n")
cat("Results saved to:", dir_out, "\n")
cat("Figures saved to:", dir_figs, "\n")
cat("\nNext step: Run 05_WGCNA_FLVCR1_FLVCR2.R\n")

