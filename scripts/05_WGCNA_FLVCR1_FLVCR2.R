# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 05: WGCNA Co-expression Network Analysis
# ============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cat("=== FLVCR Pipeline: WGCNA Co-expression Analysis ===\n\n")

# ── Configuration ──
dir_in <- "data_intermediate"
dir_out <- "data_results"
dir_figs <- "figures"

# ── Load integrated human dataset ──
cat("Loading integrated human dataset...\n")

int_file <- file.path(dir_in, "integrated_human.rds")
if (!file.exists(int_file)) {
  stop("Integrated dataset not found. Please run script 03 first.")
}

obj <- readRDS(int_file)

# ── Filter to germ cells ──
germ_cells <- c("Spermatogonia", "Spermatocytes", "RoundSpermatids", "ElongatingSpermatids")
obj_germ <- subset(obj, subset = celltype %in% germ_cells)

cat("Germ cells for WGCNA:", ncol(obj_germ), "\n")

# ── Create pseudo-bulk profiles per stage-sample ──
cat("\nCreating pseudo-bulk profiles...\n")

# Add unique sample-stage identifier
obj_germ$stage_sample <- paste(obj_germ$celltype, obj_germ$dataset, sep = "_")

# Aggregate cells by stage-sample
Idents(obj_germ) <- "stage_sample"
pseudobulk <- AverageExpression(obj_germ, assays = "SCT", slot = "data")$SCT

cat("Pseudo-bulk samples:", ncol(pseudobulk), "\n")
cat("Genes:", nrow(pseudobulk), "\n")

# ── Prepare expression matrix for WGCNA ──
# WGCNA expects samples as rows, genes as columns
datExpr <- t(as.matrix(pseudobulk))

# Filter genes: keep top variable genes
vars <- apply(datExpr, 2, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:min(5000, ncol(datExpr))]
datExpr <- datExpr[, top_genes]

cat("WGCNA input:", nrow(datExpr), "samples ×", ncol(datExpr), "genes\n")

# ── Check for missing values and outliers ──
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("Removing problematic genes/samples...\n")
  if (sum(!gsg$goodGenes) > 0) {
    datExpr <- datExpr[, gsg$goodGenes]
  }
  if (sum(!gsg$goodSamples) > 0) {
    datExpr <- datExpr[gsg$goodSamples, ]
  }
}

# ── Choose soft-thresholding power ──
cat("\nChoosing soft-thresholding power...\n")

powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free topology fit
png(file.path(dir_figs, "WGCNA_soft_threshold.png"), width = 10, height = 5, 
    units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# Select power
softPower <- ifelse(is.na(sft$powerEstimate), 8, sft$powerEstimate)
cat("Selected soft power:", softPower, "\n")

# ── Construct network and detect modules ──
cat("\nConstructing co-expression network...\n")

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# Convert labels to colors
moduleColors <- labels2colors(net$colors)
table(moduleColors)

cat("\nNumber of modules:", length(unique(moduleColors)) - 1, "\n")  # -1 for grey

# ── Plot dendrogram ──
png(file.path(dir_figs, "WGCNA_dendrogram_modules.png"), width = 12, height = 6, 
    units = "in", res = 300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ── Calculate module eigengenes ──
MEs <- net$MEs
colnames(MEs) <- gsub("ME", "", colnames(MEs))

# ── Find FLVCR1/2 modules ──
cat("\n=== FLVCR1/2 Module Assignment ===\n")

gene_modules <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

flvcr_modules <- gene_modules %>%
  filter(Gene %in% c("FLVCR1", "FLVCR2"))

if (nrow(flvcr_modules) > 0) {
  print(flvcr_modules)
  
  flvcr1_mod <- flvcr_modules$Module[flvcr_modules$Gene == "FLVCR1"]
  flvcr2_mod <- flvcr_modules$Module[flvcr_modules$Gene == "FLVCR2"]
  
  if (length(flvcr1_mod) > 0) {
    cat("\nFLVCR1 module:", flvcr1_mod, "\n")
  } else {
    cat("\nFLVCR1 not found in filtered gene set\n")
    flvcr1_mod <- NA
  }
  
  if (length(flvcr2_mod) > 0) {
    cat("FLVCR2 module:", flvcr2_mod, "\n")
  } else {
    cat("FLVCR2 not found in filtered gene set\n")
    flvcr2_mod <- NA
  }
} else {
  cat("FLVCR1/2 not found in network. They may have been filtered out.\n")
  flvcr1_mod <- NA
  flvcr2_mod <- NA
}

# ── Calculate module membership (kME) ──
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), nrow(datExpr)))

colnames(geneModuleMembership) <- paste0("MM.", colnames(geneModuleMembership))
colnames(MMPvalue) <- paste0("p.MM.", colnames(MMPvalue))

# ── Extract hub genes for each module ──
cat("\nExtracting hub genes...\n")

hub_genes_list <- list()

for (mod in unique(moduleColors)) {
  if (mod == "grey") next  # Skip unassigned genes
  
  mod_genes <- gene_modules$Gene[gene_modules$Module == mod]
  mm_col <- paste0("MM.", mod)
  
  if (mm_col %in% colnames(geneModuleMembership)) {
    kME <- geneModuleMembership[mod_genes, mm_col]
    names(kME) <- mod_genes
    
    # Top 50 hub genes
    hub_genes <- names(sort(kME, decreasing = TRUE))[1:min(50, length(kME))]
    hub_genes_list[[mod]] <- hub_genes
  }
}

# ── Analyze FLVCR2 module ──
if (!is.na(flvcr2_mod) && flvcr2_mod != "grey") {
  cat("\n=== FLVCR2 Module Analysis ===\n")
  cat("Module:", flvcr2_mod, "\n")
  
  flvcr2_hub_genes <- hub_genes_list[[flvcr2_mod]]
  cat("Hub genes (top 50):\n")
  print(head(flvcr2_hub_genes, 20))
  
  # Check for spermiogenesis markers
  spermiogenesis_markers <- c("PRM1", "PRM2", "TNP1", "TNP2", "CRISP2", 
                              "ACR", "SPACA7", "AKAP4", "ODF2", "DNAJC5B",
                              "ACRV1", "TEKT1", "SPAG16")
  
  overlap <- intersect(gene_modules$Gene[gene_modules$Module == flvcr2_mod], 
                       spermiogenesis_markers)
  cat("\nSpermiogenesis markers in FLVCR2 module:\n")
  print(overlap)
  
  # Save FLVCR2 module genes
  flvcr2_module_genes <- gene_modules %>%
    filter(Module == flvcr2_mod) %>%
    arrange(desc(Gene %in% flvcr2_hub_genes))
  
  write.csv(flvcr2_module_genes, 
            file.path(dir_out, "flvcr2_module_genes.csv"), 
            row.names = FALSE)
}

# ── Analyze FLVCR1 module ──
if (!is.na(flvcr1_mod) && flvcr1_mod != "grey") {
  cat("\n=== FLVCR1 Module Analysis ===\n")
  cat("Module:", flvcr1_mod, "\n")
  
  flvcr1_hub_genes <- hub_genes_list[[flvcr1_mod]]
  cat("Hub genes (top 50):\n")
  print(head(flvcr1_hub_genes, 20))
  
  # Check for translation/stress/heme markers
  early_markers <- c("EEF1A1", "EEF2", "RPL3", "RPS6", "HSP90AB1", 
                     "HSPA8", "FECH", "HMBS", "ALAS2")
  
  overlap <- intersect(gene_modules$Gene[gene_modules$Module == flvcr1_mod], 
                       early_markers)
  cat("\nTranslation/stress/heme markers in FLVCR1 module:\n")
  print(overlap)
  
  # Save FLVCR1 module genes
  flvcr1_module_genes <- gene_modules %>%
    filter(Module == flvcr1_mod) %>%
    arrange(desc(Gene %in% flvcr1_hub_genes))
  
  write.csv(flvcr1_module_genes, 
            file.path(dir_out, "flvcr1_module_genes.csv"), 
            row.names = FALSE)
}

# ── Save all results ──
write.csv(gene_modules, file.path(dir_out, "WGCNA_all_gene_modules.csv"), 
          row.names = FALSE)

saveRDS(list(
  net = net,
  moduleColors = moduleColors,
  MEs = MEs,
  geneModuleMembership = geneModuleMembership,
  MMPvalue = MMPvalue,
  hub_genes = hub_genes_list,
  flvcr1_module = flvcr1_mod,
  flvcr2_module = flvcr2_mod
), file.path(dir_out, "WGCNA_results.rds"))

cat("\n=== WGCNA Analysis Complete ===\n")
cat("Results saved to:", dir_out, "\n")
cat("Figures saved to:", dir_figs, "\n")
cat("\nNext step: Run 06_GO_enrichment.R\n")

