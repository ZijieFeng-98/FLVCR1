# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 00: Environment Setup & Package Installation
# ============================================================================

cat("=== FLVCR Pipeline: Environment Setup ===\n\n")

# Initialize renv for reproducibility (optional but recommended)
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Option to initialize renv project
# Uncomment if you want full reproducibility tracking:
# renv::init()

cat("Installing required packages...\n")

# Core packages list organized by category
pkgs <- c(
  # ── Data I/O & manipulation ──
  "GEOquery",
  "Matrix",
  "data.table",
  "readr",
  "jsonlite",
  "BiocManager",
  
  # ── Single-cell RNA-seq ──
  "Seurat",
  "SeuratObject",
  "SingleCellExperiment",
  "scater",
  "scran",
  "batchelor",
  "SeuratDisk",
  "uwot",
  "AUCell",
  "GSEABase",
  
  # ── Trajectory & pseudotime ──
  "slingshot",
  "tradeSeq",
  # "monocle3",  # Optional - BiocManager install
  
  # ── Bulk RNA-seq DE & deconvolution ──
  "edgeR",
  "DESeq2",
  "limma",
  "MuSiC",
  # "BisqueRNA",  # Optional
  
  # ── Network analysis & enrichment ──
  "WGCNA",
  "biomaRt",
  "clusterProfiler",
  "fgsea",
  "msigdbr",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "enrichplot",
  "GO.db",
  "DOSE",
  
  # ── Visualization ──
  "ggplot2",
  "cowplot",
  "pheatmap",
  "ComplexHeatmap",
  "ggrepel",
  "patchwork",
  "viridis",
  "RColorBrewer",
  
  # ── Utilities ──
  "tidyverse",
  "magrittr",
  "here",
  "scales"
)

# Separate Bioconductor packages
bioc_pkgs <- c(
  "GEOquery",
  "SingleCellExperiment",
  "scater",
  "scran",
  "batchelor",
  "AUCell",
  "GSEABase",
  "slingshot",
  "tradeSeq",
  "edgeR",
  "DESeq2",
  "limma",
  "biomaRt",
  "clusterProfiler",
  "fgsea",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "enrichplot",
  "GO.db",
  "DOSE",
  "ComplexHeatmap"
)

# Install BiocManager first if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install CRAN packages
cran_pkgs <- setdiff(pkgs, bioc_pkgs)
to_install_cran <- setdiff(cran_pkgs, installed.packages()[, 1])
if (length(to_install_cran) > 0) {
  cat("Installing CRAN packages:", paste(to_install_cran, collapse = ", "), "\n")
  install.packages(to_install_cran, dependencies = TRUE)
}

# Install Bioconductor packages
to_install_bioc <- setdiff(bioc_pkgs, installed.packages()[, 1])
if (length(to_install_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(to_install_bioc, collapse = ", "), "\n")
  BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)
}

# Optional: Install MuSiC from GitHub if not on CRAN/Bioc
if (!"MuSiC" %in% installed.packages()[, 1]) {
  cat("Installing MuSiC from GitHub...\n")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("xuranw/MuSiC")
}

cat("\n=== Package installation complete ===\n")

# Verify critical packages
critical <- c("Seurat", "GEOquery", "WGCNA", "clusterProfiler", "slingshot")
missing <- setdiff(critical, installed.packages()[, 1])
if (length(missing) > 0) {
  warning("Critical packages still missing: ", paste(missing, collapse = ", "))
} else {
  cat("All critical packages verified ✓\n")
}

# Optional: Create renv snapshot
# renv::snapshot()

# Session info for reproducibility
cat("\n=== Session Info ===\n")
print(sessionInfo())

cat("\nSetup complete! You can now run the analysis scripts.\n")

