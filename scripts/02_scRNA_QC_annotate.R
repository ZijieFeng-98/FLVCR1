# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 02: scRNA-seq QC & Cell Type Annotation
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(AUCell)
  library(GSEABase)
  library(patchwork)
})

cat("=== FLVCR Pipeline: scRNA-seq QC & Annotation ===\n\n")

# ── Configuration ──
dir_raw <- "data_raw"
dir_out <- "data_intermediate"
dir_figs <- "figures"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figs, showWarnings = FALSE, recursive = TRUE)

# ── Cell type marker genes ──
markers <- list(
  Spermatogonia = c("UTF1", "ID4", "GFRA1", "ZBTB16", "KIT", "UCHL1", "NANOS2", "DMRT1"),
  Spermatocytes = c("SYCP3", "SYCP1", "SYCP2", "TEX101", "PIWIL1", "HORMAD1", "DMC1"),
  RoundSpermatids = c("PRM1", "PRM2", "TNP1", "TNP2", "ACR", "SPACA7", "ACRV1"),
  ElongatingSpermatids = c("PRM2", "AKAP4", "ODF2", "DNAJC5B", "TEKT1", "SPAG16"),
  Sertoli = c("AMH", "SOX9", "GATA4", "CTSL", "CLU", "FATE1"),
  Leydig = c("STAR", "CYP17A1", "HSD3B1", "INSL3"),
  Peritubular = c("ACTA2", "MYH11", "TAGLN", "DES"),
  Endothelial = c("PECAM1", "VWF", "CDH5"),
  Macrophages = c("CD68", "CD163", "LYZ", "C1QA")
)

# ── Helper function: Load count matrix ──
load_count_matrix <- function(geo_id, data_dir = "data_raw") {
  cat("Loading", geo_id, "...\n")
  
  geo_dir <- file.path(data_dir, geo_id)
  
  # Try different file formats
  # 1. 10x format
  if (file.exists(file.path(geo_dir, "matrix.mtx")) || 
      file.exists(file.path(geo_dir, "matrix.mtx.gz"))) {
    cat("  Detected 10x format\n")
    mat <- Read10X(geo_dir)
    return(CreateSeuratObject(counts = mat, project = geo_id))
  }
  
  # 2. h5 format
  h5_files <- list.files(geo_dir, pattern = "\\.h5$", full.names = TRUE)
  if (length(h5_files) > 0) {
    cat("  Detected h5 format:", basename(h5_files[1]), "\n")
    mat <- Read10X_h5(h5_files[1])
    return(CreateSeuratObject(counts = mat, project = geo_id))
  }
  
  # 3. Text format (example - adjust based on actual format)
  txt_files <- list.files(geo_dir, pattern = "counts.*\\.(txt|tsv)", full.names = TRUE, ignore.case = TRUE)
  if (length(txt_files) > 0) {
    cat("  Detected text format:", basename(txt_files[1]), "\n")
    mat <- read.table(txt_files[1], header = TRUE, row.names = 1, sep = "\t")
    mat <- as(as.matrix(mat), "dgCMatrix")
    return(CreateSeuratObject(counts = mat, project = geo_id))
  }
  
  stop("Could not find count matrix in ", geo_dir)
}

# ── Helper function: QC & preprocessing ──
process_seurat <- function(obj, dataset_name, species = "human") {
  cat("\nProcessing", dataset_name, "...\n")
  
  # Add metadata
  obj$dataset <- dataset_name
  obj$species <- species
  
  # Calculate QC metrics
  mt_pattern <- if (species == "human") "^MT-" else "^mt-"
  ribo_pattern <- if (species == "human") "^RP[SL]" else "^Rp[sl]"
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  
  # QC plots before filtering
  p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  ggsave(file.path(dir_figs, paste0("QC_before_", dataset_name, ".png")), p1, width = 12, height = 4)
  
  cat("  Cells before QC:", ncol(obj), "\n")
  
  # Filter cells - adjust thresholds based on your data
  obj <- subset(obj, subset = 
                  nFeature_RNA >= 200 & 
                  nFeature_RNA <= 6000 & 
                  nCount_RNA >= 500 &
                  nCount_RNA <= 50000 &
                  percent.mt < 15)
  
  cat("  Cells after QC:", ncol(obj), "\n")
  
  # Normalize
  obj <- SCTransform(obj, vars.to.regress = c("percent.mt"), verbose = FALSE)
  
  # Dimensionality reduction
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
  
  # UMAP colored by clusters
  p2 <- DimPlot(obj, reduction = "umap", label = TRUE) + 
    ggtitle(paste(dataset_name, "- Clusters"))
  ggsave(file.path(dir_figs, paste0("UMAP_clusters_", dataset_name, ".png")), p2, width = 8, height = 6)
  
  return(obj)
}

# ── Helper function: AUCell-based cell type annotation ──
annotate_celltypes <- function(obj, markers, species = "human") {
  cat("Annotating cell types using AUCell...\n")
  
  # Convert mouse gene names if needed
  if (species == "mouse") {
    markers <- lapply(markers, function(genes) {
      tolower(substring(genes, 1, 1)) %>% 
        paste0(substring(genes, 2))
    })
  }
  
  # Prepare gene sets
  gene_sets <- lapply(names(markers), function(ct) {
    GeneSet(intersect(markers[[ct]], rownames(obj)), setName = ct)
  })
  names(gene_sets) <- names(markers)
  
  # Calculate AUCell scores
  cells_rankings <- AUCell_buildRankings(GetAssayData(obj, slot = "data"), 
                                         plotStats = FALSE, verbose = FALSE)
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, verbose = FALSE)
  
  # Add scores to Seurat object
  auc_mtx <- getAUC(cells_AUC)
  for (ct in rownames(auc_mtx)) {
    obj[[paste0("AUC_", ct)]] <- auc_mtx[ct, ]
  }
  
  # Assign cell type based on maximum AUC score
  max_auc <- apply(auc_mtx, 2, which.max)
  obj$celltype <- names(markers)[max_auc]
  obj$celltype_confidence <- apply(auc_mtx, 2, max)
  
  # Visualize
  p1 <- DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE) +
    ggtitle(paste(unique(obj$dataset), "- Cell Types"))
  
  # Feature plots for key markers
  key_markers <- c(
    if (species == "human") c("UTF1", "SYCP3", "PRM1", "AMH") 
    else c("Utf1", "Sycp3", "Prm1", "Amh")
  )
  key_markers <- intersect(key_markers, rownames(obj))
  p2 <- FeaturePlot(obj, features = key_markers[1:min(4, length(key_markers))], ncol = 2)
  
  combined <- p1 / p2
  ggsave(file.path(dir_figs, paste0("UMAP_celltypes_", unique(obj$dataset), ".png")), 
         combined, width = 10, height = 12)
  
  # Summary table
  cat("\nCell type distribution:\n")
  print(table(obj$celltype))
  
  return(obj)
}

# ── Main processing pipeline ──

# Process human datasets
cat("\n=== Processing Human Datasets ===\n")

# GSE120508 (example - adjust based on actual file structure)
tryCatch({
  obj_120508 <- load_count_matrix("GSE120508", dir_raw)
  obj_120508 <- process_seurat(obj_120508, "GSE120508", species = "human")
  obj_120508 <- annotate_celltypes(obj_120508, markers, species = "human")
  saveRDS(obj_120508, file.path(dir_out, "GSE120508_qc_annotated.rds"))
  cat("✓ GSE120508 complete\n")
}, error = function(e) {
  cat("✗ Error processing GSE120508:", e$message, "\n")
  cat("  This dataset may require manual download. See data_raw/ notes.\n")
})

# GSE109037
tryCatch({
  obj_109037 <- load_count_matrix("GSE109037", dir_raw)
  obj_109037 <- process_seurat(obj_109037, "GSE109037", species = "human")
  obj_109037 <- annotate_celltypes(obj_109037, markers, species = "human")
  saveRDS(obj_109037, file.path(dir_out, "GSE109037_qc_annotated.rds"))
  cat("✓ GSE109037 complete\n")
}, error = function(e) {
  cat("✗ Error processing GSE109037:", e$message, "\n")
})

# Process mouse dataset
cat("\n=== Processing Mouse Dataset ===\n")

# Convert human markers to mouse
markers_mouse <- markers

tryCatch({
  obj_107644 <- load_count_matrix("GSE107644", dir_raw)
  obj_107644 <- process_seurat(obj_107644, "GSE107644", species = "mouse")
  obj_107644 <- annotate_celltypes(obj_107644, markers_mouse, species = "mouse")
  saveRDS(obj_107644, file.path(dir_out, "GSE107644_qc_annotated.rds"))
  cat("✓ GSE107644 complete\n")
}, error = function(e) {
  cat("✗ Error processing GSE107644:", e$message, "\n")
})

# ── Summary ──
cat("\n=== Processing Complete ===\n")
cat("Annotated objects saved to:", dir_out, "\n")
cat("QC and annotation plots saved to:", dir_figs, "\n")
cat("\nNext step: Run 03_sc_integration_pseudotime.R\n")

