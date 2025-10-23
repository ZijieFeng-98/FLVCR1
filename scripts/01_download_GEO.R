# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 01: Download GEO Datasets
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(readr)
})

cat("=== FLVCR Pipeline: GEO Data Download ===\n\n")

# ── Configuration ──
dir_out <- "data_raw"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Target datasets for validation
geo_ids <- c(
  "GSE109037",  # Human testis scRNA-seq (Wang et al. 2018)
  "GSE120508",  # Human testis scRNA-seq (Guo et al. 2018)
  "GSE107644"   # Mouse testis scRNA-seq (Hermann et al. 2018)
  # Add additional bulk datasets as needed:
  # "GSE145467", # Example bulk infertile vs control
)

cat("Target datasets:", paste(geo_ids, collapse = ", "), "\n\n")

# ── Download function ──
download_geo <- function(geo_id, out_dir) {
  cat("Processing", geo_id, "...\n")
  
  tryCatch({
    # Download series matrix
    gset <- getGEO(geo_id, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = out_dir)
    
    # Save metadata
    if (length(gset) > 0) {
      saveRDS(gset, file = file.path(out_dir, paste0(geo_id, "_series_matrix.rds")))
      cat("  ✓ Downloaded series matrix for", geo_id, "\n")
    }
    
    # Get supplementary files (often contains count matrices)
    supp_files <- getGEOSuppFiles(geo_id, baseDir = out_dir, makeDirectory = TRUE)
    cat("  ✓ Downloaded", nrow(supp_files), "supplementary file(s)\n")
    
    # Extract if compressed
    supp_dir <- file.path(out_dir, geo_id)
    gz_files <- list.files(supp_dir, pattern = "\\.gz$", full.names = TRUE)
    
    for (gz in gz_files) {
      out_file <- sub("\\.gz$", "", gz)
      if (!file.exists(out_file)) {
        cat("  Extracting:", basename(gz), "\n")
        R.utils::gunzip(gz, out_file, remove = FALSE, overwrite = TRUE)
      }
    }
    
    return(TRUE)
    
  }, error = function(e) {
    cat("  ✗ Error downloading", geo_id, ":", e$message, "\n")
    return(FALSE)
  })
}

# ── Download all datasets ──
cat("Starting downloads...\n")
results <- sapply(geo_ids, download_geo, out_dir = dir_out)

# ── Summary ──
cat("\n=== Download Summary ===\n")
cat("Successful:", sum(results), "/", length(results), "\n")
cat("Failed:", sum(!results), "\n")

if (sum(!results) > 0) {
  cat("\nFailed datasets:", paste(names(results)[!results], collapse = ", "), "\n")
}

# ── Instructions for manual download ──
cat("\n=== IMPORTANT NOTES ===\n")
cat("Some large scRNA-seq datasets may require manual download:\n")
cat("1. Visit GEO website: https://www.ncbi.nlm.nih.gov/geo/\n")
cat("2. Search for dataset (e.g., GSE120508)\n")
cat("3. Download processed count matrices from 'Supplementary file' section\n")
cat("4. Place files in: data_raw/<GSE_ID>/\n")
cat("\nExpected file formats:\n")
cat("  - 10x format: matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz\n")
cat("  - h5 format: filtered_feature_bc_matrix.h5\n")
cat("  - Text: counts.txt.gz or similar\n")

# ── Save download log ──
log_df <- data.frame(
  GEO_ID = geo_ids,
  Downloaded = results,
  Timestamp = Sys.time()
)
write_tsv(log_df, file.path(dir_out, "download_log.tsv"))

cat("\nDownload log saved to:", file.path(dir_out, "download_log.tsv"), "\n")
cat("\nScript complete!\n")

