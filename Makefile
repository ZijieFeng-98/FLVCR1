# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Makefile for automated execution
# ============================================================================

.PHONY: all env download scqc integrate stageprofiles wgcna enrich bulk ortholog report clean help

# Default target
all: env download scqc integrate stageprofiles wgcna enrich bulk ortholog report

# Setup environment and install packages
env:
	@echo "=== Setting up R environment ==="
	Rscript scripts/00_setup_env.R

# Download GEO datasets
download:
	@echo "=== Downloading GEO datasets ==="
	Rscript scripts/01_download_GEO.R

# Single-cell QC and annotation
scqc:
	@echo "=== Running scRNA-seq QC and annotation ==="
	Rscript scripts/02_scRNA_QC_annotate.R

# Integration and pseudotime analysis
integrate:
	@echo "=== Running integration and pseudotime analysis ==="
	Rscript scripts/03_sc_integration_pseudotime.R

# Stage-resolved expression profiles
stageprofiles:
	@echo "=== Calculating stage-specific profiles ==="
	Rscript scripts/04_stage_profiles_FLVCR.R

# WGCNA co-expression network analysis
wgcna:
	@echo "=== Running WGCNA analysis ==="
	Rscript scripts/05_WGCNA_FLVCR1_FLVCR2.R

# GO and pathway enrichment
enrich:
	@echo "=== Running enrichment analysis ==="
	Rscript scripts/06_GO_enrichment.R

# Bulk validation with composition adjustment
bulk:
	@echo "=== Running bulk validation ==="
	Rscript scripts/07_bulk_validation_composition_adjust.R

# Cross-species ortholog analysis
ortholog:
	@echo "=== Running ortholog analysis ==="
	Rscript scripts/08_cross_species_orthologs.R

# Generate final report
report:
	@echo "=== Generating Quarto report ==="
	quarto render scripts/09_report_quarto.qmd

# Clean intermediate files (use with caution!)
clean:
	@echo "=== Cleaning intermediate files ==="
	@echo "WARNING: This will delete processed data!"
	@read -p "Are you sure? (y/N): " confirm && [ "$$confirm" = "y" ] && \
		rm -rf data_intermediate/* data_results/* figures/* || echo "Cancelled"

# Clean everything including raw data (DANGEROUS!)
clean-all:
	@echo "=== Cleaning ALL data files ==="
	@echo "WARNING: This will delete ALL data including raw downloads!"
	@read -p "Are you sure? (y/N): " confirm && [ "$$confirm" = "y" ] && \
		rm -rf data_raw/* data_intermediate/* data_results/* figures/* || echo "Cancelled"

# Display help
help:
	@echo "FLVCR Validation Pipeline - Available targets:"
	@echo ""
	@echo "  make all          - Run complete pipeline"
	@echo "  make env          - Setup R environment"
	@echo "  make download     - Download GEO datasets"
	@echo "  make scqc         - scRNA-seq QC & annotation"
	@echo "  make integrate    - Integration & pseudotime"
	@echo "  make stageprofiles- Stage expression profiles"
	@echo "  make wgcna        - Co-expression networks"
	@echo "  make enrich       - Pathway enrichment"
	@echo "  make bulk         - Bulk validation"
	@echo "  make ortholog     - Cross-species analysis"
	@echo "  make report       - Generate final report"
	@echo "  make clean        - Remove intermediate files"
	@echo "  make clean-all    - Remove ALL data files"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Example usage:"
	@echo "  make env          # First-time setup"
	@echo "  make download     # Get data"
	@echo "  make all          # Run full pipeline"

# Individual script shortcuts for debugging
script-%:
	@echo "=== Running script $* ==="
	Rscript scripts/$*.R

