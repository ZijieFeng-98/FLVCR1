# ============================================================================
# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline
# Script 06: GO & Pathway Enrichment Analysis
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(msigdbr)
  library(dplyr)
  library(DOSE)
})

cat("=== FLVCR Pipeline: GO & Pathway Enrichment ===\n\n")

# ── Configuration ──
dir_in <- "data_results"
dir_out <- "data_results"
dir_figs <- "figures"

# ── Helper function: Convert symbols to ENTREZ ──
symbols_to_entrez <- function(symbols, organism = "org.Hs.eg.db") {
  orgdb <- get(organism)
  mapping <- AnnotationDbi::select(orgdb, 
                                   keys = symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  mapping <- mapping[!is.na(mapping$ENTREZID), ]
  return(unique(mapping))
}

# ── Helper function: Run comprehensive enrichment ──
run_enrichment <- function(gene_symbols, module_name, species = "human") {
  cat("\n=== Enrichment for", module_name, "===\n")
  
  # Convert to ENTREZ
  orgdb <- if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
  gene_mapping <- symbols_to_entrez(gene_symbols, orgdb)
  entrez_ids <- unique(gene_mapping$ENTREZID)
  
  cat("Input genes:", length(gene_symbols), "\n")
  cat("Mapped to ENTREZ:", length(entrez_ids), "\n")
  
  if (length(entrez_ids) < 5) {
    warning("Too few genes for enrichment (need >= 5)")
    return(NULL)
  }
  
  results <- list()
  
  # ── GO Biological Process ──
  tryCatch({
    ego_bp <- enrichGO(
      gene = entrez_ids,
      OrgDb = orgdb,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
      results$GO_BP <- as.data.frame(ego_bp)
      cat("  GO BP terms:", nrow(results$GO_BP), "\n")
    } else {
      cat("  No significant GO BP terms\n")
    }
  }, error = function(e) {
    cat("  GO BP error:", e$message, "\n")
  })
  
  # ── GO Molecular Function ──
  tryCatch({
    ego_mf <- enrichGO(
      gene = entrez_ids,
      OrgDb = orgdb,
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    if (!is.null(ego_mf) && nrow(ego_mf) > 0) {
      results$GO_MF <- as.data.frame(ego_mf)
      cat("  GO MF terms:", nrow(results$GO_MF), "\n")
    }
  }, error = function(e) {
    cat("  GO MF error:", e$message, "\n")
  })
  
  # ── KEGG Pathway ──
  tryCatch({
    kegg_organism <- if (species == "human") "hsa" else "mmu"
    kegg <- enrichKEGG(
      gene = entrez_ids,
      organism = kegg_organism,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if (!is.null(kegg) && nrow(kegg) > 0) {
      kegg <- setReadable(kegg, orgdb, keyType = "ENTREZID")
      results$KEGG <- as.data.frame(kegg)
      cat("  KEGG pathways:", nrow(results$KEGG), "\n")
    }
  }, error = function(e) {
    cat("  KEGG error:", e$message, "\n")
  })
  
  # ── MSigDB Hallmarks ──
  tryCatch({
    species_name <- if (species == "human") "Homo sapiens" else "Mus musculus"
    hallmarks <- msigdbr(species = species_name, category = "H")
    
    hallmark_list <- split(hallmarks$entrez_gene, hallmarks$gs_name)
    
    enrich_hallmark <- enricher(
      gene = entrez_ids,
      TERM2GENE = hallmarks[, c("gs_name", "entrez_gene")],
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
      results$Hallmarks <- as.data.frame(enrich_hallmark)
      cat("  Hallmark gene sets:", nrow(results$Hallmarks), "\n")
    }
  }, error = function(e) {
    cat("  Hallmarks error:", e$message, "\n")
  })
  
  # ── MSigDB GO (C5) for specific reproductive terms ──
  tryCatch({
    species_name <- if (species == "human") "Homo sapiens" else "Mus musculus"
    c5_bp <- msigdbr(species = species_name, category = "C5", subcategory = "GO:BP")
    
    enrich_c5 <- enricher(
      gene = entrez_ids,
      TERM2GENE = c5_bp[, c("gs_name", "entrez_gene")],
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if (!is.null(enrich_c5) && nrow(enrich_c5) > 0) {
      results$MSigDB_C5 <- as.data.frame(enrich_c5)
      cat("  MSigDB C5 terms:", nrow(results$MSigDB_C5), "\n")
    }
  }, error = function(e) {
    cat("  MSigDB C5 error:", e$message, "\n")
  })
  
  return(results)
}

# ── Load WGCNA results ──
cat("Loading WGCNA results...\n")

wgcna_file <- file.path(dir_in, "WGCNA_results.rds")
if (!file.exists(wgcna_file)) {
  stop("WGCNA results not found. Please run script 05 first.")
}

wgcna_results <- readRDS(wgcna_file)

# ── Enrichment for FLVCR2 module ──
flvcr2_genes_file <- file.path(dir_in, "flvcr2_module_genes.csv")
if (file.exists(flvcr2_genes_file)) {
  flvcr2_genes <- read.csv(flvcr2_genes_file)
  flvcr2_enrich <- run_enrichment(flvcr2_genes$Gene, "FLVCR2_module", species = "human")
  
  # Save results
  if (!is.null(flvcr2_enrich)) {
    saveRDS(flvcr2_enrich, file.path(dir_out, "FLVCR2_module_enrichment.rds"))
    
    # Save as tables
    for (name in names(flvcr2_enrich)) {
      write.csv(flvcr2_enrich[[name]], 
                file.path(dir_out, paste0("FLVCR2_module_", name, ".csv")),
                row.names = FALSE)
    }
    
    # ── Visualize FLVCR2 enrichment ──
    if ("GO_BP" %in% names(flvcr2_enrich) && nrow(flvcr2_enrich$GO_BP) > 0) {
      ego_bp_obj <- new("enrichResult", result = flvcr2_enrich$GO_BP)
      
      # Dot plot
      p1 <- dotplot(ego_bp_obj, showCategory = 20) +
        ggtitle("FLVCR2 Module: GO Biological Process")
      ggsave(file.path(dir_figs, "FLVCR2_GO_BP_dotplot.png"), p1, 
             width = 10, height = 8, dpi = 300)
      
      # Bar plot
      p2 <- barplot(ego_bp_obj, showCategory = 15) +
        ggtitle("FLVCR2 Module: Top GO BP Terms")
      ggsave(file.path(dir_figs, "FLVCR2_GO_BP_barplot.png"), p2, 
             width = 10, height = 6, dpi = 300)
    }
  }
} else {
  cat("FLVCR2 module genes not found. Skipping.\n")
}

# ── Enrichment for FLVCR1 module ──
flvcr1_genes_file <- file.path(dir_in, "flvcr1_module_genes.csv")
if (file.exists(flvcr1_genes_file)) {
  flvcr1_genes <- read.csv(flvcr1_genes_file)
  flvcr1_enrich <- run_enrichment(flvcr1_genes$Gene, "FLVCR1_module", species = "human")
  
  # Save results
  if (!is.null(flvcr1_enrich)) {
    saveRDS(flvcr1_enrich, file.path(dir_out, "FLVCR1_module_enrichment.rds"))
    
    # Save as tables
    for (name in names(flvcr1_enrich)) {
      write.csv(flvcr1_enrich[[name]], 
                file.path(dir_out, paste0("FLVCR1_module_", name, ".csv")),
                row.names = FALSE)
    }
    
    # ── Visualize FLVCR1 enrichment ──
    if ("GO_BP" %in% names(flvcr1_enrich) && nrow(flvcr1_enrich$GO_BP) > 0) {
      ego_bp_obj <- new("enrichResult", result = flvcr1_enrich$GO_BP)
      
      p3 <- dotplot(ego_bp_obj, showCategory = 20) +
        ggtitle("FLVCR1 Module: GO Biological Process")
      ggsave(file.path(dir_figs, "FLVCR1_GO_BP_dotplot.png"), p3, 
             width = 10, height = 8, dpi = 300)
      
      p4 <- barplot(ego_bp_obj, showCategory = 15) +
        ggtitle("FLVCR1 Module: Top GO BP Terms")
      ggsave(file.path(dir_figs, "FLVCR1_GO_BP_barplot.png"), p4, 
             width = 10, height = 6, dpi = 300)
    }
  }
} else {
  cat("FLVCR1 module genes not found. Skipping.\n")
}

# ── Compare FLVCR1 vs FLVCR2 modules ──
if (exists("flvcr1_enrich") && exists("flvcr2_enrich")) {
  cat("\n=== Comparing FLVCR1 vs FLVCR2 Modules ===\n")
  
  if ("GO_BP" %in% names(flvcr1_enrich) && "GO_BP" %in% names(flvcr2_enrich)) {
    # Combine top terms
    flvcr1_top <- flvcr1_enrich$GO_BP %>%
      arrange(p.adjust) %>%
      head(10) %>%
      mutate(Module = "FLVCR1")
    
    flvcr2_top <- flvcr2_enrich$GO_BP %>%
      arrange(p.adjust) %>%
      head(10) %>%
      mutate(Module = "FLVCR2")
    
    combined <- rbind(
      flvcr1_top[, c("Description", "p.adjust", "Module")],
      flvcr2_top[, c("Description", "p.adjust", "Module")]
    )
    
    combined$log10padj <- -log10(combined$p.adjust)
    
    p_compare <- ggplot(combined, aes(x = reorder(Description, log10padj), 
                                      y = log10padj, fill = Module)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("FLVCR1" = "#E64B35", "FLVCR2" = "#4DBBD5")) +
      theme_minimal(base_size = 10) +
      labs(x = "", y = "-log10(adjusted p-value)",
           title = "FLVCR1 vs FLVCR2 Module: Top GO BP Terms") +
      theme(legend.position = "top")
    
    ggsave(file.path(dir_figs, "FLVCR1_vs_FLVCR2_GO_comparison.png"), p_compare,
           width = 12, height = 10, dpi = 300)
  }
}

# ── Validation criteria check ──
cat("\n=== VALIDATION CRITERIA CHECK ===\n")

# Check for spermiogenesis terms in FLVCR2
if (exists("flvcr2_enrich") && "GO_BP" %in% names(flvcr2_enrich)) {
  spermiogenesis_terms <- c("spermatid", "sperm", "acrosome", "flagellum", 
                            "cilium", "axoneme", "chromatin condensation")
  
  flvcr2_descriptions <- tolower(flvcr2_enrich$GO_BP$Description)
  
  matches <- sapply(spermiogenesis_terms, function(term) {
    any(grepl(term, flvcr2_descriptions))
  })
  
  cat("\n✓ FLVCR2 module spermiogenesis-related terms:\n")
  for (term in names(matches)) {
    status <- ifelse(matches[term], "✓ FOUND", "✗ NOT FOUND")
    cat(sprintf("  %s: %s\n", term, status))
  }
  
  # Count significant terms
  sig_count <- sum(flvcr2_enrich$GO_BP$p.adjust < 0.05)
  cat(sprintf("\nTotal significant GO BP terms (q<0.05): %d\n", sig_count))
}

cat("\n=== Enrichment Analysis Complete ===\n")
cat("Results saved to:", dir_out, "\n")
cat("Figures saved to:", dir_figs, "\n")
cat("\nNext step: Run 07_bulk_validation_composition_adjust.R\n")

