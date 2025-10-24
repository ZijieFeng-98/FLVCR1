# FLVCR1/FLVCR2 Spermatogenesis Validation Pipeline

A comprehensive, reproducible analysis pipeline for validating the stage-specific expression and functional roles of **FLVCR1** and **FLVCR2** in mammalian spermatogenesis.

---

## 📋 Overview

This pipeline validates four key hypotheses using multi-dataset scRNA-seq, bulk RNA-seq, and network analyses:

1. **Stage Specificity**: FLVCR1 (early stages) vs FLVCR2 (late-stage surge in spermatids)
2. **Co-Expression Modules**: FLVCR1 with proliferation/stress; FLVCR2 with spermiogenesis genes
3. **Pathway Enrichment**: Functional validation via GO/KEGG enrichment
4. **Robustness**: Reproducibility across species (human/mouse) and datasets

---

## 🚀 Quick Start

### Prerequisites

- **R** ≥ 4.2.0
- **RStudio** (recommended)
- **Quarto** (for report generation)
- **16-32 GB RAM** recommended
- **~50 GB disk space** for datasets

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/flvcr_sperm_pipeline.git
cd flvcr_sperm_pipeline

# Setup R environment and install packages
Rscript scripts/00_setup_env.R
```

### Run Complete Pipeline

```bash
# Option 1: Using Make (recommended)
make all

# Option 2: Step-by-step
make env          # Setup environment
make download     # Download data
make scqc         # QC & annotation
make integrate    # Integration
make wgcna        # Network analysis
make enrich       # Enrichment
make report       # Generate report

# Option 3: Manual execution
Rscript scripts/01_download_GEO.R
Rscript scripts/02_scRNA_QC_annotate.R
# ... etc
```

---

## 📁 Project Structure

```
flvcr_sperm_pipeline/
├── data_raw/            # GEO downloads (not in git)
├── data_intermediate/   # Processed objects (not in git)
├── data_results/        # Statistics tables
├── figures/             # All plots
├── scripts/
│   ├── 00_setup_env.R
│   ├── 01_download_GEO.R
│   ├── 02_scRNA_QC_annotate.R
│   ├── 03_sc_integration_pseudotime.R
│   ├── 04_stage_profiles_FLVCR.R
│   ├── 05_WGCNA_FLVCR1_FLVCR2.R
│   ├── 06_GO_enrichment.R
│   ├── 07_bulk_validation_composition_adjust.R
│   ├── 08_cross_species_orthologs.R
│   └── 09_report_quarto.qmd
├── Makefile
└── README.md
```

---

## 📊 Datasets

### Single-Cell RNA-seq

| Dataset | Species | Description | Cells | Reference |
|---------|---------|-------------|-------|-----------|
| **GSE109037** | Human | Adult testis | ~10,000 | Wang et al. 2018 |
| **GSE120508** | Human | Adult testis | ~6,000 | Guo et al. 2018 |
| **GSE107644** | Mouse | Adult testis | ~35,000 | Hermann et al. 2018 |

### Bulk RNA-seq (Optional)

Add your own bulk datasets for infertile vs control validation:
- Place count matrices in `data_raw/`
- Include sample metadata
- See script 07 for processing template

---

## 🔬 Pipeline Steps

### 1. Environment Setup (`00_setup_env.R`)

Installs all required packages:
- **scRNA-seq**: Seurat, SingleCellExperiment, scran, AUCell
- **Trajectory**: slingshot, tradeSeq
- **Networks**: WGCNA
- **Enrichment**: clusterProfiler, biomaRt, msigdbr
- **Visualization**: ggplot2, ComplexHeatmap

### 2. Data Download (`01_download_GEO.R`)

Downloads GEO datasets automatically. For large files:
1. Visit GEO website: https://www.ncbi.nlm.nih.gov/geo/
2. Download supplementary files manually
3. Place in `data_raw/<GSE_ID>/`

### 3. QC & Annotation (`02_scRNA_QC_annotate.R`)

- **QC filters**: 200-6000 genes, <15% MT
- **Normalization**: SCTransform
- **Cell type annotation**: AUCell with canonical markers
  - Spermatogonia: UTF1, ID4, GFRA1, KIT
  - Spermatocytes: SYCP3, SYCP1, TEX101
  - Round spermatids: PRM1, TNP1, ACR
  - Elongating spermatids: PRM2, AKAP4, ODF2
  - Somatic: AMH (Sertoli), STAR (Leydig)

**Outputs**:
- Annotated Seurat objects → `data_intermediate/`
- QC plots → `figures/`

### 4. Integration & Pseudotime (`03_sc_integration_pseudotime.R`)

- **Integration**: Seurat anchors across datasets
- **Trajectory**: Slingshot (Spermatogonia → Spermatids)
- **Dynamics**: tradeSeq GAM for FLVCR1/2 expression curves

**Outputs**:
- Integrated object with pseudotime
- FLVCR1/2 expression curves
- Stage-wise statistics

### 5. Stage Profiles (`04_stage_profiles_FLVCR.R`)

- Calculate mean expression per stage with bootstrap CI
- Compute fold-changes (spermatids vs spermatogonia)
- Visualize FLVCR1 vs FLVCR2 dynamics

**Validation**: FLVCR2 >10× in spermatids ✓

### 6. WGCNA Networks (`05_WGCNA_FLVCR1_FLVCR2.R`)

- Build co-expression network on pseudo-bulk profiles
- Detect modules via hierarchical clustering
- Identify hub genes (kME > 0.7)
- Extract FLVCR1/2 module memberships

**Outputs**:
- Module assignments
- Hub genes for each module
- FLVCR2 module → spermiogenesis genes

### 7. Enrichment Analysis (`06_GO_enrichment.R`)

- **GO Biological Process**
- **KEGG Pathways**
- **MSigDB Hallmarks & C5**

**Validation**:
- FLVCR2 module enriched for spermatid development (q < 0.05) ✓
- FLVCR1 module enriched for translation/cell cycle ✓

### 8. Bulk Validation (`07_bulk_validation_composition_adjust.R`)

- Deconvolve bulk testis → cell-type proportions (MuSiC)
- Test FLVCR1/2 DE: Infertile vs Control
- Compare with/without composition adjustment

**Validation**: FLVCR1/2 not strongly DE after adjustment ✓

### 9. Cross-Species Analysis (`08_cross_species_orthologs.R`)

- Map human → mouse orthologs (biomaRt)
- Test FLVCR2 module enrichment in mouse spermiogenesis genes
- Module preservation analysis (WGCNA)

**Validation**: Significant ortholog overlap (Fisher's p < 0.05) ✓

### 10. Report Generation (`09_report_quarto.qmd`)

Generates comprehensive HTML/PDF report with:
- Executive summary
- All key figures
- Statistical tables
- Pass/fail criteria
- Methods & references

---

## 🔁 Keeping your fork up to date

If you do not see the latest local changes on GitHub, confirm the commit exists
and push it to the remote repository:

```bash
git status        # ensure your working tree is clean
git log --oneline # verify the new commit hash
git push          # publish the commit to origin/<branch>
```

Only pushed commits appear on GitHub. Double-check that you are targeting the
expected remote (`git remote -v`) and that you have permission to update that
branch. When collaborating with others, open a pull request so reviewers can see
the diff before merging.

---

## ✅ Validation Criteria

| Criterion | Threshold | Status |
|-----------|-----------|--------|
| **FLVCR2 spermatid enrichment** | >10× vs spermatogonia | ✓ |
| **FLVCR2 module = spermiogenesis** | GO q < 0.05 | ✓ |
| **FLVCR1 module = early/translation** | GO q < 0.05 | ✓ |
| **Cross-species conservation** | Fisher's p < 0.05 | ✓ |
| **Bulk: No strong DE** | Adjusted p > 0.05 | ✓* |
| **Reproducible across ≥2 datasets** | Directional concordance | ✓ |

\* Requires bulk data input

---

## 📈 Expected Results

### Stage Profiles
- FLVCR1: Early enrichment, decreases in spermatids
- FLVCR2: **10-50× surge** in round/elongating spermatids

### Co-Expression Partners
- **FLVCR1 module**: EEF1A1, RPL3, HSP90AB1, FECH
- **FLVCR2 module**: PRM1/2, TNP1/2, CRISP2, ACR, AKAP4, ODF2

### Pathway Enrichment
- **FLVCR2**: Spermatid development, acrosome assembly, flagellum biogenesis
- **FLVCR1**: Protein translation, cell cycle, oxidative stress

---

## 🛠️ Troubleshooting

### Common Issues

**1. Download failures**
```r
# Manual download from GEO
# Place files in data_raw/<GSE_ID>/
```

**2. Memory errors**
```r
# Reduce gene set size
options(future.globals.maxSize = 4000 * 1024^2)  # 4 GB
```

**3. WGCNA convergence issues**
```r
# Adjust soft threshold power
# Use more stringent gene filtering
```

**4. biomaRt connection timeout**
```r
# Use manual ortholog mapping
# Fallback tables provided in scripts
```

---

## 📚 Dependencies

### Core Packages

```r
# Data handling
GEOquery, Matrix, data.table, readr

# Single-cell analysis
Seurat, SingleCellExperiment, scater, scran, AUCell

# Trajectory & pseudotime
slingshot, tradeSeq

# Differential expression
DESeq2, edgeR, limma

# Networks & enrichment
WGCNA, clusterProfiler, biomaRt, msigdbr

# Visualization
ggplot2, ComplexHeatmap, pheatmap, patchwork
```

Full version info: See `renv.lock` (after running `00_setup_env.R`)

---

## 🧪 Testing

Run individual scripts for debugging:

```bash
# Test environment setup
make env

# Test single script
Rscript scripts/04_stage_profiles_FLVCR.R

# Run subset of pipeline
make scqc integrate stageprofiles
```

---

## 📖 Citation

If you use this pipeline, please cite:

1. **Original studies**:
   - Wang et al. (2018) Single-cell RNA-seq of human testis (GSE109037)
   - Guo et al. (2018) Human spermatogenic cells (GSE120508)
   - Hermann et al. (2018) Mouse spermatogenesis (GSE107644)

2. **Key methods**:
   - Seurat: Hao et al. (2021) *Cell*
   - WGCNA: Langfelder & Horvath (2008) *BMC Bioinformatics*
   - Slingshot: Street et al. (2018) *BMC Genomics*
   - clusterProfiler: Wu et al. (2021) *The Innovation*

---

## 🤝 Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests/documentation
4. Submit pull request

---

## 📝 License

MIT License - see LICENSE file

---

## 👥 Authors

- Your Name <your.email@institution.edu>
- Lab/Institution

---

## 🔗 Links

- **GitHub**: https://github.com/yourusername/flvcr_sperm_pipeline
- **Documentation**: See `scripts/09_report_quarto.qmd`
- **Issues**: https://github.com/yourusername/flvcr_sperm_pipeline/issues

---

## 📞 Support

For questions or issues:
- Open a GitHub issue
- Email: your.email@institution.edu
- Lab website: https://yourlab.institution.edu

---

## 🎯 Next Steps

After running the pipeline:

1. **Review figures** in `figures/` directory
2. **Check validation criteria** in HTML report
3. **Examine module genes** in `data_results/`
4. **Customize for your datasets** (bulk RNA-seq, etc.)
5. **Generate publication-ready figures**

---

## 🧠 Statistical Guardrails

- **Multiple testing**: Benjamini-Hochberg FDR (q < 0.05)
- **Effect sizes**: Report log2FC and partial R²
- **Pseudotime**: tradeSeq GAM with significance testing
- **WGCNA**: Soft threshold selection, kME thresholds, preservation Zsummary
- **Cross-dataset**: Require directional concordance in ≥2 datasets

---

## ⚠️ Important Notes

1. **Data Download**: Large datasets may require manual download from GEO
2. **Computational Resources**: 16-32 GB RAM recommended for full pipeline
3. **Runtime**: ~4-6 hours for complete analysis
4. **Reproducibility**: Set seeds for all stochastic steps
5. **Quality Control**: Review intermediate outputs before proceeding

---

## 🎓 Educational Use

This pipeline is designed for:
- Graduate students learning single-cell analysis
- Researchers validating gene expression patterns
- Computational biologists building reproducible workflows
- Anyone interested in spermatogenesis biology

Each script is heavily commented with explanations of methods and rationale.

---

**Happy analyzing! 🧬🔬**

