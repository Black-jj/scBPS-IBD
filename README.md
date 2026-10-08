# scBPS-IBD

Analysis code for a study that maps microbiota-associated host genetic signals onto single-cell, bulk and spatial transcriptomes of the human colon in inflammatory bowel disease (IBD), using the single-cell Bacteria Polygenic Score (scBPS).

## Repository layout

```
R/utils.R            shared functions, loaded by every script
scripts/             analysis steps, run in numerical order
python/              scBPS scoring and Geneformer in silico deletion
python/environments/ conda environment files for the two Python steps
```

## Requirements

**R (4.3 or later)** with the following packages:

- CRAN: Seurat, harmony, dplyr, tidyr, tidyverse, data.table, Matrix, ggplot2, ggpubr, ggrepel, ggraph, ggrastr, ggsci, ggplotify, ggstatsplot, cowplot, patchwork, aplot, corrplot, RColorBrewer, scales, reshape2, stringr, igraph, tidygraph, openxlsx, readxl, yaml, httr, fields, future, furrr, R.utils, KernSmooth, sctransform, msigdbr, magrittr
- Bioconductor: edgeR, limma, clusterProfiler, enrichplot, GOplot, org.Hs.eg.db, org.Mm.eg.db, AnnotationDbi, GSEABase, GSVA, AUCell, RcisTarget, SingleR, DelayedArray, monocle, mistyR
- GitHub or other sources: CellChat, CytoTRACE, ClusterGVis, jjAnno, scRNAtoolVis, IOBR, DoubletFinder, scBPS

**Python** environments for the scBPS and Geneformer steps:

```bash
conda env create -f python/environments/scbps.yml
conda env create -f python/environments/geneformer.yml
```

The R scripts call these environments through `conda run -n scBPS` and `conda run -n geneformer`. To use a specific interpreter instead, set `SCBPS_PYTHON` or `GENEFORMER_PYTHON`. Set `CONDA_EXE` if `conda` is not on the `PATH`.

## Data

Raw data are not included in this repository. Download the following GEO series and place the files under `data/`:

| Path | Source | Content |
|---|---|---|
| `data/GSE214695/` | [GSE214695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214695) | 10x matrix, barcode and feature files for 18 colon samples |
| `data/GSE235236/GSE235236_RAW/` | [GSE235236](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235236) | RSEM gene results and the FPKM expression matrix (56 samples) |
| `data/GSE235236/GSE235236_sample_manifest.tsv` | GSE235236 sample annotations | Tab-separated table with columns `gsm` and `group3` (HC, UC or CD) |
| `data/GSE234713/` | [GSE234713](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234713) | CosMx normalized matrix, annotation file and per-sample metadata files |

Reference resources used by the scripts (placed under `data/`): scBPS gene z-scores for the 207 Dutch Microbiome Project traits, gutMGene v2.0 human microbe-host gene relations, miRWalk exports, DGIdb interactions, RcisTarget hg19 CIS-BP motif files, MSigDB Hallmark v7.5.1 gene sets, a human KEGG pathway table, immune signature gene sets for ssGSEA, and the Geneformer V2-104M model with its gene dictionaries ([Hugging Face](https://huggingface.co/ctheodoris/Geneformer)).

## Running the analysis

Run every script from the repository root, in numerical order:

```bash
for f in scripts/*.R; do Rscript "$f"; done
```

Intermediate objects are written to `work/`, figures to `results/` and session information to `logs/`.

| Script | Analysis |
|---|---|
| `00_preflight.R` | Input checks and sample manifests |
| `01_scRNA_qc_integration.R` | Quality control, normalization and Harmony integration |
| `02_scRNA_annotation.R` | Cell annotation and composition |
| `03_marker_modules.R` | Marker modules and functional terms |
| `04_scBPS_key_cell.R` | scBPS cell-type prioritization |
| `05_key_cell_DE_GSEA.R` | Enterocyte pseudobulk differential expression and GSEA |
| `06_BPS_pathway_association.R` | BPS gradient and pathway activity |
| `07_trajectory.R` | CytoTRACE and Monocle2 trajectory |
| `08_cellchat.R` | Cell-cell communication |
| `09_gutMGene_network.R` | Microbe-gene network and key genes |
| `10_key_gene_enrichment.R` | Bulk GSEA and GSVA by key-gene expression |
| `11_immune_microenvironment.R` | ssGSEA immune signatures |
| `12_immune_regulators.R` | Key gene and immunomodulator correlations |
| `13_RcisTarget.R` | Transcription factor network |
| `14_miRNA_network.R` | miRNA network |
| `15_key_gene_scRNA_expression.R` | Key-gene expression in the single-cell atlas |
| `16_key_gene_pathway_activity.R` | AUCell pathway activity |
| `17_CosMx_spatial_mapping.R` | CosMx tissue compartments and epithelial subtypes |
| `18_mistyR.R` | mistyR spatial modelling |
| `19_key_gene_spatial_expression.R` | Spatial expression of the key genes |
| `20_geneformer.R` | Geneformer in silico deletion |
| `21_geneformer_enrichment.R` | Enrichment of perturbed genes |
| `22_DGIdb.R` | Drug-gene interaction network |
| `23_final_figures.R` | Final figure layouts |

## License

This code is released under the MIT License. See [LICENSE](LICENSE).
