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

Reference resources are expected under `data/`:

| Path | Source |
|---|---|
| `data/reference/scBPS/gene_zscore_Dutch207.txt` | scBPS gene z-scores for the 207 Dutch Microbiome Project traits (Li et al., *Nat Microbiol* 2025) |
| `data/gutMGene_v2_20260802/microbe_gene.tsv` | Human gut microbe-host gene relations from [gutMGene v2.0](http://bio-computing.hrbmu.edu.cn/gutmgene) |
| `data/miRWalk_20260802/` | Human miRNA-target exports from [miRWalk](http://mirwalk.umm.uni-heidelberg.de/) for CCL20, CLDN4 and CXCL8 |
| `data/reference/interactions.tsv` | Interaction file from [DGIdb](https://www.dgidb.org/) |
| `data/reference/RcisTarget/` | hg19 CIS-BP 500 bp upstream motif ranking and motif annotation files for [RcisTarget](https://resources.aertslab.org/cistarget/) |
| `data/reference/h.all.v7.5.1.symbols.gmt` | MSigDB Hallmark gene sets, v7.5.1 |
| `data/reference/kegg.all.entrez.hsa.rds` | Human KEGG pathway table with columns `PathwayID`, `EntrezID` and `PathwayName` |
| `data/reference/ssGSEA/immune.gmt` | 29 immune signature gene sets used for ssGSEA |
| `data/reference/hallmark.gs.RData`, `data/reference/immune.gmt`, `data/reference/GenesetInfo.txt` | Gene sets used for AUCell and GSVA scoring |
| `data/reference/Immunomodulator_and_chemokines.txt` | Chemokine, receptor, MHC and immunomodulator gene list |
| `data/models/Geneformer-V2-104M/`, `data/reference/geneformer/` | Geneformer V2-104M model and gene dictionaries from [Hugging Face](https://huggingface.co/ctheodoris/Geneformer) |

## Running the analysis

Run every script from the repository root, in numerical order:

```bash
for f in scripts/*.R; do Rscript "$f"; done
```

Intermediate objects are written to `work/`, figures to `results/` and session information to `logs/`.

| Script | Analysis | Figure |
|---|---|---|
| `00_preflight.R` | Input checks and sample manifests | |
| `01_scRNA_qc_integration.R` | Quality control, normalization and Harmony integration | Supplementary Fig. 1 |
| `02_scRNA_annotation.R` | Cell annotation and composition | Fig. 1A–D |
| `03_marker_modules.R` | Marker modules and functional terms | Fig. 1E |
| `04_scBPS_key_cell.R` | scBPS cell-type prioritization | Fig. 2 |
| `05_key_cell_DE_GSEA.R` | Enterocyte pseudobulk differential expression and GSEA | Fig. 3A–D |
| `06_BPS_pathway_association.R` | BPS gradient and pathway activity | Fig. 3E |
| `07_trajectory.R` | CytoTRACE and Monocle2 trajectory | Supplementary Fig. 2 |
| `08_cellchat.R` | Cell-cell communication | Fig. 4A–E |
| `09_gutMGene_network.R` | Microbe-gene network and key genes | Fig. 4F |
| `10_key_gene_enrichment.R` | Bulk GSEA and GSVA by key-gene expression | Supplementary Fig. 3 |
| `11_immune_microenvironment.R` | ssGSEA immune signatures | Supplementary Fig. 4A–D |
| `12_immune_regulators.R` | Key gene and immunomodulator correlations | Supplementary Fig. 4E–I |
| `13_RcisTarget.R` | Transcription factor network | Fig. 5E |
| `14_miRNA_network.R` | miRNA network | Fig. 5F |
| `15_key_gene_scRNA_expression.R` | Key-gene expression in the single-cell atlas | Fig. 5A–C |
| `16_key_gene_pathway_activity.R` | AUCell pathway activity | Fig. 5D |
| `17_CosMx_spatial_mapping.R` | CosMx tissue compartments and epithelial subtypes | Fig. 6A–B |
| `18_mistyR.R` | mistyR spatial modelling | Fig. 6C–D |
| `19_key_gene_spatial_expression.R` | Spatial expression of the key genes | Fig. 6E–G |
| `20_geneformer.R` | Geneformer in silico deletion | Fig. 7A, D, G |
| `21_geneformer_enrichment.R` | Enrichment of perturbed genes | Fig. 7B, C, E, F, H, I |
| `22_DGIdb.R` | Drug-gene interaction network | Fig. 7J |
| `23_final_figures.R` | Final figure layouts | |

## License

This code is released under the MIT License. See [LICENSE](LICENSE).
