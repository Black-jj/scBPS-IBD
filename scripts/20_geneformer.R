source(file.path("R", "utils.R"))
geneformer_input_dir <- file.path(client_root, "work", "20_geneformer", "gf_input")
geneformer_result_dir <- file.path(client_root, "results", "20_Geneformer")
geneformer_enrichment_dir <- file.path(client_root, "results", "21_Geneformer_Downstream_Enrichment")
dir.create(geneformer_input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geneformer_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geneformer_enrichment_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(geneformer_result_dir, full.names = TRUE), recursive = TRUE)
unlink(list.files(geneformer_enrichment_dir, full.names = TRUE), recursive = TRUE)
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(CellChat)
})
geneformer_object <- readRDS(file.path(client_root, "work", "02_scRNA_annotation_abundance", "scRNA_annotated.rds"))
geneformer_cells <- rownames(geneformer_object@meta.data)[geneformer_object$cellType_1 == "Enterocytes"]
stopifnot(length(geneformer_cells) == 4037L)
geneformer_enterocytes <- subset(geneformer_object, cells = geneformer_cells)
geneformer_counts <- GetAssayData(geneformer_enterocytes, assay = "RNA", slot = "counts")
Matrix::writeMM(geneformer_counts, file.path(geneformer_input_dir, "matrix.mtx"))
writeLines(rownames(geneformer_counts), file.path(geneformer_input_dir, "genes.txt"))
writeLines(colnames(geneformer_counts), file.path(geneformer_input_dir, "barcodes.txt"))
geneformer_meta <- geneformer_enterocytes@meta.data[, c("orig.ident", "group3", "group", "cellType_1"), drop = FALSE]
geneformer_meta$barcode <- rownames(geneformer_meta)
write.csv(geneformer_meta, file.path(geneformer_input_dir, "meta.csv"), row.names = FALSE)
geneformer_interactions <- CellChatDB.human$interaction
geneformer_complexes <- CellChatDB.human$complex
geneformer_ligands <- unique(as.character(geneformer_interactions$ligand))
geneformer_flat_ligands <- character()
for (geneformer_ligand in geneformer_ligands) {
    if (!is.na(geneformer_ligand) && geneformer_ligand %in% rownames(geneformer_complexes)) {
        geneformer_subunits <- as.character(unlist(geneformer_complexes[geneformer_ligand, ]))
        geneformer_flat_ligands <- c(geneformer_flat_ligands, geneformer_subunits[geneformer_subunits != "" & !is.na(geneformer_subunits)])
    }
    else {
        geneformer_flat_ligands <- c(geneformer_flat_ligands, geneformer_ligand)
    }
}
geneformer_flat_ligands <- unique(trimws(geneformer_flat_ligands))
geneformer_flat_ligands <- geneformer_flat_ligands[geneformer_flat_ligands != "" & !is.na(geneformer_flat_ligands)]
writeLines(geneformer_flat_ligands, file.path(geneformer_input_dir, "ligands.txt"))
rm(geneformer_object, geneformer_enterocytes, geneformer_counts)
gc()
conda_executable <- Sys.getenv("CONDA_EXE", "conda")
geneformer_python <- Sys.getenv("GENEFORMER_PYTHON", "")
geneformer_script <- file.path(client_root, "python", "geneformer_pipeline.py")
geneformer_status <- if (nzchar(geneformer_python)) system2(geneformer_python, geneformer_script) else system2(conda_executable, 
    c("run", "--no-capture-output", "-n", "geneformer", "python", geneformer_script))
stopifnot(geneformer_status == 0L)
