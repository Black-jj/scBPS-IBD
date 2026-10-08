source(file.path("R", "utils.R"))
project_dir <- client_root
work_dir <- file.path(project_dir, "work", "01_scRNA_qc_integration")
result_dir <- file.path(project_dir, "results", "01_IBD_scRNA_Atlas")
log_dir <- file.path(project_dir, "logs")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
manifest <- data.table::fread(file.path(project_dir, "work", "00_preflight", "GSE214695_sample_manifest.tsv"))
tenx_root <- file.path(project_dir, "work", "00_preflight", "GSE214695_10x")
raw_rds <- file.path(work_dir, "scRNA_raw_min500.rds")
if (file.exists(raw_rds)) {
    scRNA <- readRDS(raw_rds)
} else {
    objects <- lapply(seq_len(nrow(manifest)), function(i) {
        x <- Seurat::Read10X(file.path(tenx_root, manifest$sample[i]))
        if (is.list(x)) 
            x <- x[["Gene Expression"]]
        z <- Seurat::CreateSeuratObject(x, project = manifest$sample[i], min.cells = 10, min.features = 500)
        z$orig.ident <- manifest$sample[i]
        z$group3 <- manifest$group3[i]
        z$group <- ifelse(manifest$group3[i] == "HC", "Control", "Disease")
        z
    })
    names(objects) <- manifest$sample
    scRNA <- merge(objects[[1]], y = objects[-1], add.cell.ids = manifest$sample, merge.data = FALSE)
    saveRDS(scRNA, raw_rds, compress = FALSE)
}
scan <- function(...) c(1, 2)
options(future.globals.maxSize = 900 * 1024^3)
set.seed(20260802)
scRNA <- scRNAAutoAnno(Path = work_dir, SeuratObject = scRNA, Multi = TRUE, min.features = 500, Idents = "RNA_snn_res.0.2", 
    DoubletFinder = FALSE, species = "human", set.resolutions = seq(0.1, 0.6, by = 0.1), PC = 20, nfeatures = 2000, npcs = 50)
saveRDS(scRNA, file.path(work_dir, "scRNA_harmony.rds"), compress = FALSE)
pdfs <- list.files(work_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
stopifnot(length(pdfs) > 0)
for (f in pdfs) {
    target <- file.path(result_dir, paste0(basename(dirname(f)), "_", basename(f)))
    stopifnot(file.copy(f, target, overwrite = TRUE))
}
stopifnot(all(tolower(tools::file_ext(list.files(result_dir))) == "pdf"))
capture.output(sessionInfo(), file = file.path(log_dir, "01_scRNA_qc_integration_sessionInfo.txt"))
