source(file.path("R", "utils.R"))
project_dir <- client_root
data_dir <- file.path(project_dir, "data")
work_dir <- file.path(project_dir, "work", "00_preflight")
log_dir <- file.path(project_dir, "logs")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
sc_dir <- file.path(data_dir, "GSE214695")
mat_files <- sort(list.files(sc_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
bar_files <- sort(list.files(sc_dir, pattern = "_barcodes\\.tsv\\.gz$", full.names = TRUE))
feat_files <- sort(list.files(sc_dir, pattern = "_features\\.tsv\\.gz$", full.names = TRUE))
stopifnot(length(mat_files) == 18, length(bar_files) == 18, length(feat_files) == 18)
samples <- sub("_matrix\\.mtx\\.gz$", "", basename(mat_files))
stopifnot(all(vapply(samples, function(z) sum(grepl(z, basename(c(mat_files, bar_files, feat_files)), fixed = TRUE)) == 3, 
    logical(1))))
read_dim <- function(f) {
    con <- gzfile(f, "rt")
    on.exit(close(con))
    repeat {
        z <- readLines(con, n = 1)
        if (!startsWith(z, "%")) 
            return(scan(text = z, quiet = TRUE))
    }
}
dims <- t(vapply(mat_files, read_dim, numeric(3)))
sc_manifest <- data.frame(sample = samples, gsm = sub("_.*$", "", samples), group3 = sub("^.*_(HC|UC|CD)-.*$", "\\1", samples), 
    genes = dims[, 1], raw_barcodes = dims[, 2], nonzero_entries = dims[, 3], matrix = mat_files, stringsAsFactors = FALSE)
stopifnot(identical(as.integer(table(factor(sc_manifest$group3, levels = c("HC", "UC", "CD")))), c(6L, 6L, 6L)))
data.table::fwrite(sc_manifest, file.path(work_dir, "GSE214695_sample_manifest.tsv"), sep = "\t")
tenx_root <- file.path(work_dir, "GSE214695_10x")
dir.create(tenx_root, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(sc_manifest))) {
    d <- file.path(tenx_root, sc_manifest$sample[i])
    dir.create(d, showWarnings = FALSE)
    src <- c(bar_files[grepl(sc_manifest$sample[i], bar_files, fixed = TRUE)], feat_files[grepl(sc_manifest$sample[i], feat_files, 
        fixed = TRUE)], mat_files[i])
    dst <- file.path(d, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
    for (j in seq_along(dst)) if (!file.exists(dst[j])) 
        stopifnot(file.symlink(src[j], dst[j]))
}
bulk_dir <- file.path(data_dir, "GSE235236", "GSE235236_RAW")
rsem_files <- sort(list.files(bulk_dir, pattern = "RSEM\\.genes\\.results$", full.names = TRUE))
stopifnot(length(rsem_files) == 56)
xlsx <- file.path(bulk_dir, "GSE235236_FPKM_expression_matrix.xlsx")
stopifnot(file.exists(xlsx))
bulk_cols <- names(readxl::read_excel(xlsx, sheet = "FPKM_matrix", n_max = 1))[-1]
bulk_gsm <- sub("_.*$", "", bulk_cols)
bulk_manifest <- data.table::fread(file.path(data_dir, "GSE235236", "GSE235236_sample_manifest.tsv"))
bulk_manifest <- bulk_manifest[match(bulk_gsm, bulk_manifest$gsm), ]
bulk_manifest$expression_column <- bulk_cols
stopifnot(!anyNA(bulk_manifest$gsm), !anyNA(bulk_manifest$group3), identical(as.integer(table(factor(bulk_manifest$group3, 
    levels = c("HC", "UC", "CD")))), c(8L, 26L, 22L)))
data.table::fwrite(bulk_manifest, file.path(work_dir, "GSE235236_sample_manifest.tsv"), sep = "\t")
sp_dir <- file.path(data_dir, "GSE234713")
expr_file <- file.path(sp_dir, "GSE234713_CosMx_normalized_matrix.txt.gz")
ann_file <- file.path(sp_dir, "GSE234713_CosMx_annotation.csv.gz")
coord_files <- sort(list.files(sp_dir, pattern = "^GSM.*_metadata_file\\.csv\\.gz$", full.names = TRUE))
stopifnot(file.exists(expr_file), file.exists(ann_file), length(coord_files) == 9)
expr_id <- data.table::fread(expr_file, skip = 5, select = 1:3, showProgress = FALSE)
expr_id[, `:=`(id, paste(gsub(" ", "_", patient), cell_id, fov, sep = "_"))]
ann <- data.table::fread(ann_file, select = c("id", "subset", "SingleR2"), showProgress = FALSE)
coords <- data.table::rbindlist(lapply(coord_files, function(f) {
    d <- data.table::fread(f, select = c("fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"), showProgress = FALSE)
    sample <- sub("_metadata_file.csv.gz$", "", sub("^GSM[0-9]+_", "", basename(f)))
    d[, `:=`(patient, gsub("_", " ", sample))]
    d[, `:=`(id, paste(gsub(" ", "_", patient), cell_ID, fov, sep = "_"))]
    d
}))
sp_audit <- data.frame(metric = c("expression_cells", "annotation_cells", "coordinate_cells", "expression_with_annotation", 
    "expression_without_annotation", "expression_with_coordinates", "annotation_without_expression", "coordinates_without_expression"), 
    value = c(nrow(expr_id), nrow(ann), nrow(coords), sum(expr_id$id %in% ann$id), sum(!expr_id$id %in% ann$id), sum(expr_id$id %in% 
        coords$id), sum(!ann$id %in% expr_id$id), sum(!coords$id %in% expr_id$id)))
stopifnot(sp_audit$value[sp_audit$metric == "expression_cells"] == 459095, sp_audit$value[sp_audit$metric == "expression_with_coordinates"] == 
    459095)
data.table::fwrite(sp_audit, file.path(work_dir, "GSE234713_join_audit.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(log_dir, "00_preflight_sessionInfo.txt"))
