source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrastr)
    library(ggpubr)
    library(yaml)
})
set.seed(20260804)
root <- client_root
data_dir <- file.path(root, "data/GSE234713")
work_dir <- file.path(root, "work/17_cosmx_spatial_mapping")
result_dir <- file.path(root, "results/18.CosMx")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
expr_files <- sort(list.files(data_dir, pattern = "^GSM[0-9]+_.+_exprMat_file[.]csv[.]gz$", full.names = TRUE))
meta_files <- sort(list.files(data_dir, pattern = "^GSM[0-9]+_.+_metadata_file[.]csv[.]gz$", full.names = TRUE))
stopifnot(length(expr_files) == 9, length(meta_files) == 9)
sample_from_file <- function(x) sub("_metadata_file[.]csv[.]gz$", "", sub("^GSM[0-9]+_", "", basename(x)))
expr_sample_from_file <- function(x) sub("_exprMat_file[.]csv[.]gz$", "", sub("^GSM[0-9]+_", "", basename(x)))
raw_ids <- rbindlist(lapply(expr_files, function(f) {
    sample <- expr_sample_from_file(f)
    x <- fread(f, select = c("fov", "cell_ID"))
    x <- x[cell_ID != 0]
    x[, `:=`(sample = sample, id = paste(sample, fov, cell_ID, sep = "_"))]
    x[, .(id, sample, fov, cell_ID)]
}))
coords <- rbindlist(lapply(meta_files, function(f) {
    sample <- sample_from_file(f)
    x <- fread(f, select = c("fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"))
    x[, `:=`(sample = sample, id = paste(sample, fov, cell_ID, sep = "_"), x_um = CenterX_global_px * 0.18, y_um = CenterY_global_px * 
        0.18)]
    x[, .(id, sample, fov, cell_ID, x_um, y_um)]
}))
stopifnot(nrow(raw_ids) == 482277, nrow(coords) == 482277, uniqueN(raw_ids$id) == nrow(raw_ids), uniqueN(coords$id) == nrow(coords), 
    setequal(raw_ids$id, coords$id))
normalized_file <- file.path(data_dir, "GSE234713_CosMx_normalized_matrix.txt.gz")
norm_ids <- fread(normalized_file, skip = 5, select = 1:3)
setnames(norm_ids, c("patient", "cell_id", "fov"))
norm_ids[, `:=`(patient, gsub("[[:space:]]+", "_", trimws(patient)))]
norm_ids[, `:=`(sample = patient, id = paste(patient, fov, cell_id, sep = "_"))]
annotation <- fread(file.path(data_dir, "GSE234713_CosMx_annotation.csv.gz"))
annotation_parts <- tstrsplit(annotation$id, "_", fixed = TRUE)
stopifnot(length(annotation_parts) == 4)
annotation[, `:=`(sample = paste(annotation_parts[[1]], annotation_parts[[2]], sep = "_"), cell_ID = as.integer(annotation_parts[[3]]), 
    fov = as.integer(annotation_parts[[4]]))]
annotation[, `:=`(join_id, paste(sample, fov, cell_ID, sep = "_"))]
annotation[, `:=`(PublishedCompartment, trimws(subset))]
annotation[, `:=`(PublishedCellType, trimws(SingleR2))]
annotation[is.na(PublishedCellType) | !nzchar(PublishedCellType), `:=`(PublishedCellType, fifelse(!is.na(PublishedCompartment) & 
    nzchar(PublishedCompartment), paste0(PublishedCompartment, " (unspecified)"), "Unresolved"))]
mapped <- merge(norm_ids[, .(id, sample, fov, cell_ID = cell_id)], coords[, .(id, x_um, y_um)], by = "id", all.x = TRUE, 
    sort = FALSE)
mapped <- merge(mapped, annotation[, .(id = join_id, published_annotation_id = id, subset, SingleR2, PublishedCompartment, 
    PublishedCellType)], by = "id", all.x = TRUE, sort = FALSE)
mapped[is.na(PublishedCompartment) | !nzchar(PublishedCompartment), `:=`(PublishedCompartment, "Unresolved")]
mapped[is.na(PublishedCellType) | !nzchar(PublishedCellType), `:=`(PublishedCellType, "Unresolved")]
mapped[, `:=`(CellType, PublishedCellType)]
mapped[, `:=`(group3, sub("_.*$", "", sample))]
published_compartments <- c("epi", "stroma", "tcells", "myeloids", "plasmas")
epithelial_subtypes <- sort(unique(mapped[PublishedCompartment == "epi" & !is.na(SingleR2) & nzchar(trimws(SingleR2)), trimws(SingleR2)]))
stopifnot(nrow(mapped) == 459095, all(is.finite(mapped$x_um)), all(is.finite(mapped$y_um)), sum(mapped$PublishedCompartment %chin% 
    published_compartments) == 458983, setequal(unique(mapped[PublishedCompartment != "Unresolved", PublishedCompartment]), 
    published_compartments), sum(mapped$PublishedCellType != "Unresolved") == 458983, length(epithelial_subtypes) == 9)
audit <- rbindlist(list(raw_ids[, .(cells = .N), by = sample][, `:=`(universe, "raw_integer_counts")], coords[, .(cells = .N), 
    by = sample][, `:=`(universe, "coordinates")], annotation[, .(cells = .N), by = sample][, `:=`(universe, "published_annotation")], 
    norm_ids[, .(cells = .N), by = sample][, `:=`(universe, "published_normalized_expression")], mapped[PublishedCompartment %chin% 
        published_compartments, .(cells = .N), by = sample][, `:=`(universe, "published_subset_annotated_expression")], mapped[PublishedCompartment == 
        "epi" & SingleR2 %chin% epithelial_subtypes, .(cells = .N), by = sample][, `:=`(universe, "published_epithelial_subtype_expression")]), 
    use.names = TRUE)
setcolorder(audit, c("sample", "universe", "cells"))
setorder(audit, sample, universe)
fwrite(audit, file.path(work_dir, "CosMx_mapping_audit.tsv"), sep = "\t")
unlink(file.path(work_dir, "CosMx_label_harmonization.tsv"))
fwrite(annotation[, .(published_subset = subset, published_label = SingleR2, PublishedCompartment, PublishedCellType)][, 
    unique(.SD)], file.path(work_dir, "CosMx_published_annotation_levels.tsv"), sep = "\t")
fwrite(mapped[, .(sample, group3, fov, cell_ID, id, published_annotation_id, x_um, y_um, published_subset = subset, published_label = SingleR2, 
    PublishedCompartment, PublishedCellType, CellType)], file.path(work_dir, "CosMx_primary_spatial_metadata.tsv.gz"), sep = "\t")
saveRDS(mapped, file.path(work_dir, "CosMx_primary_spatial_metadata.rds"), compress = FALSE)
sample_order <- c("HC_a", "HC_b", "HC_c", "UC_a", "UC_b", "UC_c", "CD_a", "CD_b", "CD_c")
stopifnot(setequal(unique(mapped$sample), sample_order))
compartment_palette <- c(epi = "#E64B35", stroma = "#4DBBD5", tcells = "#00A087", myeloids = "#3C5488", plasmas = "#F39B7F")
epithelial_palette <- setNames(grDevices::hcl.colors(length(epithelial_subtypes), palette = "Dynamic"), epithelial_subtypes)
compartment_plots <- lapply(sample_order, function(s) {
    d <- mapped[sample == s]
    colored <- d[PublishedCompartment %chin% published_compartments]
    ggplot(d, aes(x = x_um, y = y_um)) + ggrastr::geom_point_rast(color = "grey88", size = 0.22, alpha = 0.55, raster.dpi = 300) + 
        ggrastr::geom_point_rast(data = colored, aes(color = PublishedCompartment), size = 0.18, alpha = 0.9, raster.dpi = 300) + 
        scale_color_manual(values = compartment_palette, breaks = published_compartments, drop = FALSE) + scale_y_reverse() + 
        coord_fixed() + theme_void() + labs(title = s, color = "Published subset") + theme(plot.title = element_text(hjust = 0.5, 
        size = 11), legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 13)) + 
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 1))
})
pdf(file.path(result_dir, "1.Published compartment map.pdf"), width = 21, height = 17, useDingbats = FALSE)
print(ggpubr::ggarrange(plotlist = compartment_plots, ncol = 3, nrow = 3, common.legend = TRUE, legend = "top"))
dev.off()
epithelial_plots <- lapply(sample_order, function(s) {
    d <- mapped[sample == s]
    epi <- d[PublishedCompartment == "epi" & SingleR2 %chin% epithelial_subtypes]
    epi[, `:=`(EpithelialSubtype, trimws(SingleR2))]
    ggplot(d, aes(x = x_um, y = y_um)) + ggrastr::geom_point_rast(color = "grey88", size = 0.22, alpha = 0.5, raster.dpi = 300) + 
        ggrastr::geom_point_rast(data = epi, aes(color = EpithelialSubtype), size = 0.2, alpha = 0.95, raster.dpi = 300) + 
        scale_color_manual(values = epithelial_palette, breaks = epithelial_subtypes, drop = FALSE) + scale_y_reverse() + 
        coord_fixed() + theme_void() + labs(title = s, color = "Published epithelial subtype") + theme(plot.title = element_text(hjust = 0.5, 
        size = 11), legend.position = "top", legend.text = element_text(size = 10), legend.title = element_text(size = 12), 
        legend.key.width = grid::unit(0.42, "cm")) + guides(color = guide_legend(override.aes = list(size = 2.8, alpha = 1), 
        ncol = 5, byrow = TRUE))
})
pdf(file.path(result_dir, "2.Epithelial subtype map.pdf"), width = 21, height = 17, useDingbats = FALSE)
print(ggpubr::ggarrange(plotlist = epithelial_plots, ncol = 3, nrow = 3, common.legend = TRUE, legend = "top"))
dev.off()
writeLines(capture.output(sessionInfo()), file.path(root, "logs/17_cosmx_spatial_mapping_sessionInfo.txt"))
