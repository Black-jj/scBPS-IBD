source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(mistyR)
    library(furrr)
    library(R.utils)
    library(ggpubr)
    library(yaml)
})
set.seed(20260804)
root <- client_root
work_dir <- file.path(root, "work/18_mistyr_spatial_interactions")
result_dir <- file.path(root, "results/19.mistyR")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
meta <- as.data.table(readRDS(file.path(root, "work/17_cosmx_spatial_mapping/CosMx_primary_spatial_metadata.rds")))
meta <- meta[!is.na(PublishedCompartment) & nzchar(trimws(PublishedCompartment)) & PublishedCompartment != "Unresolved"]
meta[, `:=`(PublishedSubset, trimws(PublishedCompartment))]
published_features <- c("epi", "stroma", "tcells", "myeloids", "plasmas")
model_features <- make.names(published_features, unique = TRUE)
feature_lookup <- setNames(published_features, model_features)
stopifnot(nrow(meta) == 458983, setequal(unique(meta$PublishedSubset), published_features))
bin_size_um <- 20
juxta_um <- 20
para_um <- 100
para_neighbors <- 50
future::plan(future::multisession, workers = 8)
meta[, `:=`(bin_x = floor(x_um/bin_size_um), bin_y = floor(y_um/bin_size_um))]
bin_counts_long <- meta[, .(cells = .N), by = .(sample, group3, fov, bin_x, bin_y, PublishedSubset)]
bins <- dcast(bin_counts_long, sample + group3 + fov + bin_x + bin_y ~ PublishedSubset, value.var = "cells", fill = 0)
for (nm in setdiff(published_features, names(bins))) bins[, `:=`((nm), 0)]
bins[, `:=`(bin_id = paste(sample, fov, bin_x, bin_y, sep = "_"), x_um = (bin_x + 0.5) * bin_size_um, y_um = (bin_y + 0.5) * 
    bin_size_um, total_cells = rowSums(.SD)), .SDcols = published_features]
stopifnot(sum(bins$total_cells) == nrow(meta), uniqueN(bins$bin_id) == nrow(bins), all(bins$total_cells > 0))
make_bin_matrix <- function(d) {
    x <- log1p(as.matrix(d[, ..published_features]))
    rownames(x) <- d$bin_id
    colnames(x) <- model_features
    as.data.frame(x, check.names = FALSE)
}
build_sample_views <- function(sample_id) {
    d <- bins[sample == sample_id]
    view_file <- file.path(work_dir, paste0(sample_id, "_published_subset_20um_bin_views_fixedrowbind.rds"))
    if (file.exists(view_file)) 
        return(readRDS(view_file))
    x <- make_bin_matrix(d)
    fov_order <- sort(unique(d$fov))
    juxta_list <- vector("list", length(fov_order))
    para_list <- vector("list", length(fov_order))
    names(juxta_list) <- names(para_list) <- fov_order
    for (f in fov_order) {
        df <- d[fov == f]
        stopifnot(nrow(df) > 1)
        xf <- x[df$bin_id, , drop = FALSE]
        geometry <- data.frame(x = df$x_um, y = df$y_um, row.names = df$bin_id)
        initial <- mistyR::create_initial_view(xf)
        jv <- mistyR::add_juxtaview(initial, positions = geometry, neighbor.thr = juxta_um + 1e-06, cached = FALSE, verbose = FALSE)[[paste0("juxtaview.", 
            juxta_um + 1e-06)]]$data
        pv <- mistyR::add_paraview(initial, positions = geometry, l = para_um, nn = min(para_neighbors, nrow(df) - 1), cached = FALSE, 
            verbose = FALSE)[[paste0("paraview.", para_um)]]$data
        rownames(jv) <- rownames(pv) <- df$bin_id
        juxta_list[[as.character(f)]] <- as.data.frame(jv, check.names = FALSE)
        para_list[[as.character(f)]] <- as.data.frame(pv, check.names = FALSE)
    }
    juxta <- do.call(rbind, unname(juxta_list))
    para <- do.call(rbind, unname(para_list))
    juxta <- juxta[rownames(x), , drop = FALSE]
    para <- para[rownames(x), , drop = FALSE]
    stopifnot(!anyNA(juxta), !anyNA(para), nrow(juxta) == nrow(x), nrow(para) == nrow(x))
    out <- list(intra = x, juxta = juxta, para = para)
    saveRDS(out, view_file, compress = FALSE)
    out
}
sample_order <- c("HC_a", "HC_b", "HC_c", "UC_a", "UC_b", "UC_c", "CD_a", "CD_b", "CD_c")
group_order <- c("HC", "IBD")
sample_manifest <- unique(meta[, .(sample, group3)])
sample_manifest[, `:=`(group2, fifelse(group3 == "HC", "HC", "IBD"))]
sample_manifest[, `:=`(sample, factor(sample, levels = sample_order))]
setorder(sample_manifest, sample)
sample_manifest[, `:=`(sample, as.character(sample))]
stopifnot(identical(sample_manifest$sample, sample_order))
fwrite(sample_manifest, file.path(work_dir, "CosMx_published_subset_sample_manifest.tsv"), sep = "\t")
label_counts <- meta[, .(cells = .N), by = .(sample, group3, PublishedSubset)]
fwrite(label_counts, file.path(work_dir, "CosMx_published_subset_counts.tsv"), sep = "\t")
fwrite(bins[, c("sample", "group3", "fov", "bin_x", "bin_y", "bin_id", "x_um", "y_um", "total_cells", published_features), 
    with = FALSE], file.path(work_dir, "CosMx_published_subset_20um_bins.tsv.gz"), sep = "\t")
fwrite(data.table(PublishedSubset = published_features, ModelFeature = model_features), file.path(work_dir, "CosMx_published_subset_model_label_map.tsv"), 
    sep = "\t")
model_dirs <- setNames(character(nrow(sample_manifest)), sample_manifest$sample)
for (s in sample_manifest$sample) {
    views_data <- build_sample_views(s)
    out_alias <- file.path(work_dir, "models_published_subset_20um_bins_fixedrowbind", s)
    model_dirs[s] <- out_alias
    complete_model <- file.exists(file.path(out_alias, "performance.txt")) && file.exists(file.path(out_alias, "coefficients.txt")) && 
        length(list.files(out_alias, pattern = "^importances_.+_intra[.]txt$")) == length(model_features)
    if (!complete_model) {
        views <- mistyR::create_initial_view(views_data$intra)
        views <- mistyR::add_views(views, c(mistyR::create_view(paste0("juxtaview.", juxta_um), views_data$juxta, paste0("juxta.", 
            juxta_um)), mistyR::create_view(paste0("paraview.", para_um), views_data$para, paste0("para.", para_um))))
        mistyR::run_misty(views, results.folder = out_alias, seed = 20260804, cv.folds = 10, cached = FALSE)
    }
    rm(views_data)
    gc()
}
future::plan(future::multisession, workers = 3)
misty_res_slide <- mistyR::collect_results(unname(model_dirs[sample_order]))
for (nm in intersect(c("improvements", "contributions", "improvements.stats", "contributions.stats"), names(misty_res_slide))) {
    if ("target" %in% names(misty_res_slide[[nm]])) 
        misty_res_slide[[nm]]$target <- unname(feature_lookup[misty_res_slide[[nm]]$target])
}
for (nm in intersect(c("importances", "importances.aggregated"), names(misty_res_slide))) {
    if ("Predictor" %in% names(misty_res_slide[[nm]])) 
        misty_res_slide[[nm]]$Predictor <- unname(feature_lookup[misty_res_slide[[nm]]$Predictor])
    if ("Target" %in% names(misty_res_slide[[nm]])) 
        misty_res_slide[[nm]]$Target <- unname(feature_lookup[misty_res_slide[[nm]]$Target])
}
saveRDS(misty_res_slide, file.path(work_dir, "misty_res_slide.rds"))
pdf(file.path(result_dir, "1.Intra_heatmap.pdf"), width = 10, height = 8, useDingbats = FALSE)
mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
dev.off()
pdf(file.path(result_dir, "2.Intra_communities.pdf"), width = 10, height = 8, useDingbats = FALSE)
mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
dev.off()
fwrite(as.data.table(misty_res_slide$contributions.stats), file.path(work_dir, "all_samples_published_subset_misty_view_contributions.tsv"), 
    sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "logs/18_mistyr_spatial_interactions_sessionInfo.txt"))
