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
work_dir <- file.path(root, "work/19_key_gene_spatial_expression")
result_dir <- file.path(root, "results/20.Spatial key-gene expression")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
normalized_file <- file.path(data_dir, "GSE234713_CosMx_normalized_matrix.txt.gz")
required_genes <- c("CCL20", "CLDN4", "CXCL8")
key_gene_file <- file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv")
pipeline_key_genes <- unique(fread(key_gene_file)$KeyGene)
stopifnot(all(required_genes %chin% pipeline_key_genes))
panel <- setdiff(names(fread(normalized_file, skip = 5, nrows = 0)), c("patient", "cell_id", "fov"))
coverage <- data.table(KeyGene = required_genes, status = ifelse(required_genes %chin% panel, "measurable", "unavailable_not_in_980_gene_panel"))
fwrite(coverage, file.path(work_dir, "key_gene_panel_coverage.tsv"), sep = "\t")
measurable <- coverage[status == "measurable", KeyGene]
stopifnot(identical(measurable, required_genes))
expr <- fread(normalized_file, skip = 5, select = c("patient", "cell_id", "fov", measurable))
expr[, `:=`(patient, gsub("[[:space:]]+", "_", trimws(patient)))]
expr[, `:=`(id, paste(patient, fov, cell_id, sep = "_"))]
meta <- as.data.table(readRDS(file.path(root, "work/17_cosmx_spatial_mapping/CosMx_primary_spatial_metadata.rds")))
d <- merge(expr, meta[, .(id, sample, group3, fov, cell_ID, x_um, y_um, subset, SingleR2, PublishedCompartment)], by = "id", 
    all.x = TRUE, sort = FALSE)
d[, `:=`(group2, fifelse(group3 == "HC", "HC", "IBD"))]
stopifnot(nrow(d) == 459095, all(d$sample == expr$patient), all(is.finite(d$x_um)), all(is.finite(d$y_um)))
sample_order <- c("HC_a", "HC_b", "HC_c", "UC_a", "UC_b", "UC_c", "CD_a", "CD_b", "CD_c")
stopifnot(setequal(unique(d$sample), sample_order))
spatial_limits <- rbindlist(lapply(measurable, function(gene) {
    cap <- as.numeric(quantile(d[[gene]], 0.99, na.rm = TRUE))
    if (!is.finite(cap) || cap <= 0) 
        cap <- max(d[[gene]], na.rm = TRUE)
    data.table(gene = gene, lower = 0, upper = cap, quantile_cap = 0.99)
}))
fwrite(spatial_limits, file.path(work_dir, "key_gene_spatial_color_limits.tsv"), sep = "\t")
for (gene in measurable) {
    current_gene <- gene
    cap <- spatial_limits[match(current_gene, spatial_limits$gene), upper]
    gene_plots <- lapply(sample_order, function(s) {
        z <- copy(d[sample == s])
        z[, `:=`(plot_value, pmin(pmax(get(gene), 0), cap))]
        ggplot(z, aes(x = x_um, y = y_um, color = plot_value)) + ggrastr::geom_point_rast(size = 0.18, alpha = 0.9, raster.dpi = 300) + 
            scale_color_viridis_c(option = "B", limits = c(0, cap), oob = scales::squish, name = paste0(gene, "\nnormalized")) + 
            scale_y_reverse() + coord_fixed() + theme_void() + labs(title = s) + theme(plot.title = element_text(hjust = 0.5, 
            size = 10), legend.position = "right")
    })
    pdf(file.path(result_dir, paste0(match(gene, measurable), ".", gene, " spatial expression.pdf")), width = 21, height = 17, 
        useDingbats = FALSE)
    print(ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = gene_plots, ncol = 3, nrow = 3, common.legend = TRUE, legend = "right"), 
        top = ggpubr::text_grob(paste0(gene, " spatial expression"), face = "bold", size = 16)))
    dev.off()
}
epi <- d[PublishedCompartment == "epi"]
stopifnot(nrow(epi) > 0, setequal(unique(epi$sample), sample_order))
summary_dt <- rbindlist(lapply(measurable, function(gene) {
    epi[, .(mean_expression = mean(get(gene), na.rm = TRUE), median_expression = median(get(gene), na.rm = TRUE), pct_expressing = mean(get(gene) > 
        0, na.rm = TRUE) * 100, cells = .N), by = .(sample, group3, group2)][, `:=`(gene, gene)]
}))
summary_dt[, `:=`(group2, factor(group2, levels = c("HC", "IBD")))]
summary_dt[, `:=`(group3, factor(group3, levels = c("HC", "UC", "CD")))]
fwrite(summary_dt, file.path(work_dir, "key_gene_published_epi_sample_expression.tsv"), sep = "\t")
tests <- rbindlist(lapply(measurable, function(gene) {
    current_gene <- gene
    x <- summary_dt[which(summary_dt$gene == current_gene)]
    values <- x$mean_expression
    has_ties <- anyDuplicated(values) > 0
    wt <- suppressWarnings(wilcox.test(mean_expression ~ group2, data = x, exact = !has_ties, conf.int = FALSE))
    kt <- kruskal.test(mean_expression ~ group3, data = x)
    data.table(gene = gene, HC_samples = x[group2 == "HC", .N], IBD_samples = x[group2 == "IBD", .N], UC_samples = x[group3 == 
        "UC", .N], CD_samples = x[group3 == "CD", .N], wilcoxon_p = wt$p.value, wilcoxon_exact = !has_ties, kruskal_p = kt$p.value)
}))
tests[, `:=`(wilcoxon_FDR, p.adjust(wilcoxon_p, "BH"))]
tests[, `:=`(kruskal_FDR, p.adjust(kruskal_p, "BH"))]
fwrite(tests, file.path(work_dir, "key_gene_published_epi_sample_level_tests.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "logs/19_key_gene_spatial_expression_sessionInfo.txt"))
