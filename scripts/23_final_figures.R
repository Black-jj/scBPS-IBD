source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout")))
suppressPackageStartupMessages({
    library(data.table)
    library(Seurat)
    library(ggplot2)
    library(ggpubr)
    library(ggrastr)
    library(igraph)
    library(dplyr)
})
set.seed(20260802)
root <- client_root
revision_work <- file.path(root, "work/23_final_figures")
result_dir <- file.path(root, "results/final_figures")
dir.create(revision_work, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- function(folder, file) {
    dir.create(file.path(result_dir, folder), recursive = TRUE, showWarnings = FALSE)
    file.path(result_dir, folder, file)
}
valid_pdf <- function(path) {
    info <- tryCatch(system2("pdfinfo", shQuote(path), stdout = TRUE, stderr = TRUE), error = function(e) character())
    page_line <- grep("^Pages:", info, value = TRUE)
    length(page_line) == 1L && as.integer(sub("^Pages:[[:space:]]*", "", page_line)) > 0L
}
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
stopifnot(all(qc_features %in% colnames(scRNA@meta.data)))
p <- VlnPlot(scRNA, features = qc_features, group.by = "orig.ident", ncol = 4)
p <- p & theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.5))
ggsave(pdf_path("1.QCFilter", "1.VlnPlot.pdf"), p, width = 24, height = 5.2, limitsize = FALSE)
p <- DimPlot(scRNA, group.by = "RNA_snn_res.0.1", label = TRUE, label.size = 3, reduction = "umap", raster = FALSE) + NoLegend() + 
    ggtitle(NULL)
ggsave(pdf_path("3.CellAnnotate", "1.Umap.Cluster.pdf"), p, width = 5, height = 5)
ordered_markers <- c("IL7R", "CD3D", "CD3E", "MZB1", "JCHAIN", "IGKC", "KRT20", "PIGR", "FABP1", "COL1A1", "DCN", "LUM", 
    "LYZ", "TYMP", "CD14", "MS4A1", "CD79A", "CD79B", "MUC2", "TFF3", "SPINK4", "CPA3", "KIT", "TPSAB1", "VWF", "PECAM1", 
    "KDR", "MKI67", "TOP2A", "STMN1", "RGS5", "ACTA2", "MYH11", "POU2F3", "TRPM5", "IL17RB", "SOX10", "S100B", "PLP1", "MPZ")
ordered_markers <- intersect(ordered_markers, rownames(scRNA))
p <- DotPlot(scRNA, features = ordered_markers, group.by = "cellType_1") + scale_color_gradientn(colours = c("#3B4CC0", "#F7F7F7", 
    "#B40426"), name = "Average Expression") + scale_size_continuous(range = c(0, 6)) + theme_bw() + theme(panel.border = element_rect(colour = "black", 
    fill = NA), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "top", legend.key.height = grid::unit(0.3, 
    "cm"), legend.key.width = grid::unit(0.8, "cm"))
ggsave(pdf_path("3.CellAnnotate", "3.DotPlot.pdf"), p, width = 13, height = 7)
rm(scRNA)
gc()
enrichment_dir <- file.path(result_dir, "11.Enrichment")
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)
gsea_work <- file.path(revision_work, "key_gene_KEGG_GSEA")
dir.create(gsea_work, recursive = TRUE, showWarnings = FALSE)
kegg_snapshot <- file.path(client_root, "data", "reference", "kegg.all.entrez.hsa.rds")
stopifnot(file.exists(kegg_snapshot))
kegg_reference <- readRDS(kegg_snapshot)
kegg_term2gene <- unique(kegg_reference[, c("PathwayID", "EntrezID")])
kegg_term2name <- unique(kegg_reference[, c("PathwayID", "PathwayName")])
gseKEGG <- function(geneList, organism = "hsa", ...) {
    geneList <- geneList[!is.na(geneList) & !duplicated(names(geneList))]
    clusterProfiler::GSEA(geneList = geneList, TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name, pvalueCutoff = 1, 
        verbose = FALSE)
}
append_kegg_category <- function(x) {
    x@result$category <- ifelse(grepl("^(KEGG_)?hsa05", x@result$ID), "Human Diseases", "Non-disease")
    x
}
bulk_mat <- readRDS(file.path(root, "work/10_key_gene_enrichment/bulk_log2_FPKM_matrix.rds"))
key_genes <- c("CCL20", "CLDN4", "CXCL8")
stopifnot(all(key_genes %in% rownames(bulk_mat)))
gsea_display_used <- character()
gsea_display_selected <- list()
GSEA_module(mat = bulk_mat, outpath = gsea_work, species = "human", KeyGene = key_genes, isKEGG = TRUE)
selected_display <- rbindlist(gsea_display_selected, use.names = TRUE, fill = TRUE)
fwrite(selected_display, file.path(gsea_work, "selected_GSEA_display_pathways.tsv"), sep = "\t")
for (i in seq_along(key_genes)) {
    src <- file.path(gsea_work, "GSEA", paste0(i, ".", key_genes[i]), "1.GSEA_multi_pathways.pdf")
    stopifnot(file.exists(src), valid_pdf(src))
    stopifnot(file.copy(src, file.path(enrichment_dir, paste0(i, ".", key_genes[i], " KEGG GSEA.pdf")), overwrite = TRUE))
}
rm(bulk_mat)
gc()
reg_result <- fread(file.path(root, "work/12_immune_regulators/key_gene_immune_regulator_spearman.tsv"))
reg_ref <- fread(file.path(client_root, "data", "reference", "Immunomodulator_and_chemokines.txt"))
type_value <- unique(reg_ref[tolower(type) == "immunostimulator", type])
stopifnot(length(type_value) == 1L)
ids <- reg_ref[type == type_value, Id]
data2 <- as.data.frame(reg_result[group == "All" & immuneGene %in% ids])
data2$pstar <- ifelse(data2$FDR < 0.001, "***", ifelse(data2$FDR < 0.01, "**", ifelse(data2$FDR < 0.05, "*", "")))
data2 <- na.omit(data2)
p <- ggplot(data2, aes(immuneGene, gene)) + geom_tile(aes(fill = cor)) + geom_text(aes(label = pstar), color = "black", size = 4) + 
    scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c", limit = c(-1, 1), name = "*    FDR < 0.05\n\n**  FDR < 0.01\n\n*** FDR < 0.001\n\nSpearman") + 
    labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(size = 12, 
    color = "black"), axis.ticks.y = element_blank(), panel.background = element_blank()) + coord_fixed()
ggsave(pdf_path("13.Immune association", "5.Immunostimulator.pdf"), width = 15, height = 4, plot = p, limitsize = FALSE)
tf_edges <- fread(file.path(root, "work/13_rcistarget/TF_gene_edges.tsv"))
key_genes <- c("CCL20", "CLDN4", "CXCL8")
stopifnot(setequal(unique(tf_edges$enrichedGenes), key_genes))
g <- graph_from_data_frame(d = tf_edges[, .(TF_highConf, enrichedGenes)], directed = TRUE)
V(g)$type <- ifelse(V(g)$name %in% tf_edges$TF_highConf, "TF", ifelse(V(g)$name %in% key_genes, "KeyGene", "Gene"))
V(g)$color <- c(TF = "#4DBBD5", Gene = "grey70", KeyGene = "#E64B35")[V(g)$type]
V(g)$size <- c(TF = 14, Gene = 11, KeyGene = 24)[V(g)$type]
lay <- layout_with_kk(g, maxiter = 5000)
pdf(pdf_path("14.Transcription factor", "1.TF-gene network.pdf"), width = 12, height = 8, useDingbats = FALSE)
plot(g, layout = lay, vertex.label.cex = 0.85, vertex.label.color = "black", edge.color = "grey70", edge.arrow.size = 0.35, 
    main = "TF-Key Gene Regulatory Network")
dev.off()
scrna_plot <- fread(file.path(root, "work/15_key_gene_scrna_expression/key_cell_sample_average_expression.tsv"))
scrna_plot[, `:=`(contrast, factor(ifelse(group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD")))]
scrna_stats <- scrna_plot[, .(HC_n = sum(contrast == "HC"), IBD_n = sum(contrast == "IBD"), pvalue = wilcox.test(expression ~ 
    contrast, exact = TRUE)$p.value), by = gene]
scrna_stats[, `:=`(FDR, p.adjust(pvalue, "BH"))]
fwrite(scrna_stats, file.path(revision_work, "scRNA_Enterocytes_IBD_vs_HC_tests.tsv"), sep = "\t")
display_genes <- scrna_stats[FDR < 0.05, gene]
stopifnot(length(display_genes) > 0L)
scrna_display <- scrna_plot[gene %in% display_genes]
p <- ggplot(scrna_display, aes(x = contrast, y = expression, fill = contrast, color = contrast)) + geom_boxplot(alpha = 0.3, 
    outlier.shape = NA) + geom_jitter(width = 0.08, size = 2) + facet_wrap(~gene, scales = "free_y", nrow = 1) + scale_fill_manual(values = c(HC = "deepskyblue", 
    IBD = "hotpink")) + scale_color_manual(values = c(HC = "dodgerblue", IBD = "plum3")) + stat_compare_means(method = "wilcox.test", 
    label = "p.format") + theme_bw() + labs(x = "", y = "Sample-average expression in Enterocytes") + theme(legend.position = "none")
ggsave(pdf_path("16.scRNA key-gene expression", "4.Sample-level boxplot.pdf"), p, width = max(7, length(display_genes) * 
    3.3), height = 5)
mapped <- as.data.table(readRDS(file.path(root, "work/17_cosmx_spatial_mapping/CosMx_primary_spatial_metadata.rds")))
sample_order <- intersect(c("HC_a", "HC_b", "HC_c", "UC_a", "UC_b", "UC_c", "CD_a", "CD_b", "CD_c"), unique(mapped$sample))
stopifnot(length(sample_order) == 9L)
pal <- c(T_cells = "#4DBBD5", Plasma_cells = "#E64B35", Enterocytes = "#00A087", Fibroblasts = "#3C5488", Myeloid = "#F39B7F", 
    B_cells = "#8491B4", Goblet_cells = "#91D1C2", Mast_cells = "#DC0000", Endothelial = "#7E6148", Proliferating_Tcells = "#B09C85", 
    Smooth_muscle_cells = "#6A6599", Tuft_cells = "#ED6A5A", Glial_cells = "#5FAD56", Unresolved = "grey75")
pdf(pdf_path("18.CosMx", "1.Cell-type spatial map.pdf"), width = 12, height = 9, useDingbats = FALSE)
for (s in sample_order) {
    z <- mapped[sample == s]
    p <- ggplot(z, aes(x = x_um, y = y_um)) + ggrastr::geom_point_rast(color = "grey88", size = 0.62, alpha = 0.75, raster.dpi = 300) + 
        ggrastr::geom_point_rast(aes(color = CellType), size = 0.16, alpha = 0.88, raster.dpi = 300) + scale_color_manual(values = pal, 
        drop = FALSE) + scale_y_reverse() + coord_fixed() + theme_void() + labs(title = paste0(s, " - CosMx cell types with full-cell tissue footprint"), 
        color = "Cell type") + theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3, 
        alpha = 1), ncol = 1))
    print(p)
}
dev.off()
spatial_dir <- file.path(result_dir, "20.Spatial key-gene expression")
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
normalized_file <- file.path(root, "data/GSE234713/GSE234713_CosMx_normalized_matrix.txt.gz")
spatial_genes <- c("CCL20", "CLDN4", "CXCL8")
expr <- fread(normalized_file, skip = 5, select = c("patient", "cell_id", "fov", spatial_genes))
expr[, `:=`(patient, gsub("[[:space:]]+", "_", trimws(patient)))]
expr[, `:=`(id, paste(patient, fov, cell_id, sep = "_"))]
spatial <- merge(expr, mapped[, .(id, sample, group3, fov, cell_ID, x_um, y_um, CellType)], by = "id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(spatial) == 459095L)
plot_index <- 1L
for (gene in spatial_genes) {
    cap <- as.numeric(quantile(spatial[[gene]], 0.99, na.rm = TRUE))
    if (!is.finite(cap) || cap <= 0) 
        cap <- max(spatial[[gene]], na.rm = TRUE)
    for (s in sample_order) {
        plot_index <- plot_index + 1L
        z <- spatial[sample == s]
        z[, `:=`(plot_value, pmin(get(gene), cap))]
        p <- ggplot(z, aes(x = x_um, y = y_um)) + ggrastr::geom_point_rast(color = "grey90", size = 0.55, alpha = 0.7, raster.dpi = 300) + 
            ggrastr::geom_point_rast(aes(color = plot_value), size = 0.18, alpha = 0.9, raster.dpi = 300) + scale_color_viridis_c(option = "B", 
            name = paste0(gene, "\nnormalized")) + scale_y_reverse() + coord_fixed() + theme_void() + labs(title = paste0(gsub("_", 
            " ", s), " - ", gene)) + theme(plot.title = element_text(hjust = 0.5, size = 15), legend.position = "right")
        out <- file.path(spatial_dir, paste0(plot_index, ".", gene, " - ", gsub("_", " ", s), ".pdf"))
        ggsave(out, p, width = 12, height = 9, limitsize = FALSE)
    }
}
spatial_summary <- fread(file.path(root, "work/19_key_gene_spatial_expression/key_gene_sample_celltype_expression.tsv"))
spatial_key <- spatial_summary[CellType == "Enterocytes"]
spatial_key[, `:=`(contrast, factor(ifelse(group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD")))]
spatial_stats <- spatial_key[, .(HC_n = sum(contrast == "HC"), IBD_n = sum(contrast == "IBD"), pvalue = wilcox.test(mean_expression ~ 
    contrast, exact = TRUE)$p.value), by = gene]
spatial_stats[, `:=`(FDR, p.adjust(pvalue, "BH"))]
fwrite(spatial_stats, file.path(revision_work, "CosMx_Enterocytes_IBD_vs_HC_tests.tsv"), sep = "\t")
p <- ggplot(spatial_key, aes(x = contrast, y = mean_expression, fill = contrast, color = contrast)) + geom_boxplot(alpha = 0.3, 
    outlier.shape = NA) + geom_jitter(width = 0.08, size = 2) + facet_wrap(~gene, scales = "free_y", nrow = 1) + scale_fill_manual(values = c(HC = "deepskyblue", 
    IBD = "hotpink")) + scale_color_manual(values = c(HC = "dodgerblue", IBD = "plum3")) + stat_compare_means(method = "wilcox.test", 
    label = "p.format") + theme_bw() + labs(x = "", y = "Sample-average CosMx expression in Enterocytes") + theme(legend.position = "none")
ggsave(file.path(spatial_dir, "29.Key-cell sample level.pdf"), p, width = 10, height = 5)
rm(mapped, expr, spatial)
gc()
drug_edges <- fread(file.path(root, "work/21_dgidb_drug_prediction/DGIdb_key_gene_drug_edges.tsv"))
g <- graph_from_data_frame(d = drug_edges, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% drug_edges$Drug, "Drug", "Gene")
V(g)$color <- ifelse(V(g)$type == "Drug", "#4DBBD5", "#E64B35")
V(g)$size <- ifelse(V(g)$type == "Drug", 3.5, 9)
lay <- layout_with_kk(g, maxiter = 8000)
pdf(pdf_path("23.Drug prediction", "1.Drug-gene network.pdf"), width = 20, height = 10, useDingbats = FALSE)
plot(g, layout = lay, vertex.label.cex = 0.48, vertex.label.color = "black", edge.color = "grey78", main = "Drug Prediction")
dev.off()
key_genes <- c("CCL20", "CLDN4", "CXCL8")
immune_work <- file.path(root, "work/11_immune_microenvironment")
immune_result <- fread(file.path(immune_work, "ssGSEA_scores_by_sample.tsv"), check.names = FALSE)
stopifnot(names(immune_result)[1] == "sample")
manifest <- fread(file.path(root, "work/10_key_gene_enrichment/bulk_sample_manifest.tsv"))
immune_cols <- setdiff(names(immune_result), "sample")
immune_long <- melt(immune_result, id.vars = "sample", measure.vars = immune_cols, variable.name = "gene", value.name = "expression")
immune_long <- merge(immune_long, manifest[, .(sample = gsm, group3)], by = "sample")
immune_long[, `:=`(contrast, factor(ifelse(group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD")))]
immune_long[, `:=`(gene, factor(gene, levels = immune_cols))]
immune_stats <- immune_long[, {
    pvalue <- if (uniqueN(expression) < 2) 
        NA_real_
    else wilcox.test(expression ~ contrast, exact = FALSE)$p.value
    .(HC_n = sum(contrast == "HC"), IBD_n = sum(contrast == "IBD"), HC_median = median(expression[contrast == "HC"]), IBD_median = median(expression[contrast == 
        "IBD"]), pvalue = pvalue, y = max(expression) + 0.045)
}, by = gene]
immune_stats[, `:=`(FDR, p.adjust(pvalue, "BH"))]
immune_stats[, `:=`(pstar, fifelse(!is.na(pvalue) & pvalue < 0.001, "***", fifelse(!is.na(pvalue) & pvalue < 0.01, "**", 
    fifelse(!is.na(pvalue) & pvalue < 0.05, "*", ""))))]
immune_stats_export <- copy(immune_stats)
setnames(immune_stats_export, "gene", "cell")
fwrite(immune_stats_export, file.path(immune_work, "immune_IBD_vs_HC_wilcoxon.tsv"), sep = "\t")
p <- ggplot(immune_long, aes(x = gene, y = expression, fill = contrast, color = contrast)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + 
    scale_fill_manual(name = "Group", values = c(HC = "deepskyblue", IBD = "hotpink")) + scale_color_manual(name = "Group", 
    values = c(HC = "dodgerblue", IBD = "plum3")) + geom_text(data = immune_stats[nzchar(pstar)], aes(x = gene, y = y, label = pstar), 
    inherit.aes = FALSE, size = 5) + theme_bw() + labs(x = "", y = "Normalized ssGSEA score") + theme(axis.text.x = element_text(vjust = 1, 
    size = 9, hjust = 1, colour = "black"), legend.position = "top") + rotate_x_text(45) + coord_cartesian(ylim = c(0, 1.08))
ggsave(pdf_path("12.Immune microenvironment", "3.Immune comparison.pdf"), width = 16, height = 7, plot = p)
scrna_work <- file.path(root, "work/15_key_gene_scrna_expression")
scrna_plot <- fread(file.path(scrna_work, "key_cell_sample_average_expression.tsv"))
scrna_plot[, `:=`(contrast, factor(ifelse(group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD")))]
scrna_stats <- scrna_plot[, .(HC_n = sum(contrast == "HC"), IBD_n = sum(contrast == "IBD"), HC_median = median(expression[contrast == 
    "HC"]), IBD_median = median(expression[contrast == "IBD"]), pvalue = wilcox.test(expression ~ contrast, exact = TRUE)$p.value), 
    by = gene]
scrna_stats[, `:=`(FDR, p.adjust(pvalue, "BH"))]
fwrite(scrna_stats, file.path(scrna_work, "key_cell_sample_IBD_vs_HC_tests.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "logs/23_final_figures_sessionInfo.txt"))
