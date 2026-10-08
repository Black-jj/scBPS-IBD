source(file.path("R", "utils.R"))
project_dir <- client_root
work_dir <- file.path(project_dir, "work", "02_scRNA_annotation_abundance")
result_dir <- file.path(project_dir, "results", "02_Cell_Annotation_Abundance")
log_dir <- file.path(project_dir, "logs")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
scRNA <- readRDS(file.path(project_dir, "work", "01_scRNA_qc_integration", "2.Cluster", "data.merge.harmony.2000.rds"))
cluster_col <- "RNA_snn_res.0.1"
cell_map <- c(`0` = "T_cells", `1` = "Plasma_cells", `2` = "Enterocytes", `3` = "Fibroblasts", `4` = "Myeloid", `5` = "B_cells", 
    `6` = "Goblet_cells", `7` = "Mast_cells", `8` = "Endothelial", `9` = "Proliferating_Tcells", `10` = "Smooth_muscle_cells", 
    `11` = "Tuft_cells", `12` = "Glial_cells")
stopifnot(setequal(unique(as.character(scRNA@meta.data[[cluster_col]])), names(cell_map)))
scRNA$cellType_1 <- factor(unname(cell_map[as.character(scRNA@meta.data[[cluster_col]])]), levels = unname(cell_map))
write.table(data.frame(cluster = names(cell_map), cellType_1 = unname(cell_map)), file.path(work_dir, "cluster_annotation.tsv"), 
    sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(scRNA, file.path(work_dir, "scRNA_annotated.rds"), compress = FALSE)
prop <- as.data.frame(table(sample = scRNA$orig.ident, group3 = scRNA$group3, cellType_1 = scRNA$cellType_1), stringsAsFactors = FALSE)
prop <- prop[prop$Freq > 0, ]
prop$Percent <- ave(prop$Freq, prop$sample, FUN = function(z) 100 * z/sum(z))
write.table(prop, file.path(work_dir, "sample_celltype_proportions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stats <- do.call(rbind, lapply(split(prop, prop$cellType_1), function(d) {
    z <- tryCatch(kruskal.test(Percent ~ group3, data = d), error = function(e) NULL)
    data.frame(cellType_1 = as.character(d$cellType_1[1]), p_value = if (is.null(z)) 
        NA_real_
    else z$p.value)
}))
stats$FDR <- p.adjust(stats$p_value, "BH")
write.table(stats, file.path(work_dir, "sample_level_abundance_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
features <- c("MS4A1", "CD79A", "CD79B", "VWF", "PECAM1", "KDR", "KRT20", "PIGR", "FABP1", "COL1A1", "DCN", "LUM", "SOX10", 
    "S100B", "PLP1", "MUC2", "TFF3", "SPINK4", "CPA3", "KIT", "TPSAB1", "LYZ", "TYMP", "CD14", "MZB1", "JCHAIN", "IGKC", 
    "MKI67", "TOP2A", "STMN1", "IL7R", "CD3D", "CD3E", "POU2F3", "TRPM5", "IL17RB", "RGS5", "ACTA2", "MYH11", "MPZ")
features <- intersect(features, rownames(scRNA))
p <- DimPlot(scRNA, group.by = cluster_col, label = TRUE, label.size = 3, reduction = "umap", raster = FALSE) + NoLegend()
ggsave_fun(filename = file.path(result_dir, "1.RNA_snn_res.0.1_umap"), plot = p, width = 5, height = 5)
p <- DimPlot(scRNA, group.by = "cellType_1", label = TRUE, label.size = 3, reduction = "umap", raster = FALSE) + NoLegend()
ggsave_fun(filename = file.path(result_dir, "2.Anno_umap"), plot = p, width = 5, height = 5)
p <- DotPlot(scRNA, features = features, cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 6)) + 
    theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "top", 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, "cm"))
ggsave_fun(filename = file.path(result_dir, "3.Anno_Dotplot"), plot = p, width = 13, height = 7)
my_pal2 <- c("#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda", 
    "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#00B3F1", "#00538A")
p <- cellRatioPlot(object = scRNA, sample.name = "group3", celltype.name = "cellType_1", flow.curve = 0.5, fill.col = my_pal2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave_fun(filename = file.path(result_dir, "4.Anno_group_ratio"), plot = p, width = 7, height = 7)
p <- cellRatioPlot(object = scRNA, sample.name = "orig.ident", celltype.name = "cellType_1", flow.curve = 0.5, fill.col = my_pal2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave_fun(filename = file.path(result_dir, "5.Anno_sample_ratio"), plot = p, width = 12, height = 7)
stopifnot(all(tolower(tools::file_ext(list.files(result_dir))) == "pdf"))
capture.output(sessionInfo(), file = file.path(log_dir, "02_scRNA_annotation_abundance_sessionInfo.txt"))
