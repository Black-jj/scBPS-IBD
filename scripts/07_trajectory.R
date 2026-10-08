source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(CytoTRACE)
    library(monocle)
    library(dplyr)
    library(ggplot2)
    library(ggplotify)
    library(RColorBrewer)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/07_key_cell_trajectory")
result_dir <- file.path(root, "results/07_Key_Cell_Trajectory")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
key_cell <- readLines(file.path(root, "work/05_key_cell_differential_gsea/key_cell_type.txt"), n = 1)
subRNA <- subset(scRNA, subset = cellType_1 == key_cell)
mat_full <- as.matrix(GetAssayData(subRNA, assay = "RNA", slot = "counts"))
results <- CytoTRACE(mat = mat_full, enableFast = FALSE, ncores = 64)
save(results, file = file.path(work_dir, "cytotrace_result_full.Robj"))
phe <- as.character(subRNA$group3)
names(phe) <- colnames(subRNA)
emb <- Embeddings(subRNA, "umap")
cyto_dir <- file.path(work_dir, "CytoTRACE")
dir.create(cyto_dir, showWarnings = FALSE)
plotCytoTRACE(results, phenotype = phe, emb = emb, outputDir = cyto_dir)
cds <- monocle2_module(object = subRNA, OutPath = work_dir, group_by = "group3", track_gene = c("VariableFeatures", 2000), 
    ShowGene = NA, OnlyPlot = FALSE)
ct <- results$CytoTRACE
if (is.null(names(ct))) names(ct) <- colnames(subRNA)
pData(cds)$CytoTRACE <- ct[colnames(cds)]
state_score <- aggregate(pData(cds)$CytoTRACE, list(State = pData(cds)$State), median, na.rm = TRUE)
root_state <- state_score$State[which.max(state_score$x)]
cds <- orderCells(cds, root_state = root_state)
saveRDS(cds, file.path(work_dir, "Monocle2/cds_rerooted.rds"))
write.table(state_score, file.path(work_dir, "root_state_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.character(root_state), file.path(work_dir, "selected_root_state.txt"))
p <- plot_cell_trajectory(cds, color_by = "CytoTRACE", cell_size = 0.15) + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, 
    "Spectral"))) + theme(legend.position = "right", text = element_text(face = "bold"))
ggsave(file.path(result_dir, "1.CytoTRACE_on_trajectory.pdf"), plot = p, height = 3.5, width = 4.8)
p <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.1) + ggsci::scale_color_gsea() + theme(legend.position = "right", 
    text = element_text(face = "bold"))
ggsave(file.path(result_dir, "2.Trajectory.Pseudotime.pdf"), plot = p, height = 4, width = 4.5)
p <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.1) + guides(color = guide_legend(override.aes = list(alpha = 1, 
    size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right", text = element_text(face = "bold"))
ggsave(file.path(result_dir, "3.Trajectory.State.pdf"), plot = p, height = 4, width = 4.5)
p <- plot_cell_trajectory(cds, color_by = "Cluster", cell_size = 0.1) + guides(color = guide_legend(override.aes = list(alpha = 1, 
    size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right", text = element_text(face = "bold"))
ggsave(file.path(result_dir, "4.Trajectory.Group.pdf"), plot = p, height = 4, width = 4.5)
pseudo_time_diff <- differentialGeneTest(cds, cores = 64, fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudo_time_diff <- arrange(pseudo_time_diff, qval)
write.table(cbind(gene = rownames(pseudo_time_diff), pseudo_time_diff), file.path(work_dir, "pseudotime_genes.tsv"), sep = "\t", 
    quote = FALSE, row.names = FALSE)
plot_heatmap_gene <- rownames(pseudo_time_diff)[seq_len(min(50, nrow(pseudo_time_diff)))]
p <- monocle::plot_pseudotime_heatmap(cds[plot_heatmap_gene, ], num_clusters = 3, show_rownames = TRUE, return_heatmap = TRUE, 
    hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
ggsave(file.path(result_dir, "5.Pseudotime_Heatmap.pdf"), plot = ggplotify::as.ggplot(p), width = 5, height = length(plot_heatmap_gene) * 
    0.12)
writeLines(capture.output(sessionInfo()), file.path(root, "logs/07_key_cell_trajectory_sessionInfo.txt"))
