source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(patchwork)
    library(ggplot2)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/08_cellchat")
result_dir <- file.path(root, "results/08_CellChat")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
key_cell <- readLines(file.path(root, "work/05_key_cell_differential_gsea/key_cell_type.txt"), n = 1)
scRNA$group4 <- ifelse(scRNA$group3 == "HC", "HC", "IBD")
for (z in c("HC", "UC", "CD")) {
    out <- file.path(work_dir, "CellChat", z)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    setwd(out)
    cellchat_func(scRNA[, scRNA$group3 == z], celltype = "cellType_1", ShowCell = key_cell, OnlyPlot = FALSE, AllPlot = FALSE)
}
out <- file.path(work_dir, "CellChat", "IBD")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
setwd(out)
cellchat_func(scRNA[, scRNA$group4 == "IBD"], celltype = "cellType_1", ShowCell = key_cell, OnlyPlot = FALSE, AllPlot = FALSE)
for (z in c("HC", "UC", "CD", "IBD")) {
    for (i in c("1.Net_number_strength.pdf", "2.Interaction Count.pdf", "3.Bubble.pdf")) file.copy(file.path(work_dir, "CellChat", 
        z, i), file.path(result_dir, paste0(z, "_", i)), overwrite = TRUE)
}
pairs <- list(c("HC", "UC"), c("HC", "CD"), c("UC", "CD"), c("HC", "IBD"))
for (v in pairs) {
    group_by <- if ("IBD" %in% v) 
        "group4"
    else "group3"
    cellchat_module(object = scRNA, OutPath = work_dir, group_by = group_by, celltype = "cellType_1", idents = v, ShowCell = key_cell, 
        OnlyPlot = TRUE, AllPlot = FALSE)
    prefix <- paste(v, collapse = "_vs_")
    for (i in c("1.CompareInteractions.pdf", "2.NetVisual.pdf", "3.NetVisualHeatmap.pdf", "4.SignalingRole.pdf", "5.Bubble.pdf")) file.copy(file.path(work_dir, 
        "CellChat", i), file.path(result_dir, paste0(prefix, "_", i)), overwrite = TRUE)
}
writeLines(capture.output(sessionInfo()), file.path(root, "logs/08_cellchat_sessionInfo.txt"))
