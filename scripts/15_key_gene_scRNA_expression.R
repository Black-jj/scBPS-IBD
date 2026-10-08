source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(ggpubr)
    library(data.table)
    library(dplyr)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/15_key_gene_scrna_expression")
result_dir <- file.path(root, "results/15_Key_Gene_scRNA_Expression")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
KeyGene <- intersect(fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene, rownames(scRNA))
stopifnot(length(KeyGene) > 0)
GeneExpre(outpath = work_dir, SeuratObject = scRNA, KeyGene = KeyGene)
shared_dir <- file.path(work_dir, "KeyGene expression abundance")
file.copy(file.path(shared_dir, "1.Keygene.FeaturePlot.pdf"), file.path(result_dir, "1.Keygene.FeaturePlot.pdf"), overwrite = TRUE)
file.copy(file.path(shared_dir, "2.Keygene.DotPlot.pdf"), file.path(result_dir, "2.Keygene.DotPlot.pdf"), overwrite = TRUE)
Idents(scRNA) <- "group3"
p <- VlnPlot(scRNA, features = KeyGene, ncol = 2, raster = FALSE)
ggsave(file.path(result_dir, "3.Keygene.Group_Vlnplot.pdf"), plot = p, width = 10, height = max(8, ceiling(length(KeyGene)/2) * 
    4))
key_cell <- readLines(file.path(root, "work/05_key_cell_differential_gsea/key_cell_type.txt"), n = 1)
subRNA <- subset(scRNA, subset = cellType_1 == key_cell)
avg <- AverageExpression(subRNA, assays = "RNA", features = KeyGene, group.by = "orig.ident", slot = "data", verbose = FALSE)$RNA
sample_group <- unique(data.frame(sample = subRNA$orig.ident, group3 = subRNA$group3))
plotdata <- reshape2::melt(as.matrix(avg), varnames = c("gene", "sample"), value.name = "expression")
plotdata <- left_join(plotdata, sample_group, by = "sample")
write.table(plotdata, file.path(work_dir, "key_cell_sample_average_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
plotdata$contrast <- factor(ifelse(plotdata$group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD"))
sample_stats <- as.data.table(plotdata)[, .(HC_n = sum(contrast == "HC"), IBD_n = sum(contrast == "IBD"), HC_median = median(expression[contrast == 
    "HC"]), IBD_median = median(expression[contrast == "IBD"]), pvalue = wilcox.test(expression ~ contrast, exact = TRUE)$p.value), 
    by = gene]
sample_stats[, `:=`(FDR, p.adjust(pvalue, "BH"))]
write.table(sample_stats, file.path(work_dir, "key_cell_sample_IBD_vs_HC_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
display_genes <- sample_stats[FDR < 0.05, gene]
stopifnot(length(display_genes) > 0)
plotdata <- plotdata[plotdata$gene %in% display_genes, , drop = FALSE]
p <- ggplot(plotdata, aes(x = contrast, y = expression, fill = contrast, color = contrast)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + 
    geom_jitter(width = 0.08, size = 2) + facet_wrap(~gene, scales = "free_y", nrow = 1) + scale_fill_manual(values = c(HC = "deepskyblue", 
    IBD = "hotpink")) + scale_color_manual(values = c(HC = "dodgerblue", IBD = "plum3")) + stat_compare_means(method = "wilcox.test", 
    label = "p.format") + theme_bw() + labs(x = "", y = paste0("Sample-average expression in ", key_cell)) + theme(legend.position = "none")
ggsave(file.path(result_dir, "4.Keygene_Sample_Level_Boxplot.pdf"), width = max(7, length(display_genes) * 3.3), height = 5, 
    plot = p)
writeLines(capture.output(sessionInfo()), file.path(root, "logs/15_key_gene_scrna_expression_sessionInfo.txt"))
