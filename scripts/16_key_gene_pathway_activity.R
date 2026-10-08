source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(AUCell)
    library(GSVA)
    library(GSEABase)
    library(limma)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/16_key_gene_pathway_activity")
result_dir <- file.path(root, "results/16_Key_Gene_Pathway_Activity")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
KeyGene <- fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
key_cell <- readLines(file.path(root, "work/05_key_cell_differential_gsea/key_cell_type.txt"), n = 1)
subRNA <- subset(scRNA, subset = cellType_1 == key_cell)
dir.create(file.path(work_dir, "data"), showWarnings = FALSE)
file.symlink(file.path(client_root, "data", "reference", "hallmark.gs.RData"), file.path(work_dir, "data/hallmark.gs.rdata"))
file.symlink(file.path(client_root, "data", "reference", "GenesetInfo.txt"), file.path(work_dir, "data/GenesetInfo.txt"))
setwd(work_dir)
hallmark_module(outpath = work_dir, SeuratObject = subRNA, KeyGene = KeyGene)
file.copy(file.path(work_dir, "Immunometabolic pathways/1.Hallmark.pdf"), file.path(result_dir, "1.Key_Cell_AUCell_Hallmark.pdf"), 
    overwrite = TRUE)
mat <- readRDS(file.path(root, "work/10_key_gene_enrichment/bulk_log2_FPKM_matrix.rds"))
manifest <- fread(file.path(root, "work/10_key_gene_enrichment/bulk_sample_manifest.tsv"))
load(file.path(client_root, "data", "reference", "hallmark.gs.RData"))
immune_gmt <- getGmt(file.path(client_root, "data", "reference", "immune.gmt"), geneIdType = SymbolIdentifier())
immune_gs <- lapply(immune_gmt, geneIds)
score_list <- list(Hallmark = gsva(as.matrix(mat), gs), Immune = gsva(as.matrix(mat), immune_gs))
for (nm in names(score_list)) {
    score <- score_list[[nm]]
    long <- reshape2::melt(score, varnames = c("Pathway", "sample"), value.name = "score")
    long <- left_join(long, manifest[, .(sample = gsm, group3)], by = "sample")
    tests <- long %>% group_by(Pathway) %>% summarise(pvalue = kruskal.test(score ~ group3)$p.value, .groups = "drop")
    tests$FDR <- p.adjust(tests$pvalue, "BH")
    write.table(long, file.path(work_dir, paste0(nm, "_bulk_GSVA_scores.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(tests, file.path(work_dir, paste0(nm, "_bulk_group_tests.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    long$group3 <- factor(long$group3, levels = c("HC", "UC", "CD"))
    long$Pathway <- as.factor(long$Pathway)
    p <- ggplot(long, aes(x = Pathway, y = score, fill = group3, color = group3)) + geom_boxplot(alpha = 0.3) + scale_fill_manual(name = "Group", 
        values = c(HC = "deepskyblue", UC = "#E69F00", CD = "hotpink")) + scale_color_manual(name = "Group", values = c(HC = "dodgerblue", 
        UC = "#CC7900", CD = "plum3")) + theme_bw() + labs(x = "", y = paste0(nm, " GSVA score")) + theme(axis.text.x = element_text(vjust = 1, 
        size = 8, hjust = 1, colour = "black"), legend.position = "top") + rotate_x_text(45) + stat_compare_means(aes(group = group3), 
        symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), label = "p.signif", 
        method = "kruskal.test", show.legend = FALSE)
    ggsave(file.path(result_dir, paste0(ifelse(nm == "Hallmark", 2, 3), ".Bulk_", nm, "_GSVA_Group.pdf")), width = max(15, 
        length(unique(long$Pathway)) * 0.28), height = 7, plot = p, limitsize = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(root, "logs/16_key_gene_pathway_activity_sessionInfo.txt"))
