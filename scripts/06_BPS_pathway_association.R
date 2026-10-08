source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(GSEABase)
    library(GSVA)
    library(msigdbr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/06_bps_pathway_association")
result_dir <- file.path(root, "results/06_BPS_Pathway_Association")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
key_cell <- readLines(file.path(root, "work/05_key_cell_differential_gsea/key_cell_type.txt"), n = 1)
scRNA <- subset(scRNA, subset = cellType_1 == key_cell)
cell_bps <- read.table(file.path(root, "work/04_scBPS_key_cell/scBPS/output/norm_score.tsv"), header = TRUE, row.names = 1, 
    sep = "\t", check.names = FALSE)
cell_anno <- read.table(file.path(root, "work/04_scBPS_key_cell/scBPS/output/cell.annotation.txt"), header = TRUE, sep = "\t")
key_auc <- read.table(file.path(root, "work/04_scBPS_key_cell/scBPS/key_cell_BPS_AUC.txt"), header = TRUE, check.names = FALSE)
key_auc <- key_auc[key_auc$class != "Nonsig", , drop = FALSE]
key_auc <- key_auc[order(key_auc$FDR, -key_auc$AUC), , drop = FALSE]
stopifnot(nrow(key_auc) > 0)
key_trait <- as.character(key_auc$trait[1])
cell_id <- intersect(colnames(scRNA), rownames(cell_bps))
stopifnot(length(cell_id) == ncol(scRNA))
cell_bps <- cell_bps[cell_id, , drop = FALSE]
scRNA <- scRNA[, cell_id]
kegg_res <- read.table(file.path(root, "work/05_key_cell_differential_gsea/IBD_vs_HC_KEGG_GSEA.tsv"), header = TRUE, sep = "\t", 
    quote = "", fill = TRUE, check.names = FALSE)
hm_res <- read.table(file.path(root, "work/05_key_cell_differential_gsea/IBD_vs_HC_HALLMARK_GSEA.tsv"), header = TRUE, sep = "\t", 
    quote = "", fill = TRUE, check.names = FALSE)
kegg_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_sets <- split(kegg_df$gene_symbol, kegg_df$gs_name)
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
normalize_kegg <- function(x) paste0("KEGG_", toupper(gsub("[^A-Za-z0-9]+", "_", x)))
kegg_res$set <- normalize_kegg(kegg_res$Description)
hm_res$set <- hm_res$ID
cand <- rbind(data.frame(set = kegg_res$set, p.adjust = kegg_res$p.adjust, source = "KEGG"), data.frame(set = hm_res$set, 
    p.adjust = hm_res$p.adjust, source = "HALLMARK"))
cand <- cand[is.finite(cand$p.adjust) & cand$set %in% c(names(kegg_sets), names(hallmark_sets)), , drop = FALSE]
cand <- cand[order(cand$p.adjust), , drop = FALSE]
target_pathways <- unique(cand$set)[seq_len(min(4, length(unique(cand$set))))]
stopifnot(length(target_pathways) == 4)
geneset_list <- c(kegg_sets[target_pathways[target_pathways %in% names(kegg_sets)]], hallmark_sets[target_pathways[target_pathways %in% 
    names(hallmark_sets)]])
geneset_list <- geneset_list[target_pathways]
exp <- as.matrix(GetAssayData(scRNA, assay = "RNA", slot = "data"))
gsva_matrix <- gsva(expr = exp, gset.idx.list = geneset_list, kcdf = "Gaussian", method = "ssgsea", abs.ranking = TRUE)
cell_ssgsea <- as.data.frame(t(gsva_matrix))
write.table(cell_ssgsea, file.path(work_dir, "all_ssgsea_results.tsv"), sep = "\t", quote = FALSE)
merged_df <- data.frame(BPS = cell_bps[[key_trait]], cell_ssgsea, check.names = FALSE)
merged_df$sample <- scRNA$orig.ident
merged_df$bin <- ntile(merged_df$BPS, 5)
merged_df$bin_label <- c(2, 4, 6, 8, 10)[merged_df$bin]
bin_pathway_mean <- merged_df %>% group_by(bin_label) %>% summarise(across(all_of(target_pathways), mean), .groups = "drop") %>% 
    pivot_longer(cols = -bin_label, names_to = "Pathway", values_to = "Mean_Score")
sample_summary <- merged_df %>% group_by(sample) %>% summarise(BPS = mean(BPS), across(all_of(target_pathways), mean), .groups = "drop")
cor_tab <- do.call(rbind, lapply(target_pathways, function(z) {
    t <- cor.test(sample_summary$BPS, sample_summary[[z]], method = "spearman", exact = FALSE)
    data.frame(trait = key_trait, pathway = z, rho = unname(t$estimate), pvalue = t$p.value)
}))
cor_tab$FDR <- p.adjust(cor_tab$pvalue, "BH")
write.table(bin_pathway_mean, file.path(work_dir, "cell_quintile_plotdata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_summary, file.path(work_dir, "sample_level_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cor_tab, file.path(work_dir, "sample_level_spearman.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(key_cell, key_trait, target_pathways), file.path(work_dir, "selected_inputs.txt"))
color_map <- setNames(c("#6A51A3", "#D84A92", "#F0AD4E", "#1F77B4"), target_pathways)
p <- ggplot(bin_pathway_mean, aes(x = bin_label, y = Mean_Score, color = Pathway, group = Pathway)) + geom_line(linewidth = 1.2) + 
    geom_point(size = 4) + scale_color_manual(values = color_map) + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(1.5, 
    10.5), name = paste0(key_trait, " BPS in ", key_cell, " (quintile bins)")) + scale_y_continuous(name = "Mean ssGSEA pathway score") + 
    theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), 
    axis.text = element_text(size = 12), legend.title = element_blank(), legend.text = element_text(size = 11), legend.position = c(0.58, 
        0.22), panel.border = element_rect(linewidth = 1.1))
ggsave(file.path(result_dir, "1.BPS_pathway_trend_original.pdf"), p, width = 12, height = 8)
writeLines(capture.output(sessionInfo()), file.path(root, "logs/06_bps_pathway_association_sessionInfo.txt"))
