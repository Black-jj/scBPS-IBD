source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(GSVA)
    library(GSEABase)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(corrplot)
    library(aplot)
    library(IOBR)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/11_immune_microenvironment")
result_dir <- file.path(root, "results/12.Immune microenvironment")
gmt_file <- file.path(root, "data/reference/ssGSEA/immune.gmt")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
stopifnot(file.exists(gmt_file))
mat <- readRDS(file.path(root, "work/10_key_gene_enrichment/bulk_log2_FPKM_matrix.rds"))
manifest <- fread(file.path(root, "work/10_key_gene_enrichment/bulk_sample_manifest.tsv"))
key_genes <- fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene
mat <- avereps(as.matrix(mat))
mat <- mat[rowMeans(mat) > 0, , drop = FALSE]
geneSet <- getGmt(gmt_file, geneIdType = SymbolIdentifier())
gene_set_list <- lapply(geneSet, geneIds)
names(gene_set_list) <- vapply(geneSet, setName, character(1))
gene_set_coverage <- rbindlist(lapply(names(gene_set_list), function(signature) {
    supplied <- unique(gene_set_list[[signature]])
    data.table(signature = signature, supplied_genes = length(supplied), matched_genes = sum(supplied %in% rownames(mat)))
}))
fwrite(gene_set_coverage, file.path(work_dir, "ssGSEA_gene_set_coverage.tsv"), sep = "\t")
ssgseaScore <- gsva(mat, geneSet, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
normalize <- function(x) {
    return((x - min(x))/(max(x) - min(x)))
}
ssgseaOut <- normalize(ssgseaScore)
stopifnot(all(is.finite(ssgseaOut)), ncol(ssgseaOut) == nrow(manifest))
ssgsea_original_format <- rbind(id = colnames(ssgseaOut), ssgseaOut)
write.table(ssgsea_original_format, file = file.path(work_dir, "ssgseaOut.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
immune <- t(ssgseaOut)
immune <- immune[manifest$gsm, , drop = FALSE]
stopifnot(identical(rownames(immune), manifest$gsm))
fwrite(data.table(sample = rownames(immune), as.data.frame(immune, check.names = FALSE)), file.path(work_dir, "ssGSEA_scores_by_sample.tsv"), 
    sep = "\t")
datGroup <- manifest[, .(Acc = gsm, Tissue = group3)]
datGroup$Tissue <- factor(datGroup$Tissue, levels = c("HC", "UC", "CD"))
datGroup <- datGroup[order(datGroup$Tissue), ]
datGroup$Acc <- factor(datGroup$Acc, levels = datGroup$Acc)
data <- data.frame(immune, Acc = rownames(immune), check.names = FALSE)
data <- inner_join(datGroup, data, by = "Acc")
data_p <- reshape2::melt(data, id.vars = c("Acc", "Tissue"))
data_p$Acc <- factor(data_p$Acc, levels = datGroup$Acc)
data_p$variable <- factor(data_p$variable, levels = rev(colnames(immune)))
datGroup$p <- "Group"
p2 <- ggplot(datGroup, aes(Acc, p, fill = Tissue)) + geom_tile() + scale_fill_manual(values = c(HC = "#4682B4", UC = "#E69F00", 
    CD = "#CD2626")) + scale_y_discrete(position = "right") + theme_minimal() + xlab(NULL) + ylab(NULL) + theme(text = element_text(size = 15), 
    axis.text.x = element_blank()) + labs(fill = "Group")
p1 <- ggplot(data_p, aes(Acc, variable, fill = value)) + geom_tile() + labs(x = NULL, y = NULL, fill = "Normalized\nssGSEA score") + 
    scale_fill_gradientn(colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"), limits = c(0, 1)) + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 9, color = "black"), 
        panel.grid = element_blank())
p <- p1 %>% aplot::insert_top(p2, height = 0.05)
ggsave(file.path(result_dir, "1.Immune infiltration.pdf"), plot = p, width = 15, height = 10, limitsize = FALSE)
cor_data <- immune[, apply(immune, 2, sd) > 0, drop = FALSE]
pdf(file.path(result_dir, "2.Immune correlation.pdf"), height = 15, width = 15)
corrplot(corr = cor(cor_data), method = "color", order = "hclust", tl.col = "black", addCoef.col = "black", number.cex = 0.55, 
    tl.cex = 0.75, col = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
outTab <- rbindlist(lapply(colnames(immune), function(i) data.table(expression = immune[, i], group3 = manifest$group3, gene = i)))
outTab[, `:=`(contrast, factor(ifelse(group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD")))]
outTab[, `:=`(gene, factor(gene, levels = colnames(immune)))]
immune_stats <- outTab[, {
    pvalue <- if (uniqueN(expression) < 2) 
        NA_real_
    else wilcox.test(expression ~ contrast, exact = FALSE)$p.value
    .(HC_n = sum(contrast == "HC"), IBD_n = sum(contrast == "IBD"), HC_median = median(expression[contrast == "HC"]), IBD_median = median(expression[contrast == 
        "IBD"]), pvalue = pvalue, y = max(expression) + 0.045)
}, by = gene]
immune_stats[, `:=`(FDR, p.adjust(pvalue, "BH"))]
immune_stats[, `:=`(pstar, fifelse(!is.na(pvalue) & pvalue < 0.001, "***", fifelse(!is.na(pvalue) & pvalue < 0.01, "**", 
    fifelse(!is.na(pvalue) & pvalue < 0.05, "*", ""))))]
fwrite(immune_stats, file.path(work_dir, "immune_IBD_vs_HC_wilcoxon.tsv"), sep = "\t")
p <- ggplot(outTab, aes(x = gene, y = expression, fill = contrast, color = contrast)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + 
    scale_fill_manual(name = "Group", values = c(HC = "deepskyblue", IBD = "hotpink")) + scale_color_manual(name = "Group", 
    values = c(HC = "dodgerblue", IBD = "plum3")) + geom_text(data = immune_stats[nzchar(pstar)], aes(x = gene, y = y, label = pstar), 
    inherit.aes = FALSE, size = 5) + theme_bw() + labs(x = "", y = "Normalized ssGSEA score") + theme(axis.text.x = element_text(vjust = 1, 
    size = 9, hjust = 1, colour = "black"), legend.position = "top") + rotate_x_text(45) + coord_cartesian(ylim = c(0, 1.08))
ggsave(file.path(result_dir, "3.Immune comparison.pdf"), width = 16, height = 7, plot = p)
cor_tab <- rbindlist(lapply(intersect(key_genes, rownames(mat)), function(gene) rbindlist(lapply(colnames(immune), function(signature) {
    dd <- cor.test(as.numeric(immune[, signature]), as.numeric(mat[gene, rownames(immune)]), method = "spearman", exact = FALSE)
    data.table(gene = gene, immune_signature = signature, cor = unname(dd$estimate), p.value = dd$p.value)
}))))
cor_tab[, `:=`(FDR, p.adjust(p.value, "BH"))]
cor_tab[, `:=`(pstar, fifelse(FDR < 0.001, "***", fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", ""))))]
fwrite(cor_tab, file.path(work_dir, "key_gene_immune_spearman.tsv"), sep = "\t")
p <- ggplot(cor_tab, aes(immune_signature, gene)) + geom_tile(aes(fill = cor)) + geom_text(aes(label = pstar), color = "black", 
    size = 4) + scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c", limit = c(-1, 1), name = paste0("*    FDR < 0.05", 
    "\n\n", "**  FDR < 0.01", "\n\n", "*** FDR < 0.001", "\n\n", "Spearman")) + labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 8, 
    angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(size = 8, color = "black"), axis.ticks.y = element_blank(), 
    panel.background = element_blank()) + coord_fixed()
ggsave(file.path(result_dir, "4.Key-gene correlation.pdf"), plot = p, width = 11, height = max(4, length(unique(cor_tab$gene)) * 
    0.5))
writeLines(capture.output(sessionInfo()), file.path(root, "logs/11_immune_microenvironment_sessionInfo.txt"))
