source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/12_immune_regulators")
result_dir <- file.path(root, "results/12_Immune_Regulators")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
mat <- readRDS(file.path(root, "work/10_key_gene_enrichment/bulk_log2_FPKM_matrix.rds"))
manifest <- fread(file.path(root, "work/10_key_gene_enrichment/bulk_sample_manifest.tsv"))
key_genes <- intersect(fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene, rownames(mat))
reg <- fread(file.path(client_root, "data", "reference", "Immunomodulator_and_chemokines.txt"))
reg <- reg[Id %in% rownames(mat)]
groups <- list(All = manifest$gsm, HC = manifest[group3 == "HC", gsm], UC = manifest[group3 == "UC", gsm], CD = manifest[group3 == 
    "CD", gsm], IBD = manifest[group3 != "HC", gsm])
res <- rbindlist(lapply(names(groups), function(gr) rbindlist(lapply(key_genes, function(g) rbindlist(lapply(reg$Id, function(r) {
    s <- groups[[gr]]
    ct <- cor.test(as.numeric(mat[g, s]), as.numeric(mat[r, s]), method = "spearman", exact = FALSE)
    data.table(group = gr, gene = g, immuneGene = r, cor = unname(ct$estimate), pvalue = ct$p.value)
}))))))
res[, `:=`(FDR, p.adjust(pvalue, "BH")), by = group]
write.table(res, file.path(work_dir, "key_gene_immune_regulator_spearman.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capitalize_first <- function(string) {
    words <- strsplit(string, " ")[[1]]
    paste0(toupper(substring(words, 1, 1)), substring(words, 2), collapse = " ")
}
for (i in seq_along(unique(reg$type))) {
    type_i <- unique(reg$type)[i]
    ids <- reg[type == type_i, Id]
    data2 <- as.data.frame(res[group == "All" & immuneGene %in% ids])
    data2$pstar <- ifelse(data2$FDR < 0.001, "***", ifelse(data2$FDR < 0.01, "**", ifelse(data2$FDR < 0.05, "*", "")))
    data2 <- na.omit(data2)
    n_regulators <- length(unique(data2$immuneGene))
    dotheight <- if (n_regulators > 40) 
        n_regulators * 0.18
    else if (n_regulators < 20) 
        10
    else 15
    p <- ggplot(data2, aes(immuneGene, gene)) + geom_tile(aes(fill = cor)) + geom_text(aes(label = pstar), color = "black", 
        size = 4) + scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c", limit = c(-1, 1), name = paste0("*    FDR < 0.05", 
        "\n\n", "**  FDR < 0.01", "\n\n", "*** FDR < 0.001", "\n\n", "Spearman")) + labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 12, 
        angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.ticks.y = element_blank(), 
        panel.background = element_blank()) + coord_fixed()
    ggsave(file.path(result_dir, paste0(i, ".", capitalize_first(type_i), ".pdf")), width = dotheight, height = max(4, length(key_genes) * 
        0.5), plot = p)
}
writeLines(capture.output(sessionInfo()), file.path(root, "logs/12_immune_regulators_sessionInfo.txt"))
