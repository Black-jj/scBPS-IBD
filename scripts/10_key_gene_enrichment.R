source(file.path("R", "utils.R"))
stopifnot(all(tolower(tools::file_ext(list.files(result_dir))) == "pdf"))
options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout")))
suppressPackageStartupMessages({
    library(openxlsx)
    library(limma)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(GSEABase)
    library(GSVA)
    library(ggplot2)
    library(stringr)
    library(data.table)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/10_key_gene_enrichment")
result_dir <- file.path(root, "results/10_Key_Gene_Enrichment")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
key_genes <- fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene
bulk_df <- read.xlsx(file.path(root, "data/GSE235236/GSE235236_RAW/GSE235236_FPKM_expression_matrix.xlsx"), sheet = "FPKM_matrix")
symbol <- sub("^.*_", "", bulk_df$gene_id)
mat <- avereps(as.matrix(bulk_df[, -1]), ID = symbol)
manifest <- fread(file.path(root, "work/00_preflight/GSE235236_sample_manifest.tsv"))
mat <- log2(mat[, manifest$expression_column, drop = FALSE] + 1)
colnames(mat) <- manifest$gsm
saveRDS(mat, file.path(work_dir, "bulk_log2_FPKM_matrix.rds"), compress = FALSE)
write.table(manifest, file.path(work_dir, "bulk_sample_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ids <- bitr(intersect(key_genes, rownames(mat)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
bulk_key_genes <- intersect(c("CCL20", "CLDN4", "CXCL8"), rownames(mat))
kegg_snapshot <- file.path(client_root, "data", "reference", "kegg.all.entrez.hsa.rds")
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
key_gene_gsea_work <- file.path(work_dir, "key_gene_KEGG_GSEA")
dir.create(key_gene_gsea_work, recursive = TRUE, showWarnings = FALSE)
gsea_display_used <- character()
gsea_display_selected <- list()
GSEA_module(mat = mat, outpath = key_gene_gsea_work, species = "human", KeyGene = bulk_key_genes, isKEGG = TRUE)
selected_display <- rbindlist(gsea_display_selected, use.names = TRUE, fill = TRUE)
fwrite(selected_display, file.path(key_gene_gsea_work, "selected_GSEA_display_pathways.tsv"), sep = "\t")
gmt <- getGmt(file.path(client_root, "data", "reference", "h.all.v7.5.1.symbols.gmt"), geneIdType = SymbolIdentifier())
gs <- lapply(gmt, geneIds)
names(gs) <- vapply(gmt, setName, character(1))
gsva_es <- gsva(as.matrix(mat), gs)
cutoff <- 1
bulk_key_genes <- intersect(c("CCL20", "CLDN4", "CXCL8"), rownames(mat))
for (i in bulk_key_genes) {
    subexpr <- as.numeric(mat[i, ])
    names(subexpr) <- colnames(mat)
    lsam <- names(subexpr[subexpr < median(subexpr)])
    hsam <- names(subexpr[subexpr >= median(subexpr)])
    group_list <- data.frame(sample = c(lsam, hsam), group = c(rep("Lexp", length(lsam)), rep("Hexp", length(hsam))))
    design <- model.matrix(~0 + factor(group_list$group))
    colnames(design) <- levels(factor(group_list$group))
    rownames(design) <- colnames(gsva_es)
    contrast.matrix <- makeContrasts(Hexp - Lexp, levels = design)
    fit2 <- eBayes(contrasts.fit(lmFit(gsva_es[, group_list$sample], design), contrast.matrix))
    x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
    pathway <- str_replace(row.names(x), "HALLMARK_", "")
    df <- data.frame(ID = pathway, score = x$t)
    df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf), labels = c(1, 2, 3))
    sortdf <- df[order(df$score), ]
    sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
    write.table(cbind(pathway = rownames(x), x), file.path(work_dir, paste0(i, "_bulk_hallmark_GSVA.tsv")), sep = "\t", quote = FALSE, 
        row.names = FALSE)
    p <- ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = c("palegreen3", 
        "snow3", "dodgerblue4"), guide = FALSE) + geom_hline(yintercept = c(-cutoff, cutoff), color = "white", linetype = 2, 
        size = 0.3) + geom_text(data = subset(df, score < 0), aes(x = ID, y = 0, label = paste0(" ", ID), color = group), 
        size = 3, hjust = 0) + geom_text(data = subset(df, score > 0), aes(x = ID, y = -0.05, label = ID, color = group), 
        size = 3, hjust = 1) + scale_colour_manual(values = c("black", "snow3", "black"), guide = FALSE) + xlab("") + ylab(paste0("t value of GSVA score\n Hexp vs Lexp group of ", 
        i)) + theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(size = 0.6)) + theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
    ggsave(file.path(result_dir, paste0("GSVA_", i, ".pdf")), plot = p, width = 8, height = 7)
}
writeLines(capture.output(sessionInfo()), file.path(root, "logs/10_key_gene_enrichment_sessionInfo.txt"))
