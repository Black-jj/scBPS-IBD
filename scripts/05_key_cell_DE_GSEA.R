source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(edgeR)
    library(ggplot2)
    library(ggrepel)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(cowplot)
    library(GSEABase)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/05_key_cell_differential_gsea")
result_dir <- file.path(root, "results/05_Key_Cell_DE_GSEA")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
scRNA <- readRDS(file.path(root, "work/02_scRNA_annotation_abundance/scRNA_annotated.rds"))
key_tab <- read.table(file.path(root, "work/04_scBPS_key_cell/scBPS/key_cell.txt"), header = TRUE, check.names = FALSE)
key_cell <- as.character(key_tab[[1]][1])
subRNA <- subset(scRNA, subset = cellType_1 == key_cell)
meta <- unique(data.frame(sample = subRNA$orig.ident, group3 = subRNA$group3))
stopifnot(!anyDuplicated(meta$sample), all(c("HC", "UC", "CD") %in% meta$group3))
counts <- GetAssayData(subRNA, assay = "RNA", slot = "counts")
pb <- do.call(cbind, lapply(meta$sample, function(z) Matrix::rowSums(counts[, subRNA$orig.ident == z, drop = FALSE])))
colnames(pb) <- meta$sample
rownames(pb) <- rownames(counts)
saveRDS(pb, file.path(work_dir, "key_cell_pseudobulk_counts.rds"), compress = FALSE)
write.table(meta, file.path(work_dir, "sample_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
run_ql <- function(groups, contrast_name) {
    keep_samples <- if (contrast_name == "IBD_vs_HC") 
        rep(TRUE, nrow(meta))
    else meta$group3 %in% strsplit(contrast_name, "_vs_")[[1]]
    m <- meta[keep_samples, , drop = FALSE]
    x <- pb[, m$sample, drop = FALSE]
    g <- if (contrast_name == "IBD_vs_HC") 
        factor(ifelse(m$group3 == "HC", "HC", "IBD"), levels = c("HC", "IBD"))
    else factor(m$group3, levels = rev(strsplit(contrast_name, "_vs_")[[1]]))
    y <- DGEList(x, group = g)
    design <- model.matrix(~0 + g)
    keep <- filterByExpr(y, design)
    y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    fit <- glmQLFit(estimateDisp(y, design, robust = TRUE), design, robust = TRUE)
    qlf <- glmQLFTest(fit, contrast = c(-1, 1))
    tab <- topTags(qlf, n = Inf, sort.by = "none")$table
    tab$gene <- rownames(tab)
    tab
}
plot_volcano_shared <- function(tab, title, file) {
    PbmcMarkers <- data.frame(avg_log2FC = tab$logFC, p_val_adj = pmax(tab$FDR, .Machine$double.xmin), gene = tab$gene)
    PbmcMarkers <- PbmcMarkers[-log10(PbmcMarkers$p_val_adj) > 0, , drop = FALSE]
    PbmcMarkers$lab <- ""
    ord <- order(PbmcMarkers$avg_log2FC)
    label_i <- unique(c(head(ord, 10), tail(ord, 10)))
    PbmcMarkers$lab[label_i] <- PbmcMarkers$gene[label_i]
    PbmcMarkers$logP <- -log10(PbmcMarkers$p_val_adj)
    logFCfilter <- 0.5
    p_valFilter <- 0.05
    p5 <- ggplot(PbmcMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC)) + geom_point(aes(size = logP), 
        alpha = 0.9) + scale_color_gradientn(colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"), name = "logFC", 
        values = seq(0, 1, 0.2)) + scale_fill_gradientn(colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"), 
        name = "logFC", values = seq(0, 1, 0.2)) + scale_size("-log10(p)") + geom_text_repel(data = PbmcMarkers, aes(x = avg_log2FC, 
        y = -log10(p_val_adj), label = lab), size = 4, box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"), 
        segment.color = "black", show.legend = FALSE) + theme_bw() + ylab("-log10 (FDR)") + xlab("log2FC") + geom_vline(xintercept = c(-logFCfilter, 
        logFCfilter), lty = 3, col = "black", lwd = 0.5) + geom_hline(yintercept = -log10(p_valFilter), lty = 3, col = "black", 
        lwd = 0.5) + labs(title = title)
    ggsave(file, plot = p5, width = 7, height = 6)
}
collapse_rank <- function(tab) {
    ids <- bitr(tab$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    z <- merge(tab[, c("gene", "logFC")], ids, by.x = "gene", by.y = "SYMBOL")
    z <- z[order(abs(z$logFC), decreasing = TRUE), ]
    z <- z[!duplicated(z$ENTREZID), ]
    out <- z$logFC
    names(out) <- z$ENTREZID
    sort(out, decreasing = TRUE)
}
plot_gsea_shared <- function(x, file) {
    x@result$Description <- sub("^HALLMARK_", "", x@result$Description)
    tab <- as.data.frame(x)
    tab <- tab[order(tab$p.adjust, -abs(tab$NES)), , drop = FALSE]
    geneSetID <- head(tab$ID, 3)
    mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9", "#7A142C", "#5D90BA", 
        "#431A3D", "#91612D", "#6E568C")
    gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
    gsdata$gsym <- rep(names(x@geneList), length(geneSetID))
    selectedGeneID <- names(sort(abs(x@geneList), decreasing = TRUE))[1]
    p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) + 
        scale_color_manual(values = mycol) + geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + ylab("Enrichment\n Score") + 
        theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "top", legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent")) + theme(axis.text.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), plot.margin = margin(t = 0.2, 
            r = 0.2, b = 0, l = 0.2, unit = "cm"))
    p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + 
        ylab(NULL) + scale_color_manual(values = mycol) + theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "none", 
        plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), axis.text = element_blank(), axis.line.x = element_blank()) + 
        scale_y_continuous(expand = c(0, 0))
    df2 <- p.res$data
    df2$y <- p.res$data$geneList[df2$x]
    selectgenes <- merge(data.frame(gsym = selectedGeneID), df2, by = "gsym")
    selectgenes <- selectgenes[selectgenes$position == 1, , drop = FALSE]
    p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + geom_segment(data = df2, 
        aes_(x = ~x, xend = ~x, y = ~y, yend = 0), color = "grey") + geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values = mycol, 
        guide = FALSE) + scale_color_manual(values = mycol, guide = FALSE) + geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + 
        ylab("Ranked list\n metric") + xlab("Rank in ordered dataset") + theme_bw() + theme(axis.text.y = element_text(size = 12, 
        face = "bold"), panel.grid = element_blank()) + geom_text_repel(data = selectgenes, show.legend = FALSE, direction = "x", 
        ylim = c(2, NA), angle = 90, size = 2.5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) + 
        theme(plot.margin = margin(t = -0.1, r = 0.4, b = 0.2, l = 0.2, unit = "cm"))
    p <- plot_grid(plotlist = list(p.res, p2, p.pos), ncol = 1, align = "v", rel_heights = c(1.5, 0.5, 1.5))
    ggsave(file, width = 8, height = 6, plot = p)
}
hallmark <- read.gmt(file.path(client_root, "data", "reference", "h.all.v7.5.1.symbols.gmt"))
kegg_reference <- readRDS(file.path(client_root, "data", "reference", "kegg.all.entrez.hsa.rds"))
kegg_term2gene <- unique(kegg_reference[, c("PathwayID", "EntrezID")])
kegg_term2name <- unique(kegg_reference[, c("PathwayID", "PathwayName")])
gseKEGG <- function(geneList, organism = "hsa", ...) {
    geneList <- geneList[!is.na(geneList) & !duplicated(names(geneList))]
    clusterProfiler::GSEA(geneList = geneList, TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name, pvalueCutoff = 1, 
        verbose = FALSE)
}
contrasts <- c("IBD_vs_HC", "UC_vs_HC", "CD_vs_HC")
all_res <- list()
for (z in contrasts) {
    tab <- run_ql(NULL, z)
    all_res[[z]] <- tab
    write.table(tab, file.path(work_dir, paste0(z, "_pseudobulk_edgeR.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    plot_volcano_shared(tab, paste0(key_cell, " pseudobulk: ", gsub("_", " ", z)), file.path(result_dir, paste0(z, "_1.Volcano.pdf")))
    rank_id <- collapse_rank(tab)
    kk <- gseKEGG(rank_id, organism = "hsa", pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, eps = 0, verbose = FALSE)
    go <- gseGO(rank_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, 
        eps = 0, verbose = FALSE)
    sym_rank <- tab$logFC
    names(sym_rank) <- tab$gene
    sym_rank <- sort(sym_rank[!duplicated(names(sym_rank))], decreasing = TRUE)
    hm <- GSEA(sym_rank, TERM2GENE = hallmark[, c("term", "gene")], pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, 
        eps = 0, verbose = FALSE)
    write.table(as.data.frame(kk), file.path(work_dir, paste0(z, "_KEGG_GSEA.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(go), file.path(work_dir, paste0(z, "_GO_BP_GSEA.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(hm), file.path(work_dir, paste0(z, "_HALLMARK_GSEA.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    plot_gsea_shared(kk, file.path(result_dir, paste0(z, "_2.KEGG_GSEA.pdf")))
    plot_gsea_shared(go, file.path(result_dir, paste0(z, "_3.GO_BP_GSEA.pdf")))
    plot_gsea_shared(hm, file.path(result_dir, paste0(z, "_4.HALLMARK_GSEA.pdf")))
    Idents(subRNA) <- "group3"
    if (z == "IBD_vs_HC") {
        subRNA$IBD_group <- ifelse(subRNA$group3 == "HC", "HC", "IBD")
        Idents(subRNA) <- "IBD_group"
        cell_de <- FindMarkers(subRNA, ident.1 = "IBD", ident.2 = "HC", assay = "RNA", logfc.threshold = 0, min.pct = 0.1)
    }
    else {
        v <- strsplit(z, "_vs_")[[1]]
        cell_de <- FindMarkers(subRNA, ident.1 = v[1], ident.2 = v[2], assay = "RNA", logfc.threshold = 0, min.pct = 0.1)
    }
    write.table(cbind(gene = rownames(cell_de), cell_de), file.path(work_dir, paste0(z, "_cell_level_descriptive.tsv")), 
        sep = "\t", quote = FALSE, row.names = FALSE)
}
writeLines(key_cell, file.path(work_dir, "key_cell_type.txt"))
writeLines(capture.output(sessionInfo()), file.path(root, "logs/05_key_cell_differential_gsea_sessionInfo.txt"))
