source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout")))
suppressPackageStartupMessages({
    library(data.table)
})
root <- client_root
stats_dir <- file.path(root, "work/20_geneformer/gf_stats")
work_dir <- file.path(root, "work/20_geneformer/enrichment")
result_dir <- file.path(root, "results/21_Geneformer_Downstream_Enrichment")
kegg_snapshot <- file.path(client_root, "data", "reference", "kegg.all.entrez.hsa.rds")
ko_genes <- c("CCL20", "CLDN4", "CXCL8")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
stopifnot(file.exists(kegg_snapshot))
kegg_reference <- readRDS(kegg_snapshot)
stopifnot(all(c("EntrezID", "PathwayID", "PathwayName") %in% colnames(kegg_reference)))
kegg_term2gene <- unique(kegg_reference[, c("PathwayID", "EntrezID")])
kegg_term2name <- unique(kegg_reference[, c("PathwayID", "PathwayName")])
valid_pdf <- function(path) {
    info <- tryCatch(system2("pdfinfo", path, stdout = TRUE, stderr = TRUE), error = function(e) character())
    page_line <- grep("^Pages:", info, value = TRUE)
    length(page_line) == 1L && as.integer(sub("^Pages:[[:space:]]*", "", page_line)) > 0L
}
write_status_pdf <- function(path, ko_gene, message) {
    pdf(path, width = 9, height = 6, useDingbats = FALSE)
    plot.new()
    title(main = paste0("Geneformer KO ", ko_gene, ": GO/KEGG enrichment"))
    text(0.5, 0.52, message, cex = 1.1)
    text(0.5, 0.42, "Filter: FDR < 0.05 and Cosine_Shift > 0", cex = 0.9)
    dev.off()
}
status <- rbindlist(lapply(ko_genes, function(ko_gene) {
    stat_file <- file.path(stats_dir, paste0(ko_gene, "_affected_stats.csv"))
    affected <- fread(stat_file)
    required <- c("Affected_gene_name", "FDR", "Cosine_Shift")
    stopifnot(all(required %in% colnames(affected)))
    sig <- unique(affected[!is.na(Affected_gene_name) & Affected_gene_name != "" & FDR < 0.05 & Cosine_Shift > 0])
    ko_work <- file.path(work_dir, ko_gene)
    dir.create(ko_work, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(ko_work, full.names = TRUE), recursive = TRUE)
    fwrite(sig, file.path(ko_work, paste0(ko_gene, "_significant_affected_genes.tsv")), sep = "\t")
    genes <- sort(unique(sig$Affected_gene_name))
    if (!length(genes)) {
        out_pdf <- file.path(result_dir, paste0("KO_", ko_gene, "_No_Significant_Affected_Gene.pdf"))
        write_status_pdf(out_pdf, ko_gene, "No significant downstream affected gene for enrichment")
        return(data.table(KO_gene = ko_gene, affected_tested = nrow(affected), significant_affected = 0L, mapped_entrez = 0L, 
            GO_terms = 0L, KEGG_terms = 0L, output_pdfs = 1L))
    }
    fwrite(data.table(gene = genes), file.path(ko_work, "gene.txt"), sep = "\t")
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(ko_work)
    run_env <- new.env(parent = environment())
    run_env$Options <- options
    run_env$enrichKEGG <- function(gene, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, ...) {
        clusterProfiler::enricher(gene = gene, TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name, pvalueCutoff = pvalueCutoff, 
            qvalueCutoff = qvalueCutoff, ...)
    }
    evalq({
        library("org.Hs.eg.db")
        rt = read.table("gene.txt", sep = "\t", check.names = F, header = T)
        genes = as.vector(rt[, 1])
        entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
        entrezIDs <- as.character(entrezIDs)
        out = cbind(rt, entrezID = entrezIDs)
        write.table(out, file = "id.txt", sep = "\t", quote = F, row.names = F)
        library("clusterProfiler")
        library("org.Hs.eg.db")
        library("enrichplot")
        library("ggplot2")
        rt = read.table("id.txt", sep = "\t", header = T, check.names = F)
        rt = rt[is.na(rt[, "entrezID"]) == F, ]
        gene = rt$entrezID
        kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", readable = T)
        write.table(kk, file = "GO.txt", sep = "\t", quote = F, row.names = F)
        if (nrow(as.data.frame(kk)) > 0) {
            pdf(file = "GO.pdf", width = 10, height = 8)
            print(barplot(kk, drop = TRUE, showCategory = 10, split = "ONTOLOGY", label_format = 100) + facet_grid(ONTOLOGY ~ 
                ., scale = "free"))
            dev.off()
        }
        if (nrow(as.data.frame(kk)) > 0) {
            pdf(file = "GObubble.pdf", width = 10, height = 8)
            print(dotplot(kk, showCategory = 10, split = "ONTOLOGY", label_format = 100) + facet_grid(ONTOLOGY ~ ., scale = "free"))
            dev.off()
        }
        library("clusterProfiler")
        library("org.Hs.eg.db")
        library("enrichplot")
        library("ggplot2")
        library(R.utils)
        options(clusterProfiler.download.method = "curl")
        Options(download.file.extra = "-k")
        R.utils::setOption("clusterProfiler.download.method", "auto")
        rt = read.table("id.txt", sep = "\t", header = T, check.names = F)
        rt = rt[is.na(rt[, "entrezID"]) == F, ]
        gene = rt$entrezID
        kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
        write.table(kk, file = "KEGG.txt", sep = "\t", quote = F, row.names = F)
        if (nrow(as.data.frame(kk)) > 0) {
            pdf(file = "KEGG.pdf", width = 10, height = 7)
            print(barplot(kk, drop = TRUE, showCategory = 30, label_format = 100))
            dev.off()
        }
        if (nrow(as.data.frame(kk)) > 0) {
            pdf(file = "KEGGbubble.pdf", width = 10, height = 7)
            print(dotplot(kk, showCategory = 30, label_format = 100))
            dev.off()
        }
    }, envir = run_env)
    setwd(old_wd)
    id_file <- file.path(ko_work, "id.txt")
    mapped_entrez <- if (file.exists(id_file)) {
        sum(!is.na(fread(id_file)$entrezID))
    }
    else 0L
    go_terms <- if (file.exists(file.path(ko_work, "GO.txt"))) 
        nrow(fread(file.path(ko_work, "GO.txt")))
    else 0L
    kegg_terms <- if (file.exists(file.path(ko_work, "KEGG.txt"))) 
        nrow(fread(file.path(ko_work, "KEGG.txt")))
    else 0L
    pdfs <- list.files(ko_work, pattern = "[.]pdf$", full.names = TRUE)
    if (!length(pdfs)) {
        out_pdf <- file.path(result_dir, paste0("KO_", ko_gene, "_No_Significant_GO_KEGG.pdf"))
        write_status_pdf(out_pdf, ko_gene, paste0("Significant affected genes: ", length(genes), "; no GO/KEGG term passed p/q < 0.05"))
        pdfs <- out_pdf
    }
    else {
        copied <- file.path(result_dir, paste0("KO_", ko_gene, "_", basename(pdfs)))
        stopifnot(all(file.copy(pdfs, copied, overwrite = TRUE)))
        pdfs <- copied
    }
    stopifnot(all(vapply(pdfs, valid_pdf, logical(1))))
    data.table(KO_gene = ko_gene, affected_tested = nrow(affected), significant_affected = length(genes), mapped_entrez = mapped_entrez, 
        GO_terms = go_terms, KEGG_terms = kegg_terms, output_pdfs = length(pdfs))
}), fill = TRUE)
fwrite(status, file.path(work_dir, "per_KO_enrichment_status.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "logs/20_geneformer_enrichment_sessionInfo.txt"))
