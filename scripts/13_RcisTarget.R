source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(RcisTarget)
    library(data.table)
    library(igraph)
    library(dplyr)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/13_rcistarget")
result_dir <- file.path(root, "results/13_RcisTarget")
final_result_dir <- file.path(root, "results/14.Transcription factor")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(work_dir, full.names = TRUE), recursive = TRUE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
load(file.path(root, "data", "reference", "RcisTarget", "hg19_motifAnnotation_cisbpOnly.RData"))
load(file.path(root, "data", "reference", "RcisTarget", "hg19_500bpUpstream_motifRanking_cispbOnly.RData"))
key_genes <- fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene
stopifnot(setequal(key_genes, c("CCL20", "CLDN4", "CXCL8")))
motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
ranking_table <- getRanking(motifRankings)
gene_map <- data.table(KeyGene = c("CCL20", "CLDN4", "CXCL8"), DatabaseGene = c("CCL20", "CLDN4", "IL8"))
stopifnot(all(gene_map$DatabaseGene %in% colnames(ranking_table)))
fwrite(gene_map, file.path(work_dir, "key_gene_database_mapping.tsv"), sep = "\t")
motif_annotation <- unique(as.data.table(hg19_motifAnnotation_cisbpOnly)[directAnnotation == TRUE, .(motif, TF, annotationSource)])
all_edges <- rbindlist(lapply(seq_len(nrow(gene_map)), function(i) {
    key_gene <- gene_map$KeyGene[i]
    database_gene <- gene_map$DatabaseGene[i]
    motif_rank <- data.table(motif = ranking_table$features, rank = as.numeric(ranking_table[[database_gene]]))
    candidate <- merge(motif_rank[is.finite(rank) & rank <= 500], motif_annotation, by = "motif", allow.cartesian = TRUE)
    setorder(candidate, rank, motif, TF)
    candidate <- candidate[, head(.SD, 1), by = TF]
    candidate[, `:=`(KeyGene = key_gene, DatabaseGene = database_gene)]
    candidate
}), fill = TRUE)
setorder(all_edges, KeyGene, rank, TF)
fwrite(all_edges, file.path(work_dir, "TF_key_gene_all_direct_top500.tsv"), sep = "\t")
plot_edges <- all_edges[, head(.SD, 8), by = KeyGene]
result <- unique(plot_edges[, .(TF_highConf = TF, enrichedGenes = KeyGene, motif_rank = rank, motif, DatabaseGene)])
fwrite(result, file.path(work_dir, "TF_gene_edges.tsv"), sep = "\t")
stopifnot(setequal(unique(result$enrichedGenes), key_genes), nrow(result) == 24)
g <- graph_from_data_frame(d = result[, .(TF_highConf, enrichedGenes)], directed = TRUE)
V(g)$type <- ifelse(V(g)$name %in% result$TF_highConf, "TF", "KeyGene")
V(g)$color <- c(TF = "#4DBBD5", KeyGene = "#E64B35")[V(g)$type]
V(g)$size <- c(TF = 14, KeyGene = 24)[V(g)$type]
lay <- layout_with_kk(g, maxiter = 5000)
plot_network <- function(path) {
    pdf(path, width = 12, height = 8, useDingbats = FALSE)
    plot(g, layout = lay, vertex.label.cex = 0.85, vertex.label.color = "black", edge.color = "grey70", edge.arrow.size = 0.35, 
        main = "TF-Key Gene Regulatory Network")
    dev.off()
}
plot_network(file.path(result_dir, "1.TF-Gene-Network.pdf"))
plot_network(file.path(final_result_dir, "1.TF-gene network.pdf"))
writeLines(capture.output(sessionInfo()), file.path(root, "logs/13_rcistarget_sessionInfo.txt"))
