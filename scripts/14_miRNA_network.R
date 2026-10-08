source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(httr)
    library(igraph)
    library(dplyr)
})
set.seed(20260802)
root <- client_root
work_dir <- file.path(root, "work/14_mirna_network")
result_dir <- file.path(root, "results/14_miRWalk_miRNA_Network")
data_dir <- file.path(root, "data/miRWalk_20260802")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
KeyGene <- fread(file.path(root, "work/09_gutmgene_key_gene_network/key_genes.tsv"))$KeyGene
idmap <- AnnotationDbi::select(org.Hs.eg.db, keys = KeyGene, keytype = "SYMBOL", columns = "ENTREZID")
idmap <- unique(idmap[!is.na(idmap$ENTREZID), ])
stopifnot(nrow(idmap) > 0)
all_data <- list()
for (i in seq_len(nrow(idmap))) {
    gene <- idmap$SYMBOL[i]
    gene_id <- idmap$ENTREZID[i]
    dest <- file.path(data_dir, paste0(gene, "_", gene_id, ".csv"))
    if (!file.exists(dest)) {
        h <- handle("http://mirwalk.umm.uni-heidelberg.de")
        r1 <- GET(paste0("http://mirwalk.umm.uni-heidelberg.de/human/gene/", gene_id, "/"), handle = h, timeout(120))
        stop_for_status(r1)
        r2 <- GET("http://mirwalk.umm.uni-heidelberg.de/export/csv/", handle = h, write_disk(dest, overwrite = TRUE), timeout(300))
        stop_for_status(r2)
    }
    x <- fread(dest)
    x$query_gene <- gene
    all_data[[length(all_data) + 1]] <- x
}
mirwalk <- rbindlist(all_data, fill = TRUE)
validated_support <- !is.na(mirwalk$validated) & mirwalk$validated != ""
other_support <- (!is.na(mirwalk$TargetScan) & mirwalk$TargetScan == 1) | (!is.na(mirwalk$miRDB) & mirwalk$miRDB == 1)
key.mirwalk <- mirwalk[bindingp >= 0.95 & (validated_support | other_support)]
stopifnot(nrow(key.mirwalk) > 0)
network.table <- unique(data.frame(Gene = key.mirwalk$query_gene, microrna = key.mirwalk$mirnaid))
write.table(key.mirwalk, file.path(work_dir, "miRWalk_filtered_interactions.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
write.table(network.table, file.path(work_dir, "edges.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
nodes <- data.frame(id = unique(c(network.table$Gene, network.table$microrna)))
nodes$type <- ifelse(nodes$id %in% network.table$microrna, "miRNA", "Gene")
write.table(nodes, file.path(work_dir, "nodes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
pdf(file.path(result_dir, "1.miRWalk-miRNA-Gene-Network.pdf"), width = 7, height = 7)
g <- graph_from_data_frame(d = network.table, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% network.table$microrna, "miRNA", "Gene")
V(g)$color <- ifelse(V(g)$type == "miRNA", "#4DBBD5", "#E64B35")
V(g)$size <- ifelse(V(g)$type == "miRNA", 5, 10)
lay <- layout_with_fr(g)
plot(g, layout = lay, vertex.label.cex = 0.5, vertex.label.color = "black", edge.color = "grey70", main = "miRWalk miRNA Regulatory Network")
dev.off()
writeLines(capture.output(sessionInfo()), file.path(root, "logs/14_mirna_network_sessionInfo.txt"))
