source(file.path("R", "utils.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
    library(igraph)
})
set.seed(20260803)
root <- client_root
work_dir <- file.path(root, "work/09_gutmgene_key_gene_network")
result_dir <- file.path(root, "results/09_GutMGene_Key_Gene_Network")
data_dir <- file.path(root, "data/gutMGene_v2_20260802")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(result_dir, full.names = TRUE), recursive = TRUE)
rank_lookup <- c(p = "phylum", c = "class", o = "order", f = "family", g = "genus", s = "species")
normalize_taxon_name <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("_", " ", x, fixed = TRUE)
    trimws(gsub("[^a-z0-9]+", " ", x))
}
auc <- fread(file.path(root, "work/04_scBPS_key_cell/scBPS/key_cell_BPS_AUC.txt"))
moderate <- auc[class == "Moderate"]
stopifnot(nrow(moderate) == 94L)
moderate[, `:=`(rank_code, sub("_.*$", "", trait))]
moderate[, `:=`(taxon_rank, unname(rank_lookup[rank_code]))]
moderate[, `:=`(scbps_taxon_name, sub("^[kpcofgs]_+", "", trait))]
moderate[, `:=`(scbps_taxon_norm, normalize_taxon_name(scbps_taxon_name))]
stopifnot(!anyNA(moderate$taxon_rank))
taxonomy_alias <- data.table(scbps_taxon_norm = "bacteroidetes", taxon_rank = "phylum", gutmgene_taxon_norm = "bacteroidota", 
    expected_taxid = "976")
moderate[, `:=`(gutmgene_taxon_norm, scbps_taxon_norm)]
moderate[taxonomy_alias, on = .(scbps_taxon_norm, taxon_rank), `:=`(gutmgene_taxon_norm, i.gutmgene_taxon_norm)]
mg <- fread(file.path(data_dir, "microbe_gene.tsv"))
mg <- mg[tolower(humanormouse) == "human"]
mg[, `:=`(taxon_rank, tolower(rank))]
mg[, `:=`(gutmgene_taxon_norm, normalize_taxon_name(gut_microbiota))]
mg[, `:=`(gut_microbiota_ncbi_id, as.character(gut_microbiota_ncbi_id))]
taxonomy_matches <- merge(moderate, unique(mg[, .(gutmgene_taxon_norm, taxon_rank, gut_microbiota, gut_microbiota_ncbi_id)]), 
    by = c("gutmgene_taxon_norm", "taxon_rank"), allow.cartesian = TRUE)
stopifnot(uniqueN(taxonomy_matches$trait) == 8L)
alias_check <- taxonomy_matches[trait == "p_Bacteroidetes", unique(gut_microbiota_ncbi_id)]
stopifnot(identical(alias_check, "976"))
direct_candidates <- merge(moderate, mg, by = c("gutmgene_taxon_norm", "taxon_rank"), allow.cartesian = TRUE, suffixes = c("_scbps", 
    "_gutmgene"))
direct_candidates <- unique(direct_candidates)
candidate_genes <- sort(unique(direct_candidates$gene))
stopifnot(uniqueN(direct_candidates$trait) == 8L)
stopifnot(length(candidate_genes) == 20L)
deg <- fread(file.path(root, "work/05_key_cell_differential_gsea/IBD_vs_HC_pseudobulk_edgeR.tsv"))
deg_sig <- deg[PValue < 0.05 & abs(logFC) >= 0.5]
direct <- direct_candidates[gene %in% deg_sig$gene]
direct <- merge(direct, deg_sig[, .(gene, logFC, logCPM, F, PValue, gene_FDR = FDR)], by = "gene", all.x = TRUE)
direct <- unique(direct[, .(scbps_trait = trait, scbps_AUC = AUC, scbps_FDR = FDR, scbps_class = class, taxon_rank, gut_microbiota, 
    gut_microbiota_ncbi_id, gene, gene_logFC = logFC, gene_logCPM = logCPM, gene_F = F, gene_PValue = PValue, gene_FDR, alteration, 
    evidence = associative_modes, pmid)])
direct <- direct[order(gene_PValue, -abs(gene_logFC), gut_microbiota)]
stopifnot(nrow(direct) == 6L)
key_genes <- sort(unique(direct$gene))
stopifnot(identical(key_genes, sort(c("CCL20", "CLDN4", "CXCL8"))))
write.table(moderate, file.path(work_dir, "moderate_microbes_key_cell.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(taxonomy_matches, file.path(work_dir, "scbps_gutmgene_taxonomy_matches.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(KeyGene = candidate_genes), file.path(work_dir, "gutmgene_direct_candidate_genes.tsv"), sep = "\t", 
    quote = FALSE, row.names = FALSE)
write.table(direct, file.path(work_dir, "gutmgene_direct_microbe_gene_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(KeyGene = key_genes), file.path(work_dir, "key_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg[gene %in% key_genes], file.path(work_dir, "key_gene_pseudobulk_statistics.tsv"), sep = "\t", quote = FALSE, 
    row.names = FALSE)
edges <- unique(direct[, .(from = gut_microbiota, to = gene, relation = "direct_microbe_gene", alteration, evidence, pmid)])
nodes <- data.table(id = unique(c(edges$from, edges$to)))
nodes[, `:=`(type, ifelse(id %in% key_genes, "Gene", "Microbe"))]
nodes[, `:=`(label, id)]
write.table(edges, file.path(work_dir, "edges.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
write.table(nodes, file.path(work_dir, "nodes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
pdf(file.path(result_dir, "1.GutMGene_Key_Cell_DEG_Microbe_Network.pdf"), width = 10, height = 10)
g <- graph_from_data_frame(d = edges[, .(from, to)], vertices = nodes, directed = FALSE)
V(g)$type <- nodes$type[match(V(g)$name, nodes$id)]
V(g)$color <- c(Microbe = "#4DBBD5", Gene = "#E64B35")[V(g)$type]
V(g)$size <- c(Microbe = 10, Gene = 16)[V(g)$type]
V(g)$frame.color <- "grey30"
V(g)$label.color <- "black"
V(g)$label.cex <- 0.85
V(g)$label.dist <- 0.5
E(g)$color <- "grey65"
E(g)$width <- 1.4
layout <- layout_with_fr(g, niter = 2000)
plot(g, layout = layout, vertex.label.family = "sans", edge.curved = 0.08, margin = 0.15, main = "GutMGene Direct Microbe-Key Gene Network")
legend("topleft", legend = c("Enterocytes-associated microbe", "IBD vs HC key gene"), pch = 21, pt.bg = c("#4DBBD5", "#E64B35"), 
    pt.cex = c(1.8, 2.2), bty = "n", cex = 0.9)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(root, "logs/09_gutmgene_key_gene_network_sessionInfo.txt"))
