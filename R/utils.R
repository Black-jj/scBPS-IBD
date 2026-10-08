options(stringsAsFactors = FALSE)
client_root <- normalizePath(getwd(), mustWork = TRUE)
require(Seurat)
require(ggplot2)
require(harmony)
require(dplyr)
require(reshape2)
require(scRNAtoolVis)
require(ggpubr)
options(dplyr.print_max = 1e+09)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)
ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
    filenames <- paste0(filename, plot_Device)
    for (filename_tmp in filenames) {
        ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
    }
}
paramSweep_V4 <- function(seu, PCs = 1:10, sct = FALSE, num.cores = 1) {
    require(Seurat)
    require(fields)
    require(parallel)
    pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
    pN <- seq(0.05, 0.3, by = 0.05)
    min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
    pK.test <- round(pK * min.cells)
    pK <- pK[which(pK.test >= 1)]
    orig.commands <- seu@commands
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
    if (num.cores > 1) {
        require(parallel)
        cl <- makeCluster(num.cores)
        output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep, n.real.cells, real.cells, pK, pN, data, orig.commands, 
            PCs, sct, mc.cores = num.cores)
        stopCluster(cl)
    }
    else {
        output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep, n.real.cells, real.cells, pK, pN, data, orig.commands, 
            PCs, sct)
    }
    sweep.res.list <- list()
    list.ind <- 0
    for (i in 1:length(output2)) {
        for (j in 1:length(output2[[i]])) {
            list.ind <- list.ind + 1
            sweep.res.list[[list.ind]] <- output2[[i]][[j]]
        }
    }
    name.vec <- NULL
    for (j in 1:length(pN)) {
        name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_"))
    }
    names(sweep.res.list) <- name.vec
    return(sweep.res.list)
}
doubletFinder_V4 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
    require(Seurat)
    require(fields)
    require(KernSmooth)
    if (reuse.pANN != FALSE) {
        pANN.old <- seu@meta.data[, reuse.pANN]
        classifications <- rep("Singlet", length(pANN.old))
        classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
        return(seu)
    }
    if (reuse.pANN == FALSE) {
        real.cells <- rownames(seu@meta.data)
        data <- seu@assays$RNA@counts[, real.cells]
        n_real.cells <- length(real.cells)
        n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
        real.cells1 <- rep(real.cells, length.out = n_doublets)
        real.cells2 <- rep(rev(real.cells), length.out = n_doublets)
        doublets <- (data[, real.cells1] + data[, real.cells2])/2
        colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
        data_wdoublets <- cbind(data, doublets)
        if (!is.null(annotations)) {
            stopifnot(typeof(annotations) == "character")
            stopifnot(length(annotations) == length(Cells(seu)))
            stopifnot(!any(is.na(annotations)))
            annotations <- factor(annotations)
            names(annotations) <- Cells(seu)
            doublet_types1 <- annotations[real.cells1]
            doublet_types2 <- annotations[real.cells2]
        }
        orig.commands <- seu@commands
        if (sct == FALSE) {
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
            seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, margin = orig.commands$NormalizeData.RNA@params$margin)
            seu_wdoublets <- FindVariableFeatures(seu_wdoublets, selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
            seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, model.use = orig.commands$ScaleData.RNA$model.use, 
                do.scale = orig.commands$ScaleData.RNA$do.scale, do.center = orig.commands$ScaleData.RNA$do.center, scale.max = orig.commands$ScaleData.RNA$scale.max, 
                block.size = orig.commands$ScaleData.RNA$block.size, min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
            seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, verbose = FALSE)
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets)
            gc()
        }
        if (sct == TRUE) {
            require(sctransform)
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
            seu_wdoublets <- SCTransform(seu_wdoublets)
            seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets)
            gc()
        }
        dist.mat <- fields::rdist(pca.coord)
        pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
        if (!is.null(annotations)) {
            neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
        }
        rownames(pANN) <- real.cells
        colnames(pANN) <- "pANN"
        k <- round(nCells * pK)
        for (i in 1:n_real.cells) {
            neighbors <- order(dist.mat[, i])
            neighbors <- neighbors[2:(k + 1)]
            pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
            if (!is.null(annotations)) {
                for (ct in unique(annotations)) {
                  neighbors_that_are_doublets <- neighbors[neighbors > n_real.cells]
                  if (length(neighbors_that_are_doublets) > 0) {
                    neighbor_types[i, ] <- table(doublet_types1[neighbors_that_are_doublets - n_real.cells]) + table(doublet_types2[neighbors_that_are_doublets - 
                      n_real.cells])
                    neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, ])
                  }
                  else {
                    neighbor_types[i, ] <- NA
                  }
                }
            }
        }
        classifications <- rep("Singlet", n_real.cells)
        classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 1]
        seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
        if (!is.null(annotations)) {
            colnames(neighbor_types) <- levels(doublet_types1)
            for (ct in levels(doublet_types1)) {
                seu@meta.data[, paste("DF.doublet.contributors", pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, ct]
            }
        }
        return(seu)
    }
}
human2mouse <- function(x) {
    x1 <- substr(x, 1, 1)
    x2 <- substr(x, 2, nchar(x))
    return(paste0(x1, tolower(x2)))
}
getScatterplot <- function(object, gene1, gene2, cor.method = "pearson", jitter.num = 0.15, pos = TRUE) {
    if (gene1 %in% rownames(object)) {
        exp.mat <- GetAssayData(object = object, assay = "RNA") %>% .[c(gene1, gene2), ] %>% as.matrix() %>% t() %>% as.data.frame()
        if (pos) {
            if (nrow(exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]) > (nrow(exp.mat) * 0.01)) {
                exp.mat <- exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]
            }
            else {
                exp.mat <- exp.mat[which(exp.mat[, 1] > 0 | exp.mat[, 2] > 0), ]
            }
        }
        colnames(exp.mat) <- c("Var1", "Var2")
        plots <- ggplot(data = exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) + geom_smooth(method = "lm", se = T, 
            color = "red", size = 1) + stat_cor(method = cor.method) + labs(x = gene1, y = gene2) + geom_jitter(width = jitter.num, 
            height = jitter.num, color = "black", size = 1, alpha = 1) + theme_bw() + theme(panel.grid = element_blank(), 
            legend.text = element_text(colour = "black", size = 10), axis.text = element_text(colour = "black", size = 10), 
            axis.line = element_line(colour = "black"), panel.border = element_rect(size = 1, linetype = "solid", colour = "black"), 
            panel.background = element_rect(fill = "white"))
        return(plots)
    }
}
capitalize_first <- function(string) {
    words <- strsplit(string, " ")[[1]]
    capitalized_words <- toupper(substring(words, 1, 1))
    rest_of_words <- substring(words, 2)
    capitalized_string <- paste0(capitalized_words, rest_of_words, collapse = " ")
    return(capitalized_string)
}
scRNAAutoAnno <- function(Path = ".", SeuratObject = NA, Multi = TRUE, ref, labels, min.features = 200, Idents = "RNA_snn_res.0.2", 
    DoubletFinder = TRUE, species = NA, KeyCell = NA, KeyGene = NA, logFCfilter = 0.585, adjPvalFilter = 0.05, assay = "RNA", 
    set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 20, nfeatures = 2000, npcs = 50) {
    AnalysisIndex <- scan(sep = ",", quiet = TRUE)
    scRNA <- SeuratObject
    if (!is.na(species)) {
        mouse <- ifelse(species == "mouse", T, F)
    }
    else {
        mouse <- F
    }
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    if (mouse) {
        s.genes <- human2mouse(s.genes)
        g2m.genes <- human2mouse(g2m.genes)
    }
    n <- 0
    KeyCells <- KeyCell
    for (step in AnalysisIndex) {
        n <- n + 1
        if (step == 1) {
            dir <- paste0(Path, "/1.QCFilter")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            grep("^[M,m][T,t]-", rownames(scRNA), value = T)
            grep("^R[P,p][SL,sl]", rownames(scRNA), value = T)
            scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, pattern = "^[M,m][T,t]-")
            scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern = "^R[P,p][SL,sl]")
            nomt <- ifelse(max(scRNA[["percent.mt"]]) == 0, T, F)
            if (nomt) {
                MT_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
                if (mouse) {
                  MT_genes <- human2mouse(MT_genes)
                }
                MT_genes <- intersect(MT_genes, rownames(scRNA))
                scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, features = MT_genes)
            }
            mt_value <- stats::mad(scRNA$percent.mt)
            lower_mt <- median(scRNA$percent.mt) - 3 * mt_value
            upper_mt <- median(scRNA$percent.mt) + 3 * mt_value
            nF_value <- stats::mad(scRNA$nFeature_RNA)
            lower_nF <- median(scRNA$nFeature_RNA) - 3 * nF_value
            upper_nF <- median(scRNA$nFeature_RNA) + 3 * nF_value
            nC_value <- stats::mad(scRNA$nCount_RNA)
            lower_nC <- median(scRNA$nCount_RNA) - 3 * nC_value
            upper_nC <- median(scRNA$nCount_RNA) + 3 * nC_value
            nR_value <- stats::mad(scRNA$percent.ribo)
            lower_nR <- median(scRNA$percent.ribo) - 3 * nR_value
            upper_nR <- median(scRNA$percent.ribo) + 2 * nR_value
            scRNA <- subset(x = scRNA, subset = nFeature_RNA >= min.features & percent.mt <= upper_mt & percent.ribo <= upper_nR & 
                nFeature_RNA <= upper_nF & nCount_RNA <= upper_nC)
            dim(scRNA)
            if (DoubletFinder) {
                require(DoubletFinder)
                scRNA_list <- SplitObject(scRNA, split.by = "orig.ident")
                newscRMA_list <- list()
                for (i in 1:length(scRNA_list)) {
                  data <- scRNA_list[[i]]
                  data <- NormalizeData(data)
                  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
                  data <- ScaleData(data)
                  data <- RunPCA(data, npcs = 30)
                  data <- RunUMAP(data, dims = 1:20)
                  sweep.res.list <- paramSweep_V4(data, PCs = 1:20, sct = FALSE)
                  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
                  bcmvn <- find.pK(sweep.stats)
                  p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
                  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
                  Doubletrate <- ncol(data) * 8 * 1e-06
                  nExp_poi <- round(Doubletrate * ncol(data))
                  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
                  data <- doubletFinder_V4(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, 
                    sct = FALSE)
                  colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
                  newscRMA_list[[i]] <- data
                }
                saveRDS(newscRMA_list, paste0(dir, "/doubletFinder_V4.rds"))
                plot_list <- list()
                for (i in 1:length(newscRMA_list)) {
                  plot_list[[i]] <- DimPlot(newscRMA_list[[i]], group.by = "doublet_info") + ggtitle(newscRMA_list[[i]]$orig.ident[1]) + 
                    theme(plot.title = element_text(hjust = 0.5))
                }
                for (i in 1:length(plot_list)) {
                  plot_list[[i]] <- plot_list[[i]] + NoAxes() + NoLegend()
                }
                p <- CombinePlots(plot_list, legend = "right")
                ggsave_fun(filename = paste0(dir, "/1.DoubletFinder"), plot = p, width = 12, height = 10)
                if (Multi) {
                  scRNA <- merge(newscRMA_list[[1]], newscRMA_list[2:length(newscRMA_list)])
                }
                else {
                  scRNA <- newscRMA_list[[1]]
                }
                scRNA <- subset(scRNA, subset = doublet_info == "Singlet")
                p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", 
                  ncol = 4)
                ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 20, height = 4)
                plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) + NoLegend()
                plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) + NoLegend()
                plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) + NoLegend()
                plot1 + plot2 + plot3
                ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
            }
            else {
                p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", 
                  ncol = 4)
                ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 20, height = 4)
                plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) + NoLegend()
                plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) + NoLegend()
                plot3 <- FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.ribo", raster = F) + NoLegend()
                plot1 + plot2 + plot3
                ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
            }
        }
        if (step == 2) {
            dir <- paste0(Path, "/2.Cluster")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            DefaultAssay(scRNA) <- assay
            scRNA <- NormalizeData(object = scRNA)
            scRNA <- CellCycleScoring(object = scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
            scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
            top10 <- head(VariableFeatures(scRNA), 10)
            plot1 <- VariableFeaturePlot(scRNA)
            plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
            plot1 + plot2
            ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
            scRNA <- ScaleData(object = scRNA, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
            scRNA <- RunPCA(scRNA, verbose = T, npcs = npcs, features = VariableFeatures(object = scRNA))
            p <- ElbowPlot(object = scRNA, ndims = npcs)
            ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)
            if (Multi) {
                scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
                scRNA <- RunUMAP(scRNA, n.neighbors = 10, min.dist = 0.05, reduction = "harmony", dims = 1:PC, verbose = T)
                scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "harmony", verbose = T)
                scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
                p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
                ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
                p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
                ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
            }
            else {
                scRNA <- RunUMAP(scRNA, n.neighbors = 10, min.dist = 0.05, reduction = "pca", dims = 1:PC, verbose = T)
                scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "pca", verbose = T)
                scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
                p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
                ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
            }
            pdf(paste0(dir, "/7.data.merge.harmony.pdf"))
            p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = "orig.ident", raster = F) + NoLegend()
            print(p)
            merge.res <- sapply(set.resolutions, function(x) {
                p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = paste0(assay, "_snn_res.", x), 
                  raster = F) + NoLegend()
                print(p)
            })
            dev.off()
            saveRDS(scRNA, paste0(dir, "/data.merge.harmony.2000.rds"))
        }
        if (step == 3) {
            dir <- paste0(Path, "/.CellAnnotate")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            Cluster.dir <- list.files(Path, pattern = "^.Cluster$", all.files = TRUE)
            p <- DimPlot(scRNA, group.by = Idents, label = T, reduction = "umap", raster = F) + NoLegend()
            print(p)
            ggsave_fun(filename = paste0(Path, "/", Cluster.dir, "/5.Umap.Cluster"), plot = p)
            CELL <- SingleR::SingleR(test = as.matrix(scRNA@assays$RNA@data), ref = ref, labels = labels, clusters = scRNA@meta.data[[Idents]])
            cells <- CELL[, ncol(CELL):1]
            cells <- data.frame(cells)
            cell_id <- cells$pruned.labels
            names(cell_id) <- levels(scRNA@meta.data[[Idents]])
            Idents(scRNA) <- Idents
            if (anyNA(cell_id)) {
                scRNA <- RenameIdents(scRNA, cell_id)
                scRNA$cellType_1 <- Idents(scRNA)
                print(cell_id)
                p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
                print(p)
                exit <- scan(what = "character", sep = ",", quiet = TRUE)
                if (exit == "yes") {
                  cell_id[which(is.na(cell_id))] <- "unknow"
                }
                if (exit == "no") {
                  return(scRNA)
                }
            }
            if (!anyNA(cell_id)) {
                p <- DimPlot(scRNA, group.by = Idents, raster = F, label = T) + NoLegend()
                print(p)
                while (T) {
                  a <- scan(what = "character")
                  gene <- intersect(a, rownames(scRNA))
                  if (length(gene) > 0) {
                    p <- FeaturePlot(scRNA, features = gene, raster = F, label = T, cols = c("lightgrey", "red"))
                    print(p)
                  }
                  if (length(a) == 0) {
                    break
                  }
                }
                cell_id <- scan(what = "character", sep = ";", quiet = FALSE)
                cell_id <- gsub(".*:", "", cell_id)
                names(cell_id) <- levels(scRNA@meta.data[[Idents]])
                Idents(scRNA) <- Idents
                scRNA <- RenameIdents(scRNA, cell_id)
                scRNA$cellType_1 <- Idents(scRNA)
                print(table(scRNA$cellType_1))
                p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
                print(p)
                saveRDS(scRNA, paste0(Path, "/result.rds"))
                Idents(scRNA) <- "cellType_1"
                scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
                scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05, ]
                write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))
                features <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
                print((features))
                features <- features$gene
                RenameGene <- scan(what = "character", sep = ",", quiet = TRUE)
                if (RenameGene == "yes") {
                  features <- scan(what = "character", sep = ",", quiet = FALSE)
                }
                my_pal2 <- c("#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", 
                  "#9ebcda", "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A", 
                  "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda", 
                  "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A", "#D4477D", 
                  "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda", "#74a9cf", 
                  "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A", "#D4477D", "#D24B27", 
                  "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda", "#74a9cf", "#fbdf72", 
                  "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A")
                umap <- scRNA@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cellType = scRNA@meta.data$cellType_1)
                celltypepos <- umap %>% group_by(cellType) %>% summarise(umap_1 = median(UMAP_1), umap_2 = median(UMAP_2))
                p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) + 
                  ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType), fontface = "bold", 
                    data = celltypepos, box.padding = 0.5, show.legend = FALSE)
                print(p)
                ggsave_fun(filename = paste0(dir, "/1.Umap1.plot"), plot = p, width = 10, height = 7)
                p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 
                  10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, 
                  vjust = 0.5), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
                  "cm"), )
                print(p)
                if (length(unique(features)) >= 50) {
                  ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 20, height = 9)
                }
                if (length(unique(features)) < 50) {
                  ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 17, height = 9)
                }
                tryCatch({
                  cellRatio <- cellRatioPlot(object = scRNA, sample.name = "group", celltype.name = "cellType_1", flow.curve = 0.5, 
                    fill.col = my_pal2) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
                  print(cellRatio)
                  ggsave_fun(filename = paste0(dir, "/3.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 7, 
                    height = 7)
                  data_p <- scRNA@meta.data[, c("group", "cellType_1")]
                  p <- ggstatsplot::ggbarstats(data = data_p, x = group, y = cellType_1, palette = "Set3") + theme(axis.text.x = element_text(angle = 45, 
                    hjust = 1))
                  print(p)
                  leg <- length(unique(data_p$cellType_1))
                  ggsave_fun(filename = paste0(dir, "/4.Dysfunction"), plot = p, width = 2 + 1.5 * leg, height = 7)
                  EachSampleClusterDis <- lapply(unique(scRNA@meta.data$orig.ident), function(Sample) {
                    seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$orig.ident %in% Sample, ]))
                    seur@meta.data$cellType_1 %>% table() %>% data.frame() %>% magrittr::set_colnames(c("CellTypes", "Number")) %>% 
                      dplyr::mutate(Per = 100 * Number/sum(Number))
                  })
                  names(EachSampleClusterDis) <- unique(scRNA@meta.data$orig.ident)
                  EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), 
                    times = unlist(lapply(EachSampleClusterDis, nrow))))
                  p1 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + 
                    geom_col(position = "fill", width = 0.8) + scale_y_continuous(expand = c(0, 0)) + labs(y = "Percentage of cell") + 
                    scale_fill_manual(values = my_pal2, name = "") + theme(axis.text = element_text(color = "black"), panel.background = element_blank(), 
                    legend.title = element_blank(), axis.text.y = element_text(color = "black"), axis.line = element_line(color = "black"), 
                    axis.title.y = element_blank(), legend.position = "right", legend.direction = "vertical", legend.text = element_text(size = 10)) + 
                    coord_flip()
                  EachSampleClusterDis <- lapply(unique(scRNA@meta.data$group), function(Sample) {
                    seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$group %in% Sample, ]))
                    seur@meta.data$cellType_1 %>% table() %>% data.frame() %>% magrittr::set_colnames(c("CellTypes", "Number")) %>% 
                      dplyr::mutate(Per = 100 * Number/sum(Number))
                  })
                  my_pal1 <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#6387C5")
                  names(EachSampleClusterDis) <- unique(scRNA@meta.data$group)
                  EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), 
                    times = unlist(lapply(EachSampleClusterDis, nrow))))
                  p2 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + 
                    geom_col(position = "fill", width = 0.8) + scale_y_continuous(expand = c(0, 0)) + labs(y = "Percentage of cell") + 
                    scale_fill_manual(values = my_pal1, name = "") + theme(axis.text = element_text(color = "black", ), panel.background = element_blank(), 
                    panel.grid = element_blank(), legend.title = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), 
                    axis.ticks.y.left = element_blank(), axis.title.y = element_blank(), axis.line.x = element_line(color = "black"), 
                    legend.position = "right", legend.direction = "vertical", legend.text = element_text(size = 10)) + coord_flip()
                  cell_counts <- table(scRNA@meta.data$cellType_1)
                  cell_counts_df <- data.frame(CellType = names(cell_counts), Count = as.numeric(cell_counts/1000))
                  cell_counts_df <- cell_counts_df[order(cell_counts_df$Count, decreasing = T), ]
                  p3 <- ggplot(cell_counts_df, aes(x = reorder(CellType, Count), y = Count)) + scale_y_continuous(expand = c(0, 
                    0)) + labs(y = "Number of cell(10^3)") + geom_bar(stat = "identity", fill = "#6387C5") + theme(axis.text = element_text(color = "black"), 
                    panel.background = element_blank(), panel.grid = element_blank(), legend.title = element_blank(), axis.text.y = element_blank(), 
                    axis.ticks.y.left = element_blank(), axis.title.y = element_blank(), axis.line.x = element_line(color = "black"), 
                    legend.text = element_text(size = 10)) + coord_flip()
                  print(p1 + p2 + p3)
                  ggsave_fun(filename = paste0(dir, "/5.MergeRatio"), plot = p1 + p2 + p3, width = 15, height = 6)
                  pB2_df <- table(scRNA@meta.data$cellType_1, scRNA@meta.data$group) %>% melt()
                  colnames(pB2_df) <- c("Cluster", "Sample", "Number")
                  pB2_df$Cluster <- factor(pB2_df$Cluster)
                  pB2_df <- pB2_df %>% group_by(Sample) %>% mutate(Percentage = Number/sum(Number))
                  pB2_df_tumor <- pB2_df %>% filter(Sample %in% c("Disease", "Tumor", "Cancer")) %>% select(Cluster, Tumor_Number = Number)
                  pB2_df_control <- pB2_df %>% filter(Sample == c("Control", "Normal")) %>% select(Cluster, Control_Number = Number)
                  pB2_df_combined <- merge(pB2_df_tumor, pB2_df_control, by = "Cluster")
                  pB2_df_combined$rate <- (pB2_df_combined$Tumor_Number - pB2_df_combined$Control_Number)/pB2_df_combined$Control_Number
                  KeyCell <- pB2_df_combined[abs(pB2_df_combined$rate) == max(abs(pB2_df_combined$rate)), ]$Cluster %>% as.character()
                  write.csv(pB2_df_combined, paste0(dir, "/Cell.rate.csv"))
                }, error = function(e) {
                  cellRatio <- cellRatioPlot(object = scRNA, sample.name = "orig.ident", celltype.name = "cellType_1", flow.curve = 0.5, 
                    fill.col = my_pal2) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
                  print(cellRatio)
                  ggsave_fun(filename = paste0(dir, "/3.Sample.ratio.plot"), plot = cellRatio, width = 7 + 0.25 * (length(unique(scRNA$orig.ident))), 
                    height = 7)
                })
            }
        }
        if (step == 4) {
            require(CellChat)
            dir <- paste0(Path, "/.Cellchat")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            options(stringsAsFactors = FALSE)
            cellchat <- createCellChat(object = scRNA, group.by = "cellType_1")
            cellchat
            groupSize <- as.numeric(table(cellchat@idents))
            if (mouse) {
                CellChatDB <- CellChatDB.mouse
            }
            else {
                CellChatDB <- CellChatDB.human
            }
            CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
            cellchat@DB <- CellChatDB.use
            cellchat <- subsetData(cellchat)
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            cellchat <- projectData(cellchat, PPI.human)
            cellchat <- computeCommunProb(cellchat, seed.use = 123456)
            cellchat <- filterCommunication(cellchat)
            df.net <- subsetCommunication(cellchat)
            write.table(df.net, file = paste0(dir, "/net_lr.txt"), quote = F, sep = "\t", row.names = F)
            cellchat <- computeCommunProbPathway(cellchat)
            df.netp <- subsetCommunication(cellchat)
            write.table(df.netp, file = paste0(dir, "/net_pathway.txt"), quote = F, sep = "\t", row.names = F)
            cellchat <- aggregateNet(cellchat)
            saveRDS(cellchat, paste0(dir, "/cellchat.rds"))
            groupSize <- as.numeric(table(cellchat@idents))
            pdf(file = paste0(dir, "/1.Net_number_strength.pdf"), width = 10, height = 5)
            par(mfrow = c(1, 2), xpd = TRUE)
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
            netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
            dev.off()
            df.net <- read.table(paste0(dir, "/net_lr.txt"), sep = "\t", check.names = F, header = T)
            data <- as.data.frame(table(c(df.net$source, df.net$target)))
            colnames(data) <- c("Cell_Type", "all_sum")
            data <- data[order(data$all_sum, decreasing = T), ]
            data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)
            head(data)
            sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
            sample_color2 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], 
                ifelse(data$all_sum > summary(data$all_sum)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
            sample_color1 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], 
                ifelse(data$all_sum > summary(data$all_sum)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), 
                alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))
            pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) + geom_bar(stat = "identity", 
                width = 0.8) + scale_fill_manual(values = sample_color1) + scale_colour_manual(values = sample_color2) + 
                theme_bw() + theme(panel.grid = element_blank()) + labs(x = "Cell Type", y = "Counts") + theme(axis.text.y = element_text(size = 12, 
                colour = "black")) + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
            ggsave_fun(paste0(dir, "/2.Interaction Count"), plot = pB2, width = 15, height = 10)
        }
        if (step == 5) {
            dir <- paste0(Path, "/.Subpopulation contribution")
            if (!dir.exists(paste0(dir))) {
                dir.create(paste0(dir))
            }
            scRNA$group <- factor(scRNA$group, c(grep("Control|Normal", unique(scRNA$group), value = T), grep("Disease|Tumor|Cancer", 
                unique(scRNA$group), value = T)))
            print(names(table(scRNA$group)))
            if (names(table(scRNA$group))[2] %in% c("Disease", "Tumor", "Cancer") & names(table(scRNA$group))[1] %in% c("Control", 
                "Normal")) {
                Bulk.DEGs <- FindMarkers(scRNA, group.by = "group", min.pct = 0.25, logfc.threshold = logFCfilter, ident.1 = names(table(scRNA$group))[2], 
                  ident.2 = names(table(scRNA$group))[1])
                Bulk.DEGs <- Bulk.DEGs[Bulk.DEGs$p_val_adj < adjPvalFilter, ]
                Bulk.DEGs <- Bulk.DEGs[Bulk.DEGs$avg_log2FC > 0, ]
                print(nrow(Bulk.DEGs))
                Bulk.DEGs$symbol <- rownames(Bulk.DEGs)
                subset.DEGs <- lapply(SplitObject(scRNA, split.by = "cellType_1"), function(subset) {
                  deg <- FindMarkers(subset, features = Bulk.DEGs$symbol, min.pct = 0, logfc.threshold = 0, group.by = "group", 
                    ident.1 = names(table(scRNA$group))[2], ident.2 = names(table(scRNA$group))[1])
                  deg$symbol <- rownames(deg)
                  deg$celltype <- unique(subset$cellType_1)
                  return(deg)
                })
                subset.DEGs <- do.call(rbind, subset.DEGs)
                subset.DEGs$FCexp <- 2^subset.DEGs$avg_log2FC
                subset.DEGs$FCprop <- subset.DEGs$pct.1/subset.DEGs$pct.2
                subset.DEGs$FCscore <- sqrt(subset.DEGs$FCexp * subset.DEGs$FCprop)
                subset.DEGs$FCscore[is.infinite(subset.DEGs$FCscore)] <- NA
                FCscore <- dcast(subset.DEGs, symbol ~ celltype, measure.var = "FCscore")
                write.table(FCscore, paste0(dir, "/output_FCscore.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
                plot.data <- data.frame(subtypes = colnames(FCscore[, -1]), FCscore = colMeans(FCscore[, -1], na.rm = T))
                plot.data <- arrange(plot.data, plot.data$FCscore)
                plot.data$subtypes <- factor(plot.data$subtypes, levels = plot.data$subtypes)
                my_pal2 <- c("#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", 
                  "#9ebcda", "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00")
                p <- ggplot(plot.data, aes(x = subtypes, y = FCscore, fill = subtypes)) + geom_bar(stat = "identity") + scale_fill_manual(values = my_pal2[1:nrow(plot.data)]) + 
                  geom_hline(yintercept = median(plot.data$FCscore), color = "grey", linetype = "dashed") + coord_polar() + 
                  theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), 
                  axis.ticks = element_blank())
                plot.data <- plot.data[order(plot.data$FCscore, decreasing = T), ]
                plot.data$subtypes <- factor(plot.data$subtypes, levels = plot.data$subtypes)
                sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
                sample_color2 <- ifelse(plot.data$FCscore > summary(plot.data$FCscore)[2], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[3], 
                  ifelse(plot.data$FCscore > summary(plot.data$FCscore)[5], sample_color[1], sample_color[2]), sample_color[3]), 
                  sample_color[4])
                sample_color1 <- ifelse(plot.data$FCscore > summary(plot.data$FCscore)[2], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[3], 
                  ifelse(plot.data$FCscore > summary(plot.data$FCscore)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 
                    0.9)), alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))
                pB2 <- ggplot(data = plot.data, aes(x = subtypes, y = FCscore, fill = subtypes, colour = subtypes)) + geom_bar(stat = "identity", 
                  width = 0.8) + scale_fill_manual(values = sample_color1) + scale_colour_manual(values = sample_color2) + 
                  theme_bw() + theme(panel.grid = element_blank()) + labs(x = "Cell Type", y = "Counts") + theme(axis.text.y = element_text(size = 12, 
                  colour = "black")) + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
                ggsave_fun(filename = paste0(dir, "/1.CombinedContribution"), plot = p + pB2, width = 12, height = 6)
                write.csv(plot.data, paste0(dir, "/plot.data.csv"), row.names = F)
            }
        }
        if (step == 6) {
            require(monocle)
            dir <- paste0(Path, "/.Pseudo-time")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            pbmcSub <- scRNA
            Idents(pbmcSub) <- pbmcSub@meta.data$cellType_1
            if (!is.na(KeyCell)) {
                pbmcSub <- subset(pbmcSub, cellType_1 %in% KeyCell)
                assay <- "RNA"
                set.resolutions <- seq(0.2, 1.2, by = 0.1)
                PC <- 25
                nfeatures <- 2000
                npcs <- 50
                pbmcSub <- NormalizeData(object = pbmcSub)
                pbmcSub <- FindVariableFeatures(object = pbmcSub, selection.method = "vst", nfeatures = nfeatures)
                pbmcSub <- ScaleData(object = pbmcSub, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
                pbmcSub <- RunPCA(pbmcSub, verbose = T, npcs = npcs, features = VariableFeatures(object = pbmcSub))
                if (Multi) {
                  pbmcSub <- RunHarmony(object = pbmcSub, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
                  pbmcSub <- RunUMAP(pbmcSub, reduction = "harmony", dims = 1:PC, verbose = T)
                  pbmcSub <- FindNeighbors(pbmcSub, dims = 1:PC, reduction = "harmony", verbose = T)
                  pbmcSub <- FindClusters(pbmcSub, resolution = set.resolutions, verbose = T)
                }
                else {
                  pbmcSub <- RunUMAP(pbmcSub, reduction = "pca", dims = 1:PC, verbose = T)
                  pbmcSub <- FindNeighbors(pbmcSub, dims = 1:PC, reduction = "pca", verbose = T)
                  pbmcSub <- FindClusters(pbmcSub, resolution = set.resolutions, verbose = T)
                }
                saveRDS(pbmcSub, paste0(dir, "/pbmcSub.rds"))
                pbmcSub[["celltype"]] <- as.data.frame(Idents(pbmcSub))
                pbmc.markers <- FindAllMarkers(object = pbmcSub, only.pos = FALSE, min.pct = 0.25, logfc.threshold = logFCfilter)
                sig.markers <- pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) > logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj)) < 
                  adjPvalFilter), ]
                write.table(cbind(Symbol = rownames(sig.markers), sig.markers), file = paste0(dir, "/ClusterMarkers.txt"), 
                  sep = "\t", row.names = F, quote = F)
            }
            if (is.na(KeyCell)) {
                pbmc.markers <- read.csv(paste0(CellAnnotate.dir, "/Findall.markers.cellType_1.csv"), header = T, row.names = 1)
                sig.markers <- pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) > logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj)) < 
                  adjPvalFilter), ]
            }
            else {
                pbmcSub <- subset(scRNA, subset = cellType_1 == KeyCell)
            }
            monocle.matrix <- GetAssayData(object = pbmcSub, slot = "data", assay = "RNA")
            monocle.sample <- pbmcSub@meta.data
            monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
            data <- as(as.matrix(monocle.matrix), "sparseMatrix")
            pd <- new("AnnotatedDataFrame", data = monocle.sample)
            fd <- new("AnnotatedDataFrame", data = monocle.geneAnn)
            cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
            if ("group" %in% names(pData(cds))) {
                names(pData(cds))[names(pData(cds)) == "group"] <- "Cluster"
            }
            else {
                names(pData(cds))[names(pData(cds)) == "orig.ident"] <- "Cluster"
            }
            pData(cds)[, "Cluster"] <- paste0("cluster", pData(cds)[, "Cluster"])
            cds <- estimateSizeFactors(cds)
            cds <- estimateDispersions(cds)
            disp_table <- dispersionTable(cds)
            disp_table <- arrange(disp_table, -dispersion_empirical)
            track_gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
            cds <- setOrderingFilter(cds, track_gene)
            cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
            cds <- orderCells(cds)
            p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
            ggsave(paste0(dir, "/1.Trajectory.Pseudotime.pdf"), plot = p)
            p <- plot_cell_trajectory(cds, color_by = "State")
            ggsave(paste0(dir, "/2.Trajectory.State.pdf"), plot = p)
            p <- plot_cell_trajectory(cds, color_by = "Cluster")
            ggsave(paste0(dir, "/3.Trajectory.Cluster.pdf"), plot = p)
            save(cds, file = paste0(dir, "/cds.rda"))
            if (!is.null(dev.list())) {
                dev.off()
            }
            if (!class(KeyGene) == "logical") {
                cds_subset <- cds[KeyGene, ]
                p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Cluster")
                p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
                p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
                p1 | p2 | p3
                ggsave(filename = paste0(dir, "/5.KeyGene.pseudotime.pdf"), plot = p1 | p2 | p3, width = 12, height = length(KeyGene) * 
                  2)
            }
        }
        if (step == 7) {
            dir <- paste0(Path, "/.KeyGene expression abundance")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            KeyGene <- intersect(rownames(scRNA), KeyGene)
            leg <- length(unique(KeyGene))
            p <- FeatureCornerAxes(object = scRNA, reduction = "umap", groupFacet = NULL, relLength = 0.5, relDist = 0.2, 
                features = unique(KeyGene), nLayout = 2, aspect.ratio = 1, pSize = 0.1)
            if (leg <= 3) {
                ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = leg * 5, height = 5)
            }
            if (leg > 3) {
                ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = 15, height = ceiling(leg/3) * 
                  5)
            }
            p <- DotPlot(scRNA, group.by = "cellType_1", features = unique(KeyGene), cols = "RdYlBu", ) + scale_size_continuous(range = c(0, 
                10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, 
                vjust = 0.5), legend.position = "right", text = element_text(size = 10))
            if (leg <= 3) {
                ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p, width = 5 + ceiling(leg/3), height = 7)
            }
            if (leg > 3) {
                ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p, width = 5 + 0.5 * leg, height = 7)
            }
        }
        if (step == 8) {
            dir <- paste0(Path, "/.Disease gene co-expression")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            showGenes <- intersect(rownames(scRNA), KeyGene)
            geneCard <- scan(sep = "\n", what = "character", quiet = FALSE)
            geneCard <- intersect(geneCard, rownames(scRNA))[1:3]
            for (i in showGenes) {
                PlotNum <- 0
                path <- paste0(dir, "/", which(showGenes == i), ".", i)
                if (!dir.exists(path)) {
                  dir.create(path)
                }
                for (j in geneCard) {
                  PlotNum <- PlotNum + 1
                  p1 <- FeaturePlot(scRNA, features = c(i, j), blend = TRUE, cols = c("gray80", "red", "green"), pt.size = 0.5, 
                    raster = F, reduction = "umap") + theme(aspect.ratio = 1)
                  p2 <- getScatterplot(scRNA, gene1 = j, gene2 = i, jitter.num = 0.15, pos = TRUE) + theme(aspect.ratio = 1)
                  p <- CombinePlots(plots = list(p1, p2), ncol = 2, rel_widths = c(4, 1))
                  ggsave_fun(paste0(path, "/", PlotNum, ".", j, " ~ ", i), plot = p, width = 15, height = 4)
                }
            }
        }
        if (step == 9) {
            dir <- paste0(Path, "/.Immunometabolic pathways")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            ShowGen <- intersect(rownames(scRNA), KeyGene)
            Idents(scRNA) <- "cellType_1"
            geneset <- clusterProfiler::read.gmt("refdata/Homo/h.all.v7.5.1.symbols.gmt")
            if (mouse) {
                geneset <- clusterProfiler::read.gmt("refdata/Mus/mh.all.v2023.2.Mm.symbols.gmt")
            }
            geneset <- split(geneset$gene, geneset$term)
            genesetInfo <- read.delim("refdata/GenesetInfo.txt", sep = ",")
            genesetInfo <- subset(genesetInfo, Classification != "")
            levels <- c("Immune", "Metabolism", "Signaling", "Proliferation")
            genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
            genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
            genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
            mat <- GetAssayData(object = scRNA, assay = "RNA", slot = "data")
            cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 1, plotStats = F)
            score <- AUCell::AUCell_calcAUC(geneset, cells_rankings, nCores = 1, aucMaxRank = nrow(cells_rankings) * 0.05)
            score <- AUCell::getAUC(score)
            DP.list <- base::lapply(ShowGen, function(cellType) {
                DatGroup <- FetchData(scRNA, vars = c(cellType), slot = "data")
                group <- setNames(object = ifelse(DatGroup[, 1] > median(DatGroup[, 1]), "Hexp", "Lexp"), nm = rownames(DatGroup))
                if (length(unique(group)) > 1) {
                  design <- model.matrix(~0 + factor(group))
                  colnames(design) <- levels(factor(group))
                  rownames(design) <- names(group)
                  contrast.matrix <- limma::makeContrasts("Hexp-Lexp", levels = design)
                  fit <- limma::lmFit(score[, names(group)], design)
                  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
                  fit2 <- limma::eBayes(fit2)
                  DPs <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "bonferroni")
                  DPs$KeyGenes <- cellType
                  DPs$Pathway <- rownames(DPs)
                  return(DPs)
                }
            })
            DPs <- do.call(rbind, DP.list)
            write.table(DPs, file = paste0(dir, "/output_hallmark.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
            plot.data <- DPs
            plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
            plot.data$KeyGenes <- factor(plot.data$KeyGenes)
            plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset)
            plot.data$FDR <- cut(plot.data$adj.P.Val, breaks = c(0, 1e-125, 1e-75, 1e-25, 1), include.lowest = T)
            plot.data$FDR <- factor(as.character(plot.data$FDR), levels = rev(levels(plot.data$FDR)))
            levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))
            color <- c("#4682B4", "#FFFFFF", "#CD2626")
            class.color <- c(Immune = "#D58986", Metabolism = "#80554C", Signaling = "#71AC7A", Proliferation = "#E8D4B4")
            p1 <- ggplot(plot.data, aes(x = Pathway, y = KeyGenes, color = logFC, size = FDR)) + geom_point() + scale_color_gradient2(low = color[1], 
                mid = color[2], high = color[3]) + geom_hline(yintercept = seq(min(as.numeric(plot.data$KeyGenes)) - 0.5, 
                max(as.numeric(plot.data$KeyGenes)) + 0.5), color = "grey80") + geom_vline(xintercept = seq(min(as.numeric(plot.data$Pathway)) - 
                0.5, max(as.numeric(plot.data$Pathway)) + 0.5), color = "grey80") + theme_classic() + theme(axis.text.x = element_text(angle = 90, 
                hjust = 1, vjust = 0.5), axis.line = element_blank(), legend.position = "top", legend.key.height = unit(0.3, 
                "cm"), legend.key.width = unit(0.8, "cm"), )
            p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) + geom_tile() + theme_classic() + scale_fill_manual(values = class.color) + 
                theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
                  legend.position = "bottom")
            p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(10, 2))
            ggsave_fun(filename = paste0(dir, "/1.Hallmark"), width = 15, height = 9, plot = p)
        }
        if (step == 10) {
            require(GSVA)
            require(GSEABase)
            require(limma)
            dir <- paste0(Path, "/.Metabolic pathway")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            gmtFile <- "refdata/Homo/immune.gmt"
            Metagroup <- scan(what = "character")
            rt <- scRNA@assays$RNA@data
            group <- data.frame(ID = rownames(scRNA@meta.data), Group = scRNA@meta.data[[Metagroup]])
            exp <- as.matrix(rt)
            dimnames <- list(rownames(exp), colnames(exp))
            mat <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
            mat <- avereps(mat)
            mat <- mat[rowMeans(mat) > 0, ]
            geneSet <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
            ssgseaScore <- gsva(mat, geneSet, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
            normalize <- function(x) {
                return((x - min(x))/(max(x) - min(x)))
            }
            ssgseaOut <- normalize(ssgseaScore)
            ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut)
            write.table(ssgseaOut, file = paste0(dir, "/ssgseaOut.txt"), sep = "\t", quote = F, col.names = F)
            plotdata <- read.table(paste0(dir, "/ssgseaOut.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
            up <- group
            up <- up[order(up$Group), ]
            plotdata <- plotdata[, up$ID]
            up$group <- up$Group
            down <- read.table("refdata/Homo/type.txt", sep = "\t", header = T, check.names = F, stringsAsFactors = F)
            compath <- intersect(rownames(plotdata), down$pathway)
            down <- down[down$pathway %in% compath, ]
            plotdata <- plotdata[down$pathway, ]
            annCol <- data.frame(Group = up$group, row.names = up$ID, stringsAsFactors = F)
            annRow <- data.frame(Direct = down$type, row.names = down$pathway, stringsAsFactors = F)
            annGroup <- c("#E7B800", "#2E9FDF")
            names(annGroup) <- unique(group$Group)
            annColors <- list(annGroup, Direct = c(`1Amino acid metabolism relevant signatures` = "#247BA0", `2lipid metabolism relevant signatures` = "#70C1B3", 
                `3Drug metabolism relevant signatures` = "#B2DBBF", `4Other metabolism signatures` = "#F3FFBD", `5C3 specific metabolism signatures` = "#6CD3A7"))
            plotdata <- t(scale(t(plotdata)))
            plotdata[plotdata > 1] <- 1
            plotdata[plotdata < -1] <- -1
            pdf(file = paste0(dir, "/1.Heatmap.pdf"), width = 25, height = 16)
            pheatmap(plotdata, scale = "none", annotation_row = annRow, annotation_col = annCol, annotation_colors = annColors, 
                color = colorRampPalette(c("#009BC7", "#F3F3F1", "#F15E4C"))(20), fontsize_row = 12, fontsize_col = 8, fontsize = 12, 
                cluster_cols = FALSE, cluster_rows = FALSE, gaps_row = c(13, 22, 25, 38), show_colnames = F)
            dev.off()
        }
        if (step == 11) {
            require(GSVA)
            require(GSEABase)
            require(limma)
            dir <- paste0(Path, "/.Hotspot mechanism")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            normalize <- function(x) {
                return((x - min(x))/(max(x) - min(x)))
            }
            mat <- GetAssayData(object = scRNA, assay = "RNA")
            mat <- as(mat[rowMeans(mat) > 0.3, ], "matrix")
            geneSet <- GSEABase::getGmt("ssGSEA.gmt", geneIdType = GSEABase::SymbolIdentifier())
            ssgseaScore <- gsva(mat, geneSet, method = "ssgsea", parallel.sz = 10, kcdf = "Gaussian", abs.ranking = TRUE)
            ssgseaOut <- normalize(ssgseaScore)
            save(ssgseaOut, file = paste0(dir, "/ssgseaOut.rda"))
            ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut)
            write.table(ssgseaOut, file = paste0(dir, "/ssgseaOutAll.txt"), sep = "\t", quote = F, col.names = F)
            load(paste0(dir, "/ssgseaOut.rda"))
            ssgseaOut <- as.data.frame(t(ssgseaOut))
            ScoGroup <- ifelse(ssgseaOut[rownames(scRNA[["orig.ident"]]), names(geneSet)] > median(ssgseaOut[rownames(scRNA[["orig.ident"]]), 
                names(geneSet)]), "Hsco", "Lsco")
            scRNA$ScoGroup <- data.frame(row.names = rownames(scRNA[["orig.ident"]]), ScoGroup)
            load(paste0(dir, "/ssgseaOut.rda"))
            meta <- as.data.frame(t(ssgseaOut))
            scRNA <- AddMetaData(scRNA, meta)
            saveRDS(scRNA, paste0(dir, "/result.rds"))
            tryCatch({
                for (i in rownames(ssgseaOut)) {
                  tmp <- scRNA@meta.data[, c("cellType_1", "group", i)]
                  colnames(tmp)[3] <- "val"
                  p1 <- ggplot(tmp, aes(x = cellType_1, y = val, fill = group, color = group)) + geom_boxplot(notch = F, 
                    alpha = 0.95, outlier.shape = 16, outlier.size = 0.65) + xlab("") + ylab("") + scale_fill_manual(values = c("#D5EBFB", 
                    "#FBEEB7", "#B4FBCD", "#F5B3FC")) + scale_color_manual(values = c("#0073C2", "#EFC000", "#00C244", "#C501D7")) + 
                    ggtitle("") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, 
                    size = 12), axis.title.y = element_text(angle = 90, size = 15)) + theme(legend.position = "top") + stat_compare_means(method = "wilcox.test", 
                    hjust = 0.5, vjust = 0, hide.ns = T, label = "p.signif", show.legend = F)
                  p2 <- FeaturePlot(object = scRNA, ncol = 1, features = i, cols = c("green", "red"))
                  p <- CombinePlots(plots = list(p2, p1), rel_widths = c(1, 1.6))
                  ggsave_fun(paste0(dir, "/1.", i), width = 12, height = 4, plot = p)
                }
            }, error = function(e) {
                for (i in rownames(ssgseaOut)) {
                  tmp <- scRNA@meta.data[, c("cellType_1", i)]
                  colnames(tmp)[2] <- "val"
                  p1 <- ggplot(tmp, aes(x = cellType_1, y = val, fill = cellType_1)) + geom_boxplot(notch = F, alpha = 0.95, 
                    outlier.shape = 16, outlier.size = 0.65) + xlab("") + ylab("") + ggtitle("") + theme_classic() + theme(axis.text.x = element_text(angle = 45, 
                    hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, size = 12), axis.title.y = element_text(angle = 90, 
                    size = 15)) + theme(legend.position = "top")
                  p2 <- FeaturePlot(object = scRNA, ncol = 1, features = i, cols = c("green", "red"))
                  p <- CombinePlots(plots = list(p2, p1), rel_widths = c(1, 1.6))
                  ggsave_fun(paste0(dir, "/1.", i), width = 12, height = 4, plot = p)
                }
            })
            isSlect <- scan(what = "character", quiet = T)
            if (isSlect == "yes") {
                sig.Cell <- scan(what = "character", quiet = T)
                pbmc <- subset(scRNA, subset = cellType_1 %in% sig.Cell)
            }
            if (isSlect == "no") {
                pbmc <- scRNA
            }
            PbmcMarkers <- FindMarkers(pbmc, ident.1 = "Hsco", min.pct = 0.1, group.by = "ScoGroup", assay = "RNA", logfc.threshold = 0.1)
            write.csv(PbmcMarkers, file = paste0(dir, "/DiffGene.csv"), row.names = T, quote = F)
            PbmcMarkers <- read.csv(paste0(dir, "/DiffGene.csv"), header = T, row.names = 1, check.names = F)
            require(ggrepel)
            cols <- c("#DC143C", "#00008B", "#808080")
            names(cols) <- c("Up", "Down", "NoSignifi")
            PbmcMarkers <- PbmcMarkers[-log10(PbmcMarkers$p_val_adj) > 0, ]
            PbmcMarkers$gene <- rownames(PbmcMarkers)
            PbmcMarkers$lab <- ""
            PbmcMarkers[order(PbmcMarkers$avg_log2FC), ][c(1:10, (nrow(PbmcMarkers) - 9):nrow(PbmcMarkers)), ]$lab <- PbmcMarkers[order(PbmcMarkers$avg_log2FC), 
                ][c(1:10, (nrow(PbmcMarkers) - 9):nrow(PbmcMarkers)), ]$gene
            PbmcMarkers$logP <- -log10(PbmcMarkers$p_val_adj)
            p <- ggplot(PbmcMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC)) + geom_point(aes(size = logP), 
                alpha = 0.9) + scale_color_gradient2(low = "#0500FF", high = "#FF0000", mid = "#EFEFEF", midpoint = 0, space = "Lab", 
                name = "logFC") + scale_size("-log10(p)") + geom_text_repel(data = PbmcMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                label = lab), size = 4, box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"), segment.color = "black", 
                show.legend = F) + theme_bw() + ylab("-log10 (p_val_adj)") + xlab("avg_log2FC") + geom_vline(xintercept = c(-logFCfilter, 
                logFCfilter), lty = 3, col = "black", lwd = 0.5) + geom_hline(yintercept = -log10(adjPvalFilter), lty = 3, 
                col = "black", lwd = 0.5)
            p
            ggsave_fun(filename = paste0(dir, "/2.VlnPlot"), plot = p, width = 10, height = 8)
        }
        if (step == 12) {
            dir <- paste0(Path, "/.Enrichment analysis")
            if (!dir.exists(dir)) {
                dir.create(dir)
            }
            require(enrichplot)
            require(clusterProfiler)
            genes <- KeyGene
            go <<- kk <<- NULL
            if (mouse) {
                require(org.Mm.eg.db)
                entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound = NA)
                entrezIDs <- as.character(entrezIDs)
                out <- cbind(genes, entrezID = entrezIDs)
                write.table(out, file = paste0(dir, "/entrezIDs.txt"), sep = "\t", quote = F, row.names = F)
                go <<- enrichGO(gene = entrezIDs, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", 
                  readable = T)
                write.csv(go, file = paste0(dir, "/GO.csv"), quote = F, row.names = F)
                kk <<- enrichKEGG(gene = entrezIDs, organism = "mmu", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
                write.csv(kk, file = paste0(dir, "/KEGG.csv"), quote = F, row.names = F)
            }
            else {
                require(org.Hs.eg.db)
                entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
                entrezIDs <- as.character(entrezIDs)
                out <- cbind(genes, entrezID = entrezIDs)
                write.table(out, file = paste0(dir, "/entrezIDs.txt"), sep = "\t", quote = F, row.names = F)
                go <<- enrichGO(gene = entrezIDs, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", 
                  readable = T)
                write.csv(go, file = paste0(dir, "/GO.csv"), quote = F, row.names = F)
                kk <<- enrichKEGG(gene = entrezIDs, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
                write.csv(kk, file = paste0(dir, "/KEGG.csv"), quote = F, row.names = F)
            }
            if (nrow(go@result) > 0) {
                go.bar <- barplot(go, drop = TRUE, showCategory = 10, split = "ONTOLOGY", label_format = 100) + facet_grid(ONTOLOGY ~ 
                  ., scale = "free")
                ggsave_fun(filename = paste0(dir, "/1.GOBarplot"), plot = go.bar, width = 12, height = 15)
                go.dot <- dotplot(go, showCategory = 10, split = "ONTOLOGY", label_format = 100) + facet_grid(ONTOLOGY ~ 
                  ., scale = "free")
                ggsave_fun(filename = paste0(dir, "/2.GODotplot"), plot = go.dot, width = 12, height = 15)
            }
            if (nrow(kk@result) > 0) {
                if (mouse) {
                  kk@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kk@result$Description)
                }
                kk.bar <- barplot(kk, drop = TRUE, showCategory = 30, label_format = 100)
                ggsave_fun(filename = paste0(dir, "/3.KEGGBarplot"), plot = kk.bar, width = 12, height = 15)
                kk.dot <- dotplot(kk, showCategory = 30, label_format = 100)
                ggsave_fun(filename = paste0(dir, "/4.KEGGBarplot"), plot = kk.dot, width = 12, height = 15)
            }
        }
    }
    return(scRNA)
}
require(ClusterGVis)
require(dplyr)
require(jjAnno)
require(Seurat)
plan("multicore", workers = 8)
options(future.globals.maxSize = 5e+06 * 1024^2)
set.seed(123456)
ClusterGVis_module <- function(object, OutPath, scRNA.markers = NULL, celltype = "cellType_1", species = NA, OnlyPlot = F) {
    if (!dir.exists(paste0(OutPath, "/ClusterGvis"))) 
        dir.create(paste0(OutPath, "/ClusterGvis"))
    setwd(OutPath)
    if (!OnlyPlot) {
        if (!is.na(species)) {
            mouse <- ifelse(species == "mouse", T, F)
        }
        else {
            mouse <- F
        }
        Idents(object) <- celltype
        if (is.null(scRNA.markers)) {
            scRNA.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
            write.csv(scRNA.markers, paste0("ClusterGvis/Findall.markers.", celltype, ".csv"))
        }
        else {
            scRNA.markers <- scRNA.markers
        }
        scRNA.top.markers <- scRNA.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 5, wt = avg_log2FC)
        write.csv(scRNA.top.markers, "ClusterGvis/TopMarkers.csv")
        st.data <- prepareDataFromscRNA(object = object, diffData = scRNA.top.markers, showAverage = T)
        if (mouse) {
            require(org.Mm.eg.db)
            enrich <- enrichCluster(object = st.data, OrgDb = org.Mm.eg.db, type = "BP", organism = "mmu", pvalueCutoff = 0.05, 
                topn = 5, seed = 123456)
        }
        else {
            require(org.Hs.eg.db)
            enrich <- enrichCluster(object = st.data, OrgDb = org.Hs.eg.db, type = "BP", organism = "hsa", pvalueCutoff = 0.05, 
                topn = 5, seed = 123456)
        }
        saveRDS(st.data, paste0("ClusterGvis/Heatmap.plotdata.rds"))
        write.csv(enrich, paste0("ClusterGvis/Enrich.plotdata.csv"))
    }
    st.data <- readRDS(paste0("ClusterGvis/Heatmap.plotdata.rds"))
    enrich <- read.csv(paste0("ClusterGvis/Enrich.plotdata.csv"), row.names = 1)
    scRNA.top.markers <- read.csv("ClusterGvis/TopMarkers.csv")
    pdf("ClusterGvis/1.Heatmap.pdf", height = 15, width = 17, onefile = F)
    visCluster(object = st.data, plot.type = "both", column_title_rot = 45, markGenes = unique(scRNA.top.markers$gene), markGenes.side = "left", 
        annoTerm.data = enrich, line.side = "left", cluster.order = c(1:length(levels(object@meta.data[[celltype]]))), add.bar = T)
    dev.off()
}
require(Seurat)
require(monocle)
require(dplyr)
require(ggplot2)
require(scales)
require(ClusterGVis)
set.seed(5201314)
my_pal2 <- c("#df4a86", "#746ea3", "#009ecb", "#00827b", "#3d4d7c", "#ad341d", "#0056a0", "#d77b1c", "#077140", "#d34132", 
    "#698da5", "#ffd597", "#f1977f", "#828bae", "#82bfbd", "#c1000a", "#FF8E00", "#00B3F1", "#354270", "#85b38f")
monocle2_module <- function(object, OutPath, group_by, track_gene = NULL, ShowGene = NA, OnlyPlot = F) {
    setwd(OutPath)
    if (!dir.exists("Monocle2")) 
        dir.create("Monocle2")
    if (!OnlyPlot) {
        scRNASub <- object
        monocle.matrix <- GetAssayData(object = scRNASub, slot = "data", assay = "RNA")
        monocle.sample <- scRNASub@meta.data
        monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
        data <- as(as.matrix(monocle.matrix), "sparseMatrix")
        pd <- new("AnnotatedDataFrame", data = monocle.sample)
        fd <- new("AnnotatedDataFrame", data = monocle.geneAnn)
        cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
        names(pData(cds))[names(pData(cds)) == group_by] <- "Cluster"
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        DefaultAssay(scRNASub) <- "RNA"
        if (is.null(track_gene)) {
            track_gene = c("VariableFeatures", "all")
            genenubmers <- "all"
        }
        else if (length(track_gene) == 1) {
            genenubmers <- "all"
            if (!is.na(as.numeric(as.character(track_gene[1])))) {
                genenubmers <- as.numeric(as.character(track_gene[1]))
                track_gene <- c("VariableFeatures", genenubmers)
            }
            else if (track_gene %in% c("dispersion", "cluster", "VariableFeatures")) {
                track_gene <- c(track_gene, 2000)
            }
        }
        else if (length(track_gene) == 2 && !is.na(as.numeric(as.character(track_gene[2])))) {
            genenubmers <- as.numeric(as.character(track_gene[2]))
            if (!track_gene[1] %in% c("dispersion", "cluster", "VariableFeatures")) {
                track_gene <- c("VariableFeatures", 2000)
            }
        }
        if (length(track_gene) >= 200) {
            genenubmers <- length(track_gene)
            track_gene <- track_gene
        }
        else if (track_gene[1] == "dispersion") {
            disp_table <- dispersionTable(cds)
            disp_table <- arrange(disp_table, -dispersion_empirical)
            track_gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
        }
        else if (track_gene[1] == "cluster") {
            Idents(scRNASub) <- "seurat_clusters"
            deg.cluster <- FindAllMarkers(scRNASub)
            deg.cluster <- arrange(deg.cluster, p_val)
            track_gene <- subset(deg.cluster, p_val_adj < 0.05)$gene
        }
        else if (track_gene[1] == "VariableFeatures") {
            track_gene <- VariableFeatures(scRNASub)
        }
        if (length(track_gene) <= genenubmers | genenubmers == "all") {
            track_gene <- track_gene
        }
        else if (length(track_gene) >= genenubmers) {
            track_gene <- track_gene[1:genenubmers]
        }
        cds <- setOrderingFilter(cds, track_gene)
        cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
        cds <- orderCells(cds)
        saveRDS(cds, "Monocle2/cds.rds")
    }
    cds <- readRDS("Monocle2/cds.rds")
    p <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.1) + ggsci::scale_color_gsea() + theme(legend.position = "right", 
        text = element_text(face = "bold"))
    ggsave(paste0("Monocle2/1.Trajectory.Pseudotime.pdf"), plot = p, height = 4, width = 4.5)
    p <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.1) + guides(color = guide_legend(override.aes = list(alpha = 1, 
        size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right", text = element_text(face = "bold"))
    ggsave(paste0("Monocle2/2.Trajectory.State.pdf"), plot = p, height = 4, width = 4.5)
    p <- plot_cell_trajectory(cds, color_by = "Cluster", cell_size = 0.1) + guides(color = guide_legend(override.aes = list(alpha = 1, 
        size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right", text = element_text(face = "bold"))
    ggsave(paste0("Monocle2/3.Trajectory.Cluster.pdf"), plot = p, height = 4, width = 4.5)
    if (!class(ShowGene) == "logical") {
        cds_subset <- cds[ShowGene, ]
        p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Cluster") + scale_color_manual(values = my_pal2)
        p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "State") + scale_color_manual(values = my_pal2)
        p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") + ggsci::scale_color_gsea()
        p1 | p2 | p3
        ggsave(filename = paste0("Monocle2/4.ShowGene.pseudotime.pdf"), plot = p1 | p2 | p3, width = 12, height = length(ShowGene) * 
            2)
    }
    return(cds)
}
require(data.table)
require(igraph)
require(CellChat)
require(dplyr)
require(Seurat)
require(patchwork)
require(ggsci)
require(RColorBrewer)
require(ggpubr)
require(scales)
cellchat_func <- function(object, species = NA, celltype = "cellType_1", ShowCell = NA, OnlyPlot = F, AllPlot = F) {
    if (!OnlyPlot) {
        if (!is.na(species)) {
            mouse <- ifelse(species == "mouse", T, F)
        }
        else {
            mouse <- F
        }
        if (mouse) {
            CellChatDB <- CellChatDB.mouse
            PPI.use <- PPI.mouse
        }
        else {
            CellChatDB <- CellChatDB.human
            PPI.use <- PPI.human
        }
        data.input <- object@assays$RNA@data
        identity <- data.frame(labels = object[[celltype]], row.names = colnames(data.input))
        cellchat <- createCellChat(data.input, meta = identity, group.by = celltype)
        CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor"))
        cellchat@DB <- CellChatDB.use
        cellchat <- subsetData(cellchat)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        cellchat <- projectData(cellchat, PPI.use)
        cellchat <- computeCommunProb(cellchat)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        df.net <- subsetCommunication(cellchat)
        write.table(df.net, file = paste0("net_lr.txt"), quote = F, sep = "\t", row.names = F)
        cellchat <- computeCommunProbPathway(cellchat)
        df.netp <- subsetCommunication(cellchat)
        write.table(df.netp, file = paste0("net_pathway.txt"), quote = F, sep = "\t", row.names = F)
        cellchat <- aggregateNet(cellchat)
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
        pathways.shows <- cellchat@netP$pathways
        head(cellchat@LR$LRsig)
        saveRDS(cellchat, file = paste0("cellchat.rds"))
    }
    cellchat <- readRDS("cellchat.rds")
    groupSize <- as.numeric(table(cellchat@idents))
    pathways.shows <- cellchat@netP$pathways
    {
        network_overview_dir <- c("./network_overview")
        dir.create(network_overview_dir, recursive = T)
        pdf(file = paste0(network_overview_dir, "/Net_number_strength.pdf"), width = 10, height = 5)
        par(mfrow = c(1, 2))
        netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
        netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
        dev.off()
        file.copy(paste0(network_overview_dir, "/Net_number_strength.pdf"), "1.Net_number_strength.pdf")
        df.net <- read.table(paste0("net_lr.txt"), sep = "\t", check.names = F, header = T)
        data <- as.data.frame(table(c(df.net$source, df.net$target)))
        colnames(data) <- c("Cell_Type", "all_sum")
        data <- data[order(data$all_sum, decreasing = T), ]
        data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)
        head(data)
        sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
        sample_color2 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], 
            ifelse(data$all_sum > summary(data$all_sum)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
        sample_color1 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], 
            ifelse(data$all_sum > summary(data$all_sum)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), alpha(sample_color[3], 
                0.9)), alpha(sample_color[4], 0.9))
        pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) + geom_bar(stat = "identity", 
            width = 0.8) + scale_fill_manual(values = sample_color1) + scale_colour_manual(values = sample_color2) + theme_bw() + 
            theme(panel.grid = element_blank()) + labs(x = "Cell Type", y = "Counts") + theme(axis.text.y = element_text(size = 12, 
            colour = "black")) + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
        ggsave(paste0(network_overview_dir, "/Interaction Count.pdf"), plot = pB2, width = 10, height = 8)
        file.copy(paste0(network_overview_dir, "/Interaction Count.pdf"), "2.Interaction Count.pdf")
        if (!is.na(ShowCell)) {
            ident1 <- which(levels(cellchat@idents) %in% ShowCell)
            ident2 <- which(!levels(cellchat@idents) %in% ShowCell)
            gg1 <- netVisual_bubble(cellchat, sources.use = ident1, targets.use = ident2, angle.x = 45, remove.isolate = F, 
                font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling from ", ShowCell))
            gg2 <- netVisual_bubble(cellchat, sources.use = ident2, targets.use = ident1, angle.x = 45, remove.isolate = F, 
                font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling to ", ShowCell))
            ggsave(paste0(network_overview_dir, "/Bubble.pdf"), plot = gg1$gg.obj + gg2$gg.obj, width = 5 + length(levels(cellchat@idents)) * 
                1, height = 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05)
            file.copy(paste0(network_overview_dir, "/Bubble.pdf"), "3.Bubble.pdf")
        }
    }
    if (AllPlot) {
        {
            network_circle_dir <- c("./network_circle")
            dir.create(network_circle_dir, recursive = T)
            mat <- cellchat@net$weight
            for (i in 1:nrow(mat)) {
                mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
                mat2[i, ] <- mat[i, ]
                pdf(file = paste0(network_circle_dir, "/netVisual_weight_", rownames(mat)[i], ".pdf"), width = 10, height = 9)
                netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
                dev.off()
            }
        }
        {
            network_aggregate_dir <- c("./network_aggregate")
            dir.create(network_aggregate_dir, recursive = T)
            for (pathways.show in pathways.shows) {
                pdf(file = paste0(network_aggregate_dir, "/netVisual_aggregate_", pathways.show, ".pdf"), width = 10, height = 9)
                netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", pt.title = 50, vertex.label.cex = 1)
                dev.off()
            }
        }
        {
            network_heatmap_dir <- c("./network_heatmap")
            dir.create(network_heatmap_dir, recursive = T)
            for (pathways.show in pathways.shows) {
                pdf(file = paste0(network_heatmap_dir, "/netVisual_heatmap_", pathways.show, ".pdf"), width = 10, height = 9)
                print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", font.size = 14, font.size.title = 20))
                dev.off()
            }
        }
        {
            pathway_contribution_dir <- c("./pathway_contribution")
            dir.create(pathway_contribution_dir, recursive = T)
            for (pathways.show in pathways.shows) {
                p <- netAnalysis_contribution(cellchat, signaling = pathways.show, font.size = 15, font.size.title = 15, 
                  title = paste0("Contribution of each L-R pair in ", pathways.show)) + theme(text = element_text(face = "bold"))
                ggsave(filename = paste0(pathway_contribution_dir, "/netAnalysis_contribution_", pathways.show, ".pdf"), 
                  plot = p)
            }
        }
        {
            communication_bubble_dir <- c("./communication_bubble")
            dir.create(communication_bubble_dir, recursive = T)
            for (sources.use in levels(cellchat@idents)) {
                targets.use <- levels(cellchat@idents)
                p <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, angle.x = 45, remove.isolate = F, 
                  font.size = 10, font.size.title = 15, return.data = T, title.name = paste0("Signaling from ", sources.use))
                ggsave(filename = paste0(communication_bubble_dir, "/netVisual_bubble_from_", sources.use, ".pdf"), plot = p$gg.obj, 
                  width = 8, height = 8 + nrow(p$communication) * 0.05)
            }
            for (targets.use in levels(cellchat@idents)) {
                sources.use <- levels(cellchat@idents)
                p <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, angle.x = 45, remove.isolate = F, 
                  font.size = 10, font.size.title = 15, return.data = T, title.name = paste0("Signaling to ", targets.use))
                ggsave(filename = paste0(communication_bubble_dir, "/netVisual_bubble_to_", targets.use, ".pdf"), plot = p$gg.obj, 
                  width = 8, height = 8 + nrow(p$communication) * 0.05)
            }
        }
        {
            signaling_role_dir <- c("./signaling_role")
            dir.create(signaling_role_dir, recursive = T)
            gg1 <- netAnalysis_signalingRole_scatter(cellchat)
            ggsave(filename = paste0(signaling_role_dir, "/netAnalysis_signalingRole_scatter", ".pdf"), plot = gg1, width = 8, 
                height = 8)
            ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 7 + length(cellchat@netP$pathways) * 
                0.125)
            ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 7 + length(cellchat@netP$pathways) * 
                0.125)
            pdf(paste0(signaling_role_dir, "/netAnalysis_signalingRole_heatmap.pdf"), width = 15, height = 7 + length(cellchat@netP$pathways) * 
                0.125)
            print(ht1 + ht2)
            dev.off()
        }
        {
            gene_expression_dir <- c("./gene_expression")
            dir.create(gene_expression_dir, recursive = T)
            for (pathways.show in pathways.shows) {
                res <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE, enriched.only = T)
                p <- plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = TRUE, type = "violin")
                ggsave(filename = paste0(gene_expression_dir, "/gene_expression_", pathways.show, ".pdf"), plot = p, width = 10, 
                  height = 7 + length(res$geneLR) * 0.225)
            }
        }
    }
    return(cellchat)
}
cellchat_module <- function(object, OutPath, group_by = "all", celltype = "cellType_1", idents = NULL, ShowCell = NA, species = NA, 
    OnlyPlot = F, AllPlot = F) {
    if (group_by == "all" & is.null(idents)) {
        out_path <- paste0(OutPath, "/CellChat")
        if (!dir.exists(out_path)) 
            dir.create(out_path)
        setwd(out_path)
        cellchat <- cellchat_func(object = object, species = species, celltype = celltype, ShowCell = ShowCell, OnlyPlot = OnlyPlot, 
            AllPlot = AllPlot)
    }
    else if (group_by != "all" & length(idents) == 2) {
        for (ident in idents) {
            out_path <- paste0(OutPath, "/CellChat/", ident)
            if (!dir.exists(out_path)) 
                dir.create(out_path, recursive = T)
            setwd(out_path)
            sub.object <- object[, object@meta.data[[group_by]] %in% ident]
            cellchat <- cellchat_func(object = sub.object, species = species, celltype = celltype, ShowCell = ShowCell, OnlyPlot = OnlyPlot, 
                AllPlot = AllPlot)
        }
        setwd(paste0(OutPath, "/CellChat/"))
        cellchat.NL <- readRDS(paste0(idents[1], "/cellchat.rds"))
        cellchat.TL <- readRDS(paste0(idents[2], "/cellchat.rds"))
        object.list <- setNames(list(cellchat.NL, cellchat.TL), idents[1:2])
        cellchat <- mergeCellChat(object.list, add.names = names(object.list))
        save(cellchat, file = paste0("cellchat_merged_", idents[1], "_", idents[2], ".RData"))
        gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
        gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
        ggsave("1.CompareInteractions.pdf", plot = gg1 + gg2, width = 5, height = 5)
        pdf("2.NetVisual.pdf", width = 10, height = 5)
        par(mfrow = c(1, 2), xpd = TRUE)
        netVisual_diffInteraction(cellchat, weight.scale = T)
        netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
        dev.off()
        pdf("3.NetVisualHeatmap.pdf", width = 10, height = 5)
        gg1 <- netVisual_heatmap(cellchat)
        gg2 <- netVisual_heatmap(cellchat, measure = "weight")
        print(gg1 + gg2)
        dev.off()
        weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
        par(mfrow = c(1, 2), xpd = TRUE)
        for (i in 1:length(object.list)) {
            netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2], 
                edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
        }
        num.link <- sapply(object.list, function(x) {
            rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
        })
        weight.MinMax <- c(min(num.link), max(num.link))
        gg <- list()
        for (i in 1:length(object.list)) {
            gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
        }
        ggsave("4.SignalingRole.pdf", plot = patchwork::wrap_plots(plots = gg), width = 10, height = 5)
        if (!is.na(ShowCell)) {
            ident1 <- which(levels(object@meta.data[[celltype]]) %in% ShowCell)
            ident2 <- which(!levels(object@meta.data[[celltype]]) %in% ShowCell)
            gg1 <- netVisual_bubble(cellchat, sources.use = ident1, targets.use = ident2, comparison = c(1, 2), angle.x = 45, 
                remove.isolate = F, font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling from ", 
                  ShowCell))
            gg2 <- netVisual_bubble(cellchat, sources.use = ident2, targets.use = ident1, comparison = c(1, 2), angle.x = 45, 
                remove.isolate = F, font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling to ", 
                  ShowCell))
            ggsave("5.Bubble.pdf", plot = gg1$gg.obj + gg2$gg.obj, width = 5 + length(levels(object@meta.data[[celltype]])) * 
                1, height = 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05)
        }
    }
}
library(Seurat)
library(harmony)
library(ggplot2)
library(scRNAtoolVis)
library(dplyr)
human2mouse <- function(x) {
    x1 <- substr(x, 1, 1)
    x2 <- substr(x, 2, nchar(x))
    return(paste0(x1, tolower(x2)))
}
GeneExpre <- function(outpath = "./result/", SeuratObject = NA, KeyGene = NA) {
    d <- list.files(outpath, pattern = "KeyGene expression abundance", full.names = T)
    if (length(d) == 0) {
        Expr.path <- paste0(outpath, "/KeyGene expression abundance/")
        dir.create(Expr.path)
    }
    else {
        Expr.path <- d
    }
    scRNA <- SeuratObject
    KeyGene <- intersect(rownames(scRNA), KeyGene)
    leg <- length(unique(KeyGene))
    p <- featureCornerAxes(object = scRNA, reduction = "umap", groupFacet = NULL, relLength = 0.5, relDist = 0.2, features = unique(KeyGene), 
        nLayout = 3, aspect.ratio = 1, pSize = 0.1)
    if (leg <= 3) {
        ggsave(filename = paste0(Expr.path, "/1.Keygene.FeaturePlot.pdf"), plot = p, width = leg * 5, height = 5)
    }
    if (leg > 3) {
        ggsave(filename = paste0(Expr.path, "/1.Keygene.FeaturePlot.pdf"), plot = p, width = 15, height = ceiling(leg/3) * 
            5)
    }
    p <- DotPlot(scRNA, group.by = "cellType_1", features = unique(KeyGene), cols = "RdYlBu", ) + scale_size_continuous(range = c(0, 
        10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = "right", text = element_text(size = 10))
    if (leg <= 3) {
        ggsave(filename = paste0(Expr.path, "/2.Keygene.DotPlot.pdf"), plot = p, width = 5 + ceiling(leg/3), height = 7)
    }
    if (leg > 3) {
        ggsave(filename = paste0(Expr.path, "/2.Keygene.DotPlot.pdf"), plot = p, width = 5 + 0.5 * leg, height = 7)
    }
}
hallmark_module <- function(outpath = "./result/", SeuratObject = NA, KeyGene = NA) {
    e <- list.files(outpath, pattern = "Immunometabolic pathways", full.names = T)
    if (length(e) == 0) {
        hallmark.path <- paste0(outpath, "/Immunometabolic pathways/")
        dir.create(hallmark.path)
    }
    else {
        hallmark.path <- e
    }
    scRNA <- SeuratObject
    ShowGen <- intersect(rownames(scRNA), KeyGene)
    Idents(scRNA) <- "cellType_1"
    load("./data/hallmark.gs.rdata")
    genesetInfo <- read.delim("./data/GenesetInfo.txt", sep = ",")
    geneset <- gs
    rownames(genesetInfo) <- genesetInfo$geneset
    genesetInfo <- genesetInfo[names(geneset), ]
    genesetInfo <- subset(genesetInfo, Classification != "")
    levels <- c("Immune", "Metabolism", "Signaling", "Proliferation")
    genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
    genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
    genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
    mat <- GetAssayData(object = scRNA, assay = "RNA", slot = "data")
    cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 1, plotStats = F)
    score <- AUCell::AUCell_calcAUC(geneset, cells_rankings, nCores = 1, aucMaxRank = nrow(cells_rankings) * 0.05)
    score <- AUCell::getAUC(score)
    DP.list <- base::lapply(ShowGen, function(cellType) {
        DatGroup <- FetchData(scRNA, vars = c(cellType), slot = "data")
        group <- setNames(object = ifelse(DatGroup[, 1] > median(DatGroup[, 1]), "Hexp", "Lexp"), nm = rownames(DatGroup))
        if (length(unique(group)) > 1) {
            design <- model.matrix(~0 + factor(group))
            colnames(design) <- levels(factor(group))
            rownames(design) <- names(group)
            contrast.matrix <- limma::makeContrasts("Hexp-Lexp", levels = design)
            fit <- limma::lmFit(score[, names(group)], design)
            fit2 <- limma::contrasts.fit(fit, contrast.matrix)
            fit2 <- limma::eBayes(fit2)
            DPs <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "bonferroni")
            DPs$KeyGenes <- cellType
            DPs$Pathway <- rownames(DPs)
            return(DPs)
        }
    })
    DPs <- do.call(rbind, DP.list)
    plot.data <- DPs
    plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
    plot.data$KeyGenes <- factor(plot.data$KeyGenes)
    plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset)
    plot.data$FDR <- cut(plot.data$adj.P.Val, breaks = c(0, 1e-125, 1e-75, 1e-25, 1), include.lowest = T)
    plot.data$FDR <- factor(as.character(plot.data$FDR), levels = rev(levels(plot.data$FDR)))
    levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))
    color <- c("#4682B4", "#FFFFFF", "#CD2626")
    class.color <- c(Immune = "#D58986", Metabolism = "#80554C", Signaling = "#71AC7A", Proliferation = "#E8D4B4")
    p1 <- ggplot(plot.data, aes(x = Pathway, y = KeyGenes, color = logFC, size = FDR)) + geom_point() + scale_color_gradient2(low = color[1], 
        mid = color[2], high = color[3]) + geom_hline(yintercept = seq(min(as.numeric(plot.data$KeyGenes)) - 0.5, max(as.numeric(plot.data$KeyGenes)) + 
        0.5), color = "grey80") + geom_vline(xintercept = seq(min(as.numeric(plot.data$Pathway)) - 0.5, max(as.numeric(plot.data$Pathway)) + 
        0.5), color = "grey80") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_blank(), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
            "cm"), )
    p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) + geom_tile() + theme_classic() + scale_fill_manual(values = class.color) + 
        theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
            legend.position = "bottom")
    p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(10, 2))
    ggsave(filename = paste0(hallmark.path, "/1.Hallmark.pdf"), width = 15, height = 9, plot = p)
}
library(AUCell)
library(ggplot2)
library(reshape2)
library(dplyr)
BPS_plot <- function(BPS_AUC_path, BPS_AUC_P_path) {
    BPS_AUC <- read.table(BPS_AUC_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, 
        fill = TRUE, quote = "")
    BPS_AUC_P <- read.table(BPS_AUC_P_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, 
        fill = TRUE, quote = "")
    auc_cols <- setdiff(colnames(BPS_AUC), c("missing", "ncells"))
    BPS_AUC <- BPS_AUC[, auc_cols, drop = FALSE]
    auc_mat <- matrix(as.numeric(as.matrix(t(BPS_AUC))), nrow = 1)
    rownames(auc_mat) <- "BPS_AUC"
    cells_assignment <- AUCell_exploreThresholds(auc_mat, plotHist = FALSE, assign = TRUE)
    cutoff <- suppressWarnings(as.numeric(cells_assignment$BPS_AUC$aucThr$selected))
    pdf("1.BPS_AUC_histogram.pdf", width = 6, height = 4)
    hist_obj <- hist(as.vector(auc_mat), breaks = 207, col = "#adc7d9", border = NA, main = "", xlab = expression(BPS[AUC] ~ 
        " histogram"))
    abline(v = cutoff, col = "red", lwd = 3)
    text(cutoff, max(hist_obj$counts), labels = paste0("Cutoff = ", round(cutoff, 3)), pos = 4, col = "#000000")
    dev.off()
    BPS_FDR <- apply(BPS_AUC_P, 2, p.adjust, method = "fdr")
    BPS_AUC$celltype <- rownames(BPS_AUC)
    AUC_long <- reshape2::melt(BPS_AUC, id.vars = "celltype", variable.name = "trait", value.name = "AUC")
    FDR_long <- reshape2::melt(BPS_FDR, varnames = c("celltype", "trait"), value.name = "FDR")
    merged <- merge(AUC_long, FDR_long, by = c("celltype", "trait"))
    df <- merged
    df$neglogFDR <- -log10(df$FDR)
    df$class <- with(df, ifelse(FDR < 0.01 & AUC > 0.5, "Stringent", ifelse(FDR < 0.01 & AUC > cutoff & AUC < 0.5, "Moderate", 
        ifelse(FDR > 0.01 & FDR < 0.05 & AUC > cutoff & AUC < 0.5, "Lenient", "Nonsig"))))
    N_sig <- sum(df$class %in% c("Stringent", "Moderate", "Lenient"))
    n_stringent <- sum(df$class == "Stringent")
    n_moderate <- sum(df$class == "Moderate")
    n_lenient <- sum(df$class == "Lenient")
    cols <- c(Stringent = "#681c21", Moderate = "#427497", Lenient = "#619581", Nonsig = "grey80")
    table(df$class)
    p_scatter <- ggplot(df, aes(x = AUC, y = neglogFDR, color = class)) + geom_point(alpha = 0.8, size = 2) + scale_color_manual(values = cols) + 
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "#722c2c") + geom_vline(xintercept = cutoff, linetype = "dashed", 
        color = "#722c2c") + geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#722c2c") + geom_hline(yintercept = -log10(0.05), 
        linetype = "dashed", color = "#722c2c") + annotate("text", x = min(df$AUC, na.rm = TRUE) + 0.2, y = max(df$neglogFDR, 
        na.rm = TRUE) * 1.06, label = paste0("Moderate (N=", format(n_moderate, big.mark = ","), ")"), hjust = 0, vjust = 1, 
        color = cols["Moderate"], size = 4, fontface = "italic") + annotate("text", x = min(df$AUC, na.rm = TRUE) + 0.55, 
        y = max(df$neglogFDR, na.rm = TRUE) * 1.06, label = paste0("Stringent (N=", format(n_stringent, big.mark = ","), 
            ")"), hjust = 0, vjust = 1, color = cols["Stringent"], size = 4, fontface = "italic") + annotate("text", x = min(df$AUC, 
        na.rm = TRUE) + 0.45, y = max(df$neglogFDR, na.rm = TRUE) * 0.5, label = paste0("Lenient (N=", format(n_lenient, 
        big.mark = ","), ")"), hjust = 0, vjust = 1, color = cols["Lenient"], size = 4, fontface = "italic") + labs(x = expression(BPS[AUC] ~ 
        " scores"), y = expression(-log[10](FDR)), color = "Class") + theme_bw(base_size = 16) + theme(panel.grid = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(color = "black"))
    p_scatter <- p_scatter + theme(legend.position = "none")
    ggsave("2.BPS_AUC_scatter.pdf", p_scatter, width = 5, height = 5)
    bar_df <- df %>% filter(class %in% c("Stringent", "Moderate", "Lenient")) %>% group_by(celltype, class) %>% summarise(n = n(), 
        .groups = "drop") %>% ungroup()
    order_df <- bar_df %>% group_by(celltype) %>% summarise(total = sum(n), .groups = "drop") %>% arrange(desc(total))
    bar_df$celltype <- factor(bar_df$celltype, levels = order_df$celltype)
    bar_df$class <- factor(bar_df$class, levels = c("Stringent", "Moderate", "Lenient"))
    p_bar <- ggplot(bar_df, aes(x = celltype, y = n, fill = class)) + geom_bar(stat = "identity", width = 0.6) + scale_fill_manual(values = cols) + 
        labs(x = "", y = "Number", fill = "") + theme_bw(base_size = 16) + theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(color = "black"), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = c(0.55, 
            1.2), legend.justification = c(0, 1), legend.background = element_rect(fill = alpha("white", 0), color = "NA"), 
        legend.key = element_rect(fill = alpha("white", 0), color = NA)) + scale_y_continuous(expand = c(0, 0))
    width <- 2.5 + 0.25 * (length(unique(bar_df$celltype)))
    ggsave("3.BPS_AUC_bar.pdf", p_bar, width = width, height = 4)
    key_cell <- bar_df %>% group_by(celltype) %>% summarise(total = sum(n), .groups = "drop") %>% arrange(desc(total)) %>% 
        slice_head(n = 5) %>% pull(celltype) %>% head(1)
    write.table(key_cell, "key_cell.txt", row.names = FALSE)
    df_key <- df %>% filter(celltype == key_cell) %>% arrange(desc(AUC))
    write.table(df_key, "key_cell_BPS_AUC.txt", row.names = FALSE)
    pie_df <- df_key %>% filter(class %in% c("Stringent", "Moderate", "Lenient", "Nonsig")) %>% dplyr::count(class, name = "n")
    pie_df$class <- factor(pie_df$class, levels = c("Stringent", "Moderate", "Lenient", "Nonsig"))
    pie_df <- pie_df[order(pie_df$class), ]
    pie_df$label <- paste0("n = ", pie_df$n)
    total_n <- sum(pie_df$n)
    p_pie <- ggplot(pie_df, aes(x = 2, y = n, fill = class)) + geom_col(color = "black", width = 1) + coord_polar(theta = "y") + 
        xlim(0.5, 2.5) + scale_fill_manual(values = cols) + geom_text(aes(label = label), position = position_stack(vjust = 0.5), 
        color = "white", size = 5) + annotate("text", x = 0.55, y = 0, label = paste0("Total: ", total_n), color = "red", 
        size = 6, fontface = "bold") + labs(fill = "") + theme_void(base_size = 16) + theme(legend.position = "right")
    ggsave("4.BPS_AUC_pie.pdf", p_pie, width = 6, height = 4)
}
calculate_scbps_auc <- function(norm_score, annotation, auc_pvalue, auc_output) {
    library(scBPS)
    library(DelayedArray)
    library(data.table)
    library(reshape2)
    random_dir <- file.path(dirname(normalizePath(auc_pvalue, mustWork = FALSE)), "tmp")
    dir.create(random_dir, recursive = TRUE, showWarnings = FALSE)
    cell_annotation <- fread(annotation, header = TRUE, stringsAsFactors = FALSE)
    annotation_ids <- setdiff(as.character(unique(cell_annotation$cell_annotation)), "other")
    score <- data.frame(fread(norm_score, header = TRUE, stringsAsFactors = FALSE))
    rownames(score) <- score$cell_id
    score$cell_id <- NULL
    rank_score <- buildRankings(score)
    auc <- calc.AUC2(annotation_ids, cell_annotation, rank_score, 0.05)
    cell_counts <- reshape2::dcast(cell_annotation, cell_annotation ~ "num", fun.aggregate = length)
    rownames(cell_counts) <- cell_counts[, 1]
    lapply(annotation_ids, function(cell_type) {
        perm_cal2(cell_type, cell_counts, 1000, cell_annotation, rank_score, 0.05, mydir = random_dir)
    })
    p_values <- pvalue_AUC(auc, random_dir)
    fwrite(p_values, file = auc_pvalue, sep = "\t", row.names = TRUE, col.names = TRUE)
    fwrite(auc, file = auc_output, sep = "\t", row.names = TRUE, col.names = TRUE)
    invisible(TRUE)
}
run_scBPS_pipeline <- function(scrna_obj, celltype_col = "cellType_1", gwas = "Dutch207", output_dir = "./", Only_Plot = FALSE) {
    old_directory <- getwd()
    on.exit(setwd(old_directory), add = TRUE)
    setwd(output_dir)
    dir.create("scBPS", showWarnings = FALSE)
    setwd("scBPS")
    dir.create("output", showWarnings = FALSE)
    auc_output <- "output/BPS_AUC.txt"
    auc_pvalue <- "output/pvalue_AUC.txt"
    if (!Only_Plot) {
        if (!file.exists("output/sc_adata.h5ad")) {
            scRNA <- scrna_obj
            scRNA[["RNA"]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
            SaveH5Seurat(scRNA, filename = "output/sc_adata.h5Seurat", overwrite = TRUE)
            Convert("output/sc_adata.h5Seurat", dest = "h5ad", overwrite = TRUE)
            cell_annotation <- data.frame(cell_id = colnames(scRNA), cell_annotation = scRNA[[celltype_col]][, 1], stringsAsFactors = FALSE)
            write.table(cell_annotation, "output/cell.annotation.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        }
        temp_dir <- file.path("output", "tmp")
        dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
        zscore_file <- file.path(client_root, "data", "reference", "scBPS", paste0("gene_zscore_", gwas, ".txt"))
        conda_executable <- Sys.getenv("CONDA_EXE", "conda")
        scbps_python <- Sys.getenv("SCBPS_PYTHON", "")
        scbps_code <- file.path(client_root, "python", "scbps")
        compute_arguments <- c(file.path(scbps_code, "compute-score.py"), zscore_file, paste0(normalizePath(temp_dir, mustWork = TRUE), 
            "/"), normalizePath("output/sc_adata.h5ad", mustWork = TRUE))
        compute_status <- if (nzchar(scbps_python)) system2(scbps_python, compute_arguments) else system2(conda_executable, c("run", 
            "--no-capture-output", "-n", "scBPS", "python", compute_arguments))
        stopifnot(compute_status == 0L)
        traits <- fread(zscore_file, select = 1)[[1]]
        score_files <- normalizePath(file.path(temp_dir, paste0(traits, ".score")), mustWork = TRUE)
        writeLines(score_files, file.path(temp_dir, "list"))
        merge_arguments <- c(file.path(scbps_code, "generate_norm_score.py"), normalizePath(file.path(temp_dir, "list"), mustWork = TRUE), 
            normalizePath("output/norm_score.tsv", mustWork = FALSE))
        merge_status <- if (nzchar(scbps_python)) system2(scbps_python, merge_arguments) else system2(conda_executable, c("run", 
            "--no-capture-output", "-n", "scBPS", "python", merge_arguments))
        stopifnot(merge_status == 0L)
        calculate_scbps_auc("output/norm_score.tsv", "output/cell.annotation.txt", auc_pvalue, auc_output)
    }
    BPS_plot(BPS_AUC_path = auc_output, BPS_AUC_P_path = auc_pvalue)
    invisible(TRUE)
}
library(clusterProfiler)
library(ggrepel)
library(cowplot)
library(ggraph)
library(limma)
gather_graph_node <- function(df, index = NULL, value = tail(colnames(df), 1), root = NULL) {
    require(dplyr)
    if (!length(index) < 2) {
        list <- lapply(seq_along(index), function(i) {
            dots <- index[1:i]
            df %>% group_by(.dots = dots) %>% summarise(node.size = sum(.data[[value]]), node.level = index[[i]], node.count = n()) %>% 
                mutate(node.short_name = as.character(.data[[dots[[length(dots)]]]]), node.branch = as.character(.data[[dots[[1]]]])) %>% 
                tidyr::unite(node.name, dots, sep = "/")
        })
        data <- do.call("rbind", list) %>% as_tibble()
        data$node.level <- factor(data$node.level, levels = index)
        if (is.null(root)) {
            return(data)
        }
        else {
            root_data <- data.frame(node.name = root, node.size = sum(df[[value]]), node.level = root, node.count = 1, node.short_name = root, 
                node.branch = root, stringsAsFactors = F)
            data <- rbind(root_data, data)
            data$node.level <- factor(data$node.level, levels = c(root, index))
            return(data)
        }
    }
}
gather_graph_edge <- function(df, index = NULL, root = NULL) {
    require(dplyr)
    if (!length(index) < 2) 
        if (length(index) == 2) {
            data <- df %>% mutate(from = .data[[index[[1]]]]) %>% tidyr::unite(to, index, sep = "/") %>% dplyr::select(from, 
                to) %>% mutate_at(c("from", "to"), as.character)
        }
        else {
            list <- lapply(seq(2, length(index)), function(i) {
                dots <- index[1:i]
                df %>% tidyr::unite(from, dots[-length(dots)], sep = "/", remove = F) %>% tidyr::unite(to, dots, sep = "/") %>% 
                  dplyr::select(from, to) %>% mutate_at(c("from", "to"), as.character)
            })
            data <- do.call("rbind", list)
        }
    data <- as_tibble(data)
    if (is.null(root)) {
        return(data)
    }
    else {
        root_data <- df %>% group_by(.dots = index[[1]]) %>% summarise(count = n()) %>% mutate(from = root, to = as.character(.data[[index[[1]]]])) %>% 
            dplyr::select(from, to)
        rbind(root_data, data)
    }
}
GSEA_module <- function(mat, outpath = "./result/", species = NA, KeyGene = NA, isKEGG = TRUE) {
    dir <- list.files(outpath, pattern = "GSEA", full.names = T)
    if (length(dir) == 0) {
        dir <- paste0(outpath, "/GSEA/")
        dir.create(dir)
    }
    if (!is.na(species)) {
        mouse <- ifelse(species == "mouse", T, F)
    }
    else {
        mouse <- F
    }
    for (GeneIndex in 1:length(KeyGene)) {
        gene <- KeyGene[GeneIndex]
        low <- mat[gene, ] <= median(mat[gene, ])
        high <- mat[gene, ] > median(mat[gene, ])
        lowRT <- mat[, low]
        highRT <- mat[, high]
        conNum <- ncol(lowRT)
        treatNum <- ncol(highRT)
        mat <- cbind(lowRT, highRT)
        Type <- c(rep("con", conNum), rep("treat", treatNum))
        design <- model.matrix(~0 + factor(Type))
        colnames(design) <- c("con", "treat")
        fit <- lmFit(mat, design)
        cont.matrix <- makeContrasts(treat - con, levels = design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)
        allDiff <- topTable(fit2, adjust = "fdr", number = 2e+05)
        gsym.fc <- data.frame(SYMBOL = rownames(allDiff), logFC = allDiff$logFC)
        if (!dir.exists(paste0(dir, "/", GeneIndex, ".", gene))) {
            dir.create(paste0(dir, "/", GeneIndex, ".", gene))
        }
        selectedGeneID <- gene
        mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9", "#7A142C", "#5D90BA", 
            "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
        if (isKEGG) {
            if (mouse) {
                gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
                gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                id.fc <- gsym.fc.id.sorted$logFC
                names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                kk <- gseKEGG(id.fc, organism = "mmu")
                kk <- append_kegg_category(kk)
                kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
                sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
            }
            else {
                gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                id.fc <- gsym.fc.id.sorted$logFC
                names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                kk <- gseKEGG(id.fc, organism = "hsa")
                kk <- append_kegg_category(kk)
                kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
                sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
            }
            write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/gsea_output.csv"), quote = TRUE, row.names = F)
        }
        else {
            if (mouse) {
                gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
                gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                id.fc <- gsym.fc.id.sorted$logFC
                names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                kk <- gseGO(id.fc, OrgDb = "org.Mm.eg.db", pAdjustMethod = "none")
                kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
                sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
            }
            else {
                gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                id.fc <- gsym.fc.id.sorted$logFC
                names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                kk <- gseGO(id.fc, OrgDb = "org.Hs.eg.db", pAdjustMethod = "none")
                kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
                sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
            }
            write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/gsea_output.csv"), quote = TRUE, row.names = F)
            if (nrow(sortkk) < 3) {
                Sys.sleep(3)
                next
            }
        }
        if (isKEGG) {
            display_df <- as.data.frame(sortkk)
            immune_ids <- grepl("^KEGG_hsa046", display_df$ID) | display_df$ID %in% c("KEGG_hsa04060", "KEGG_hsa04062", "KEGG_hsa04064")
            display_df <- display_df[display_df$p.adjust < 0.05 & immune_ids & !grepl("Viral protein interaction", display_df$Description, 
                ignore.case = TRUE), , drop = FALSE]
            display_df <- display_df[order(display_df$p.adjust, -abs(display_df$NES), display_df$Description), , drop = FALSE]
            unused_df <- display_df[!display_df$ID %in% gsea_display_used, , drop = FALSE]
            selected_df <- head(unused_df, 3)
            if (nrow(selected_df) < 3) 
                selected_df <- head(rbind(selected_df, display_df[!display_df$ID %in% selected_df$ID, , drop = FALSE]), 3)
            geneSetID <- selected_df$ID
            description.grep <- selected_df$Description
            selected_df$KeyGene <- gene
            selected_df$DisplayOrder <- seq_len(nrow(selected_df))
            selected_df$SelectionRule <- "FDR < 0.05; KEGG immune-system or cytokine/chemokine/NF-kappa B signaling; ranked by FDR then absolute NES; previously displayed pathways avoided where possible"
            gsea_display_used <<- unique(c(gsea_display_used, geneSetID))
            gsea_display_selected[[gene]] <<- selected_df
        }
        else {
            geneSetID <- sortkk$ID[1:3]
            description.grep <- sortkk[sortkk$ID %in% geneSetID, ]$Description
        }
        x <- kk
        geneList <- position <- NULL
        gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
        gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL, 3)
        p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) + 
            scale_color_manual(values = mycol) + geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + ylab("Enrichment\n Score") + 
            theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "top", legend.title = element_blank(), 
            legend.background = element_rect(fill = "transparent")) + theme(axis.text.y = element_text(size = 12, face = "bold"), 
            axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), plot.margin = margin(t = 0.2, 
                r = 0.2, b = 0, l = 0.2, unit = "cm"))
        rel_heights <- c(1.5, 0.5, 1.5)
        p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + 
            ylab(NULL) + scale_color_manual(values = mycol) + theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "none", 
            plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), axis.text = element_blank(), 
            axis.line.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
        df2 <- p.res$data
        df2$y <- p.res$data$geneList[df2$x]
        df2$gsym <- p.res$data$gsym[df2$x]
        selectgenes <- data.frame(gsym = selectedGeneID)
        selectgenes <- merge(selectgenes, df2, by = "gsym")
        selectgenes <- selectgenes[selectgenes$position == 1, ]
        p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + geom_segment(data = df2, 
            aes_(x = ~x, xend = ~x, y = ~y, yend = 0), color = "grey") + geom_bar(position = "dodge", stat = "identity") + 
            scale_fill_manual(values = mycol, guide = FALSE) + scale_color_manual(values = mycol, guide = FALSE) + geom_hline(yintercept = 0, 
            lty = 2, lwd = 0.2) + ylab("Ranked list\n metric") + xlab("Rank in ordered dataset") + theme_bw() + theme(axis.text.y = element_text(size = 12, 
            face = "bold"), panel.grid = element_blank()) + geom_text_repel(data = selectgenes, show.legend = FALSE, direction = "x", 
            ylim = c(2, NA), angle = 90, size = 2.5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) + 
            theme(plot.margin = margin(t = -0.1, r = 0.4, b = 0.2, l = 0.2, unit = "cm"))
        plotlist <- list(p.res, p2, p.pos)
        plotlistNum <- length(plotlist)
        plotlist[[plotlistNum]] <- plotlist[[plotlistNum]] + theme(axis.line.x = element_line(), axis.ticks.x = element_line(), 
            axis.text.x = element_text(size = 12, face = "bold"))
        p <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
        ggsave(paste0(dir, "/", GeneIndex, ".", gene, "/", "/1.GSEA_multi_pathways.pdf"), width = 8, height = 6, plot = p)
        sortkk <- kk.gsym[kk.gsym@result$Description %in% description.grep[1] | kk.gsym@result$Description %in% description.grep[2] | 
            kk.gsym@result$Description %in% description.grep[3], ]
        go <- data.frame(Category = "KEGG", ID = sortkk$ID, Term = sortkk$Description, Genes = gsub("/", ", ", sortkk$core_enrichment), 
            adj_pval = sortkk$p.adjust)
        genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)
        df <- GOplot::circle_dat(go, genelist)[, c(3, 5, 6)]
        nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root = "all")
        edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
        nodes <- nodes %>% mutate_at(c("node.level", "node.branch"), as.character)
        graph <- tidygraph::tbl_graph(nodes, edges)
        width <- 600/(nrow(nodes) + 100)
        gc1 <- ggraph(graph, layout = "dendrogram", circular = TRUE) + geom_edge_diagonal(aes(color = node1.node.branch, 
            filter = node1.node.level != "all"), alpha = 0.5, edge_width = 2.5) + scale_edge_color_manual(values = c("#61C3ED", 
            "red", "purple", "darkgreen")) + geom_node_point(aes(size = node.size, filter = node.level != "all"), color = "#61C3ED") + 
            scale_size(range = c(width, 10)) + theme(legend.position = "none") + geom_node_text(aes(x = 1.05 * x, y = 1.05 * 
            y, label = node.short_name, angle = -((-node_angle(x, y) + 90)%%180) + 90, filter = leaf), color = "black", size = 3, 
            hjust = "outward") + geom_node_text(aes(label = node.short_name, filter = !leaf & (node.level != "all")), color = "black", 
            fontface = "bold", size = 6, family = "sans") + theme(panel.background = element_rect(fill = NA)) + coord_cartesian(xlim = c(-1.3, 
            1.3), ylim = c(-1.3, 1.3))
        ggsave(paste0(dir, "/", GeneIndex, ".", gene, "/2.Ccgraph.pdf"), width = 14, height = 14, plot = gc1)
    }
}
