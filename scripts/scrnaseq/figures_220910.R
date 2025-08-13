#!/usr/bin/env Rscript

# load necessary libraries
library(Seurat)
library(SeuratDisk)
library(tidyverse)


path_lib.map = rbind(data.frame(lib = paste0("L", 5:22),
                                counts = paste0("../../2022_05_05_NovaSeq_L5-L22/count/tk_2022_05_05_L", 5:22),
                                demux = paste0("../../2022_05_05_NovaSeq_L5-L22/popscle/L", 5:22, "/L", 5:22, "_demu.best")),
                     data.frame(lib = paste0("L", 23:43),
                                counts = paste0("../../Ning_NovASeq_L23-L43_0812/L23-L43_count/HWTLLDSX3_L", 23:43),
                                demux = paste0("../../Ning_NovASeq_L23-L43_0812/L23-L43_demuxlet/L", 23:43, "/L", 23:43, "_demu.best")),
                     data.frame(lib = paste0("L", 44:64),
                                counts = paste0("../../Ning_NovASeq_L44-L64_0822/L44-L64_count/HYG2CDSX3_mL", 44:64),
                                demux = paste0("../../Ning_NovASeq_L44-L64_0822/L44-L64_demuxlet/L", 44:64, "/L", 44:64, "_demu.best")),
                     data.frame(lib = paste0("L", 65:85),
                                counts = paste0("../../Ning_NovASeq_L65-85_0829/L65-L85_count/HYCWWDSX3_mL", 65:85),
                                demux = paste0("../../Ning_NovASeq_L65-85_0829/L65-L85_demuxlet/L", 65:85, "/L", 65:85, "_demu.best")))


lib_samp_cond.map1 = read_tsv("../../sample_info/10xfang_topmed.ID.tsv",
                              col_names = c("lib", paste0("foo", 1:5), paste0("samp_", 1:8))) %>%
    select(-starts_with("foo")) %>%
    gather(key = samp_n, value = samp_cond, -lib) %>%
    select(-samp_n) %>%
    separate(samp_cond, into = c("samp", "cond"), sep = "_")

lib_samp_cond.map2 = read_tsv("../../sample_info/10xsample.topmed.id.tsv",
                             col_names = c("lib", paste0("foo", 1:8), paste0("samp_", 1:8))) %>%
    select(-starts_with("foo")) %>%
    gather(key = samp_n, value = samp_cond, -lib) %>%
    select(-samp_n) %>%
    separate(samp_cond, into = c("samp", "cond"), sep = "_")

lib_samp_cond.map = rbind(lib_samp_cond.map1, lib_samp_cond.map2)
lib_samp_cond.map = lib_samp_cond.map[complete.cases(lib_samp_cond.map),]


meta_info = read_tsv("../../sample_info/10x2topmed.id.tsv",
                     col_names = c("field_id", "lab_id", "sex", "age", "pop", "topmed_id"))


info = path_lib.map[1,]

cells.demux = read_tsv(info$demux) %>%
    filter(DROPLET.TYPE == "SNG") %>%
    select(BARCODE, BEST.GUESS) %>%
    separate(BEST.GUESS, into = c("samp", "foo", "pp")) %>%
    select(barcode = BARCODE, samp)

cells.demux = cells.demux[sort(sample(nrow(cells.demux), 10000)),]
    

cells.raw = Read10X(paste0(info$counts, "/outs/filtered_feature_bc_matrix")) %>%
    CreateSeuratObject %>%
    .[,cells.demux$barcode]

cells.info = data.frame(barcode = colnames(cells.raw)) %>%
    merge(cells.demux) %>%
    mutate(lib = info$lib) %>%
    merge(lib_samp_cond.map) %>%
    arrange(barcode)

cells.raw$lib = info$lib
cells.raw$samp = cells.info$samp
cells.raw$cond = cells.info$cond

cells.merge <- cells.raw


## read in count data
for (i in 3:nrow(path_lib.map)){

    info = path_lib.map[i,]

    cells.demux = read_tsv(info$demux) %>%
        filter(DROPLET.TYPE == "SNG") %>%
        select(BARCODE, BEST.GUESS) %>%
        separate(BEST.GUESS, into = c("samp", "foo", "pp")) %>%
        select(barcode = BARCODE, samp)

    if (nrow(cells.demux) >= 10000){
        cells.demux = cells.demux[sort(sample(nrow(cells.demux), 10000)),]
    }

    cells.raw = Read10X(paste0(info$counts, "/outs/filtered_feature_bc_matrix")) %>%
        CreateSeuratObject %>%
        .[,cells.demux$barcode]

    cells.info = data.frame(barcode = colnames(cells.raw)) %>%
        merge(cells.demux) %>%
        mutate(lib = info$lib) %>%
        merge(lib_samp_cond.map) %>%
        arrange(barcode)

    cells.raw$lib = info$lib
    cells.raw$samp = cells.info$samp
    cells.raw$cond = cells.info$cond
    
    cells.merge <- merge(cells.merge, y = cells.raw)

}


## features <- SelectIntegrationFeatures(object.list = cells.list)
## cells.list <- lapply(X = cells.list, FUN = function(x) {
##     x <- ScaleData(x, features = features, verbose = FALSE)
##     x <- RunPCA(x, features = features, verbose = FALSE)
## })

## anchors <- FindIntegrationAnchors(object.list = cells.list, reference = c(60),
##                                   reduction = "rpca", dims = 1:50)
## cells.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

## cells.integrated <- ScaleData(cells.integrated, verbose = FALSE)
## cells.integrated <- RunPCA(cells.integrated, verbose = FALSE)
## cells.integrated <- RunUMAP(cells.integrated, dims = 1:50)


SaveH5Seurat(cells.merge, overwrite = TRUE)


cell_sub = data.frame(cell_ids = colnames(cells.merge),
                      lib = cells.merge$lib) %>%
    group_by(lib) %>%
    sample_n(2000)

cells.merge = cells.merge[,cell_sub$cell_ids]

## cells.merge <- NormalizeData(cells.merge, verbose = FALSE)
## cells.merge <- FindVariableFeatures(cells.merge, verbose = FALSE)
## cells.merge <- ScaleData(cells.merge, verbose = FALSE)
## cells.merge <- RunPCA(cells.merge, verbose = FALSE)
## cells.merge <- RunUMAP(cells.merge, dims = 1:50)


cells.merge <- LoadH5Seurat("./tmp_data/SeuratProject.h5Seurat")


cells.merge.info = data.frame(barcode = colnames(cells.merge),
                              topmed_id = cells.merge$samp) %>%
    merge(meta_info) %>%
    select(-topmed_id) %>%
    column_to_rownames("barcode")

cells.merge <- AddMetaData(
    object = cells.merge,
    metadata = cells.merge.info)


## DimPlot(cells.merge, group.by = "cond")
## DimPlot(cells.merge, group.by = "sex")
## DimPlot(cells.merge, group.by = "age")
## DimPlot(cells.merge, group.by = "pop")


## ElbowPlot(cells.merge, 1:100)

## ## ~20 dimensions

cells.merge <- FindNeighbors(cells.merge, dims = 1:20)
cells.merge <- FindClusters(cells.merge, resolution = 0.5)


## select 500 cells per cluster
clust_cells = data.frame(barcode = colnames(cells.merge),
           clust = cells.merge$seurat_clusters) %>%
    group_by(clust) %>%
    sample_n(500)

jimmie.markers = c("CD14","S100A9","CTSS","FCGR3A", "CDKN1C", "FCER1A", "CD1C", "HLA-DQA1",
                   "JCHAIN", "PLD4", "SERPINF1", "LILRA4", "PRSS57", "SOX4", "STMN1", "TOP2A",
                   "CD3E", "CD3D", "CD3G", "CD3A", "CD4", "CD8A", "CD8B2", "CD8B", "CD40LG", "CCR7",
                   "IL7R", "TNFRSF4", "RTKN2", "FOXP3", "PRF1", "CCL5", "GZMH", "GZMB", "GZMK", "KLRB1", "GNLY",
                   "PDCD1", "CTLA4", "LAG3", "CD16", "TIGIT", "HAVCR2", "CD244", "NKG7", "CD79A", "MZB1", "IGLL5",
                   "BANK1", "CD19", "MS4A1")

DoHeatmap(cells.merge[,clust_cells$barcode], features = jimmie.markers)


## ## add library groups
## group1 = paste0("L", 5:22)
## group2 = paste0("L", 23:43)
## group3 = paste0("L", 44:64)
## group4 = paste0("L", 65:85)

## group.df = rbind(data.frame(lib = group1, grp = "Fang"),
##                  data.frame(lib = group2, grp = "Ning/Yuanqing 1"),
##                  data.frame(lib = group3, grp = "Ning/Yuanqing 2"),
##                  data.frame(lib = group4, grp = "Ning/Yuanqing 3"))

## group.info = data.frame(barcode = colnames(cells.merge),
##                         lib = cells.merge$lib) %>%
##     merge(group.df) %>%
##     column_to_rownames("barcode") %>%
##     select(-lib)

## cells.merge <- AddMetaData(
##     object = cells.merge,
##     metadata = group.info)


## celltype map
ctype.map = rbind(data.frame(clust = c(0, 1, 2, 5, 6, 8, 12, 16, 17),
                             ctype = "T"),
                  data.frame(clust = c(3, 4, 7, 15),
                             ctype = "NK"),
                  data.frame(clust = c(9, 11, 13),
                             ctype = "B"),
                  data.frame(clust = c(10, 14, 18),
                             ctype = "Monocytes"))


cells.merge.mono.fulani = subset(cells.merge, subset = (seurat_clusters %in% c(10, 14, 18)) & (pop == "Fulani"))
cells.merge.mono.tikari = subset(cells.merge, subset = (seurat_clusters %in% c(10, 14, 18)) & (pop == "TIKARI"))
cells.merge.mono.baka = subset(cells.merge, subset = (seurat_clusters %in% c(10, 14, 18)) & (pop == "BAKA"))


cells.merge.lym.fulani = subset(cells.merge, subset = (seurat_clusters %in% c(0:9, 11:13, 15:17)) & (pop == "Fulani"))
cells.merge.lym.tikari = subset(cells.merge, subset = (seurat_clusters %in% c(0:9, 11:13, 15:17)) & (pop == "TIKARI"))
cells.merge.lym.baka = subset(cells.merge, subset = (seurat_clusters %in% c(0:9, 11:13, 15:17)) & (pop == "BAKA"))


cells.merge.b.fulani = subset(cells.merge, subset = (seurat_clusters %in% c(9,11,13)) & (pop == "Fulani"))
cells.merge.b.tikari = subset(cells.merge, subset = (seurat_clusters %in% c(9,11,13)) & (pop == "TIKARI"))
cells.merge.b.baka = subset(cells.merge, subset = (seurat_clusters %in% c(9,11,13)) & (pop == "BAKA"))



Idents(cells.merge.mono.fulani) = "cond"
cells.merge.mono.fulani.lps = FindMarkers(cells.merge.mono.fulani, ident.1 = "LPS", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.mono.tikari) = "cond"
cells.merge.mono.tikari.lps = FindMarkers(cells.merge.mono.tikari, ident.1 = "LPS", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.mono.baka) = "cond"
cells.merge.mono.baka.lps = FindMarkers(cells.merge.mono.baka, ident.1 = "LPS", ident.2 = "Ctrl", logfc.threshold = 0)


Idents(cells.merge.mono.fulani) = "cond"
cells.merge.mono.fulani.ifn = FindMarkers(cells.merge.mono.fulani, ident.1 = "IFN", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.mono.tikari) = "cond"
cells.merge.mono.tikari.ifn = FindMarkers(cells.merge.mono.tikari, ident.1 = "IFN", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.mono.baka) = "cond"
cells.merge.mono.baka.ifn = FindMarkers(cells.merge.mono.baka, ident.1 = "IFN", ident.2 = "Ctrl", logfc.threshold = 0)


Idents(cells.merge.lym.fulani) = "cond"
cells.merge.lym.fulani.lps = FindMarkers(cells.merge.lym.fulani, ident.1 = "LPS", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.lym.tikari) = "cond"
cells.merge.lym.tikari.lps = FindMarkers(cells.merge.lym.tikari, ident.1 = "LPS", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.lym.baka) = "cond"
cells.merge.lym.baka.lps = FindMarkers(cells.merge.lym.baka, ident.1 = "LPS", ident.2 = "Ctrl", logfc.threshold = 0)


Idents(cells.merge.b.fulani) = "cond"
cells.merge.b.fulani.ifn = FindMarkers(cells.merge.b.fulani, ident.1 = "IFN", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.b.tikari) = "cond"
cells.merge.b.tikari.ifn = FindMarkers(cells.merge.b.tikari, ident.1 = "IFN", ident.2 = "Ctrl", logfc.threshold = 0)

Idents(cells.merge.b.baka) = "cond"
cells.merge.b.baka.ifn = FindMarkers(cells.merge.b.baka, ident.1 = "IFN", ident.2 = "Ctrl", logfc.threshold = 0)


mono.de.df = cells.merge.mono.tikari.lps %>%
    rownames_to_column("gene_name") %>%
    select(gene_name, tikari_fc = avg_log2FC) %>%
    merge(cells.merge.mono.fulani.lps %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, fulani_fc = avg_log2FC)) %>%
    merge(cells.merge.mono.baka.lps %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, baka_fc = avg_log2FC)) %>%
    gather(key = pop, value = log_fc, -c(gene_name, tikari_fc))


mono.ifn.de.df = cells.merge.mono.tikari.ifn %>%
    rownames_to_column("gene_name") %>%
    select(gene_name, tikari_fc = avg_log2FC) %>%
    merge(cells.merge.mono.fulani.ifn %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, fulani_fc = avg_log2FC)) %>%
    merge(cells.merge.mono.baka.ifn %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, baka_fc = avg_log2FC)) %>%
    gather(key = pop, value = log_fc, -c(gene_name, tikari_fc))


lym.de.df = cells.merge.lym.tikari.lps %>%
    rownames_to_column("gene_name") %>%
    select(gene_name, tikari_fc = avg_log2FC) %>%
    merge(cells.merge.lym.fulani.lps %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, fulani_fc = avg_log2FC)) %>%
    merge(cells.merge.lym.baka.lps %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, baka_fc = avg_log2FC))

b.ifn.de.df = cells.merge.b.tikari.ifn %>%
    rownames_to_column("gene_name") %>%
    select(gene_name, tikari_fc = avg_log2FC) %>%
    merge(cells.merge.b.fulani.ifn %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, fulani_fc = avg_log2FC)) %>%
    merge(cells.merge.b.baka.ifn %>%
          rownames_to_column("gene_name") %>%
          select(gene_name, baka_fc = avg_log2FC)) %>%
    gather(key = pop, value = log_fc, -c(gene_name, tikari_fc))




ggplot(mono.de.df, aes(x = tikari_fc, y = log_fc)) +
    geom_point(alpha = 0.75) +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    theme_bw() +
    geom_label_repel(data = ful_de_labels, aes(x = tikari_fc, y = log_fc, label = gene_name), max.overlaps = 10000) +
    facet_grid(. ~ pop)

ggplot(mono.ifn.de.df, aes(x = tikari_fc, y = log_fc)) +
    geom_point(alpha = 0.75) +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    theme_bw(base_size = 20) +
    facet_grid(. ~ pop)


ggplot(b.ifn.de.df, aes(x = tikari_fc, y = log_fc)) +
    geom_point(alpha = 0.75) +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    theme_bw(base_size = 20) +
    facet_grid(. ~ pop)


DimPlot(cells.merge[, cells.merge$pop == "TIKARI"], split.by = "cond", group.by = "cond")




library(gprofiler2)

mono.de.go = gost(query = mono.de.df %>% arrange(-abs(tikari_fc - fulani_fc)) %>% pull(gene_name),
                  custom_bg = mono.de.df$gene_name,
                  organism = "hsapiens",
                  ordered_query = T)
