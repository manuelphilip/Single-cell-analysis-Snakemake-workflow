log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)
library(devtools)
install_github('immunogenomics/presto')

seurat_obj <- readRDS(snakemake@input[["sleuth_object"]])

feature_plots <- list()
pos_markers <- list()
all_markers <- list()
top10_markers <- list()
heatmaps <- list()
vln_plots <- list()

for (i in 1: length(seurat_obj)){

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
  all_markers <- FindAllMarkers(seurat_obj[[i]])
  pos_markers[[i]] <- FindAllMarkers(seurat_obj[[i]], only.pos = TRUE)
  pos_markers[[i]] %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10_markers[[i]]
  heatmaps[[i]] <-
    DoHeatmap(seurat_obj[[i]], features = top10_markers$gene) + NoLegend()
}

pdf(file = snakemake@output[["heatmap"]])
heatmaps
dev.off()
saveRDS(seurat_obj, file = snakemake@output[["sleuth_object"]])