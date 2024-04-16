log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["sleuth_object"]])

heatmaps <- list()

for (i in 1: length(seurat_obj)){
  
  pos_markers[[i]] %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10_markers[[i]]

  write.csv(top10_markers[[i]],
            file = snakemake@output[["top_10_markers"]], quote = FALSE)
  pdf(file = snakemake@output[["heatmap"]])
  heatmaps[[i]] <-
    DoHeatmap(seurat_obj[[i]], features = top10_markers$gene) + NoLegend()
  dev.off()
}
