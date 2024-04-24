log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)
library(devtools)
install_github('immunogenomics/presto')

seurat_obj <- readRDS(snakemake@input[["seurat_object"]])

genes_of_interest <- c(snakemake@params[["genes_of_interest"]])
sample <- snakemake@wildcards$sample

all_markers <- FindAllMarkers(seurat_obj[[sample]])
write.csv(all_markers,
          file = snakemake@output[["all_markers"]], quote = FALSE)
  pos_markers <- FindAllMarkers(seurat_obj[[sample]], only.pos = TRUE)
  pos_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10_markers
write.csv(top10_markers,
          file = snakemake@output[["top_10_markers"]], quote = FALSE)

heatmap <-
  DoHeatmap(seurat_obj[[sample]], features = top10_markers$gene) + NoLegend()

pdf(file = snakemake@output[["heatmap"]])
heatmap
dev.off()

feature_plot <-
  FeaturePlot(seurat_obj[[sample]],
              features = genes_of_interest, label = TRUE)
pdf(file = snakemake@output[["feature_plot"]])
feature_plot
dev.off()
vln_plot <- VlnPlot(seurat_obj[[sample]], features = genes_of_interest)
pdf(file = snakemake@output[["Vln_plot"]])
vln_plot
dev.off()