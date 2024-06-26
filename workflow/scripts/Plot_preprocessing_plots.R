log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(scCustomize)
library(ggplot2)
model <- snakemake@params[["model"]]



seurat_obj <- readRDS(snakemake@input[["seurat_object"]])

vln_plot <- list()
variable_features <- list()
elbow_plot <- list()

for (i in 1:length(seurat_obj)) {
  vln_plot[[i]] <-
    VlnPlot(seurat_obj[[i]],
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3
    )
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seurat_obj[[i]]), 10)

  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(seurat_obj[[i]])
  plot2 <- LabelPoints(
    plot = plot1,
    points = top10, repel = TRUE
  ) +
    labs(x = levels(seurat_obj[[i]]))
  variable_features[[i]] <- plot1 + plot2
  elbow_plot[[i]] <- ElbowPlot(seurat_obj[[i]]) +
    labs(x = levels(seurat_obj[[i]]))
}

pdf(file = snakemake@output[["QC_vln_plot"]])
vln_plot
dev.off()
pdf(file = snakemake@output[["variable_features"]])
variable_features
dev.off()
pdf(file = snakemake@output[["elbow_plot"]])
elbow_plot
dev.off()
