log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["sleuth_object"]])

variable_features <- list()
elbow_plot <- list()

for (i in 1: length(seurat_obj)){

  seurat_obj[[i]] <- FindVariableFeatures(seurat_obj[[i]],
                                          selection.method = "vst",
                                          nfeatures = 2000)
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seurat_obj[[i]]), 10)

  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(seurat_obj[[i]])
  plot2 <- LabelPoints(plot = plot1,
                        points = top10, repel = TRUE) +
    labs(x = levels(seurat_obj[[i]]))
  variable_features[[i]] <- plot1 + plot2
  seurat_obj[[i]] <-
    RunPCA(seurat_obj[[i]], 
           features = VariableFeatures(object = seurat_obj[[i]]))
  elbow_plot[[i]] <- ElbowPlot(seurat_obj[[i]]) +
    labs(x = levels(seurat_obj[[i]]))
}

pdf(file = snakemake@output[["variable_features"]])
variable_features
dev.off()
pdf(file = snakemake@output[["elbow_plot"]])
elbow_plot
dev.off()
saveRDS(seurat_obj, file = snakemake@output[["sleuth_object"]])