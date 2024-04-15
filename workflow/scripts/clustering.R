log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["sleuth_object"]])

dim_plot <- list()

for (i in 1: length(seurat_obj)){

  seurat_obj[[i]] <- RunUMAP(seurat_obj[[i]],
                             dims = 1:snakemake@params[["dims"]],
                             verbose = FALSE)
  seurat_obj[[i]] <-
    FindNeighbors(seurat_obj[[i]], dims = 1:snakemake@params[["dims"]], verbose = FALSE)
  seurat_obj[[i]] <-
    FindClusters(seurat_obj[[i]], verbose = FALSE,
                 resolution = snakemake@params[["resolution"]])
  dim_plot[[i]] <- DimPlot(seurat_obj[[i]], label = TRUE) +
    labs(x = levels((seurat_obj[[i]]@meta.data$orig.ident)))

}

pdf(file = snakemake@output[["dim_plot"]])
dim_plot
dev.off()
saveRDS(seurat_obj, file = snakemake@output[["sleuth_object"]])