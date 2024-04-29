log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["seurat_object"]])


dim_plot <- list()
for (i in 1: length(seurat_obj)){


  dim_plot[[i]] <- DimPlot(seurat_obj[[i]], label = TRUE) +
    labs(x = levels((seurat_obj[[i]]@meta.data$orig.ident)))

}

pdf(file = snakemake@output[["dim_plot"]])
dim_plot
dev.off()