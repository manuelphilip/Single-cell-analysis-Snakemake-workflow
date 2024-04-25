log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["intergrated_seurat_object"]])

column_name <- snakemake@params[["column_name"]]

merge_dim_plot <- DimPlot(seurat_obj,
                          reduction = "umap.unintegrated",
                          group.by = column_name)
pdf(file = snakemake@output[["merge_dim_plot"]])
merge_dim_plot
dev.off()
