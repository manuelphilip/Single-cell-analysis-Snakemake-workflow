log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["sleuth_object"]])


do_clustering <- function(obj) {
  obj <- RunUMAP(obj, dims = 1:snakemake@params[["dims"]],
                 verbose = FALSE)
  obj <- FindNeighbors(obj,
                    dims = 1:snakemake@params[["dims"]], verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE,
                      esolution = snakemake@params[["resolution"]])

}

for (i in 1:length(seurat_obj)) {

  seurat_obj[[i]] <- do_clustering(seurat_obj[[i]])
}

saveRDS(seurat_obj, file = snakemake@output[["sleuth_object"]])