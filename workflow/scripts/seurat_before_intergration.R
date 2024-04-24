log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)
library(tibble)

model <- snakemake@params[["model"]]

print(model)

print(str(model))
samples <- read.csv(snakemake@input[["samples"]], header = TRUE, sep = "\t")

print(snakemake@params[["groups"]])
print(snakemake@params[["base_level"]])

intergrated_seurat_object <-
  readRDS(snakemake@input[["intergrated_seurat_object"]])

merged_seurat_object <- list()

do_seurat <- function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj,
                       dims = 1:snakemake@params[["dims"]], verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE,
                      resolution = snakemake@params[["resolution"]])
  obj <- RunUMAP(obj, dims = 1:snakemake@params[["dims"]],
                 verbose = FALSE, reduction = "pca",
                 reduction.name = "umap.unintegrated")
}



merge_seurat_objects <- function(x, y) {
  # Combine objects (modify for your specific merging logic)
  merged_seurat <- merge(x, y)
  return(merged_seurat)
}

merged_seurat <- Reduce(merge_seurat_objects, intergrated_seurat_object)

merged_seurat <- JoinLayers(merged_seurat)

samples <- rename(samples,  orig.ident = sample)

groups <- snakemake@params[["groups"]]

accessor_string <- paste0("merged_seurat$", groups)

merged_seurat@meta.data <-
  merged_seurat@meta.data %>%
  rownames_to_column("cells") %>%
  full_join(samples, by = "orig.ident") %>%
  column_to_rownames("cells")
merged_seurat[["RNA"]] <-
  split(merged_seurat[["RNA"]], f = eval(parse(text = accessor_string)))

merged_seurat <- do_seurat(merged_seurat)

saveRDS(merged_seurat, file = snakemake@output[["merge_seurat_object"]])




















