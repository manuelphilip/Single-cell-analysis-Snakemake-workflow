log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)
library(tibble)
library(harmony)
library(devtools)
install_github('immunogenomics/presto')

model <- snakemake@params[["model"]]


merge_seurat_object <-
  readRDS(snakemake@input[["merge_seurat_object"]])

print(merge_seurat_object)

do_integration <- function(obj) {
  obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                         orig.reduction = "pca", new.reduction = "harmony",
                         assay = "RNA", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony",
                       dims = 1:snakemake@params[["dims"]])
  obj <- FindClusters(obj, resolution = snakemake@params[["resolution"]],
                      cluster.name = "harmony_clusters")
  obj <- RunUMAP(obj, reduction = "harmony",
                 dims = 1:snakemake@params[["dims"]],
                 reduction.name = "umap.harmony")
}

merge_seurat_object <- do_integration(merge_seurat_object)

merge_seurat_object <- JoinLayers(merge_seurat_object)

pos_markers <- FindAllMarkers(merge_seurat_object, only.pos = TRUE)
pos_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_markers
write.csv(top10_markers,
          file = snakemake@output[["top_10_markers"]], quote = FALSE)


saveRDS(merge_seurat_object,
        file = snakemake@output[["intergrated_seurat_object"]])
