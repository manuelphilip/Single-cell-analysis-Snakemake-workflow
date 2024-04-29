log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)

model <- snakemake@params[["model"]]

expression_matrices <-
  Read10X_Multi_Directory(base_path = snakemake@params[["path"]])

samples <- read.csv(snakemake@input[["samples"]], header = TRUE, sep = "\t")


do_preprocess <- function(obj) {
  obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
}


seurat_obj <- list() # create an empty seurat_object list
intr_seurat_obj <- list()

for (i in 1:length(samples[[1]])) {
  sam_seurat_objt <- names(expression_matrices[i])
  sample_name <- names(expression_matrices[i])
  sam_seurat_objt <-
    CreateSeuratObject(counts = expression_matrices[i], project = sample_name)
  intr_seurat_obj[[sample_name]] <-
    PercentageFeatureSet(sam_seurat_objt,
                         pattern = "^mt-", col.name = "percent.mt")
  seurat_obj[[sample_name]] <- do_preprocess(sam_seurat_objt)
}


saveRDS(intr_seurat_obj, file = snakemake@output[["intergrated_seurat_object"]])
saveRDS(seurat_obj, file = snakemake@output[["seurat_object"]])
