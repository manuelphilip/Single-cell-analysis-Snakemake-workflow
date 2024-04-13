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

# samples<-read.csv("sample.csv", header = T,sep = "\t")
samples <- read.csv(snakemake@input[["samples"]], header = TRUE, sep = "\t")
print(samples[[1]])

do_preprocess <- function(obj) {
  obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  obj
}

seurat_obj <- list() # create an empty seurat_object list
vln_plot <- list()
for (i in 1:length(samples[[1]])) {
  sam_seurat_objt <- names(expression_matrices[i])
  sample_name <- names(expression_matrices[i])
  print(sample_name)
  sam_seurat_objt <-
    CreateSeuratObject(counts = expression_matrices[i], project = sample_name)
  seurat_obj[[sample_name]] <- do_preprocess(sam_seurat_objt)
  vln_plot[[i]] <-
    VlnPlot(seurat_obj[[i]],
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3
    )
}
pdf(file = snakemake@output[["QC_vln_plot"]])
vln_plot
dev.off()
saveRDS(seurat_obj, file = snakemake@output[["sleuth_object"]])
