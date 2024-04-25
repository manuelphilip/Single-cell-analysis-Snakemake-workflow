log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(stringr)
library(scCustomize)
library(ggplot2)
library(tidyr)

seurat_obj <- readRDS(snakemake@input[["intergrated_seurat_object"]])

file_celltype <-
  read.csv(snakemake@input[["file_celltype"]], sep = "\t", header = FALSE)

top_10_markers <-
  read.csv(snakemake@input[["top_10_celltype_markers"]], sep = ",", header = TRUE)

names(file_celltype) <- c("gene", "cell_type")

file_celltype$gene <- str_to_title(file_celltype$gene)
top_10_markers <-
  left_join(top_10_markers, file_celltype, by = "gene")

top_10_markers$cell_type <- top_10_markers$cell_type %>% replace_na("unknown")
cell_type_markers <- top_10_markers %>%
  group_by(cluster, cell_type) %>%
  summarise(count = n()) %>%
  filter(count == max(count)) %>%
  distinct(cluster, .keep_all = TRUE)

new_cluster_ids <- c(cell_type_markers$cell_type)


names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
seurat_obj<-AddMetaData(seurat_obj, Idents(seurat_obj), 
col.name = 'Cell_type')

pdf(file = snakemake@output[["dim_plot"]])
DimPlot(seurat_obj, label = TRUE) +
labs(x = levels((seurat_obj@meta.data$orig.ident)))
dev.off()
saveRDS(seurat_obj,
    file = snakemake@output[["celltype_seurat_object"]])