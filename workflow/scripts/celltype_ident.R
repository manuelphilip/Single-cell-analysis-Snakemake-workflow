log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(stringr)
library(scCustomize)
library(ggplot2)
library(tidyr)

seurat_obj <- readRDS(snakemake@input[["seurat_object"]])
sample <- snakemake@wildcards$sample

file_celltype <-
  read.csv(snakemake@input[["file_celltype"]], sep = "\t", header = FALSE)

top_10_markers <-
  read.csv(snakemake@input[["top_10_markers"]], sep = ",", header = TRUE)

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

names(new_cluster_ids) <- levels(seurat_obj[[sample]])
seurat_obj[[sample]] <- RenameIdents(seurat_obj[[sample]], new_cluster_ids)


pdf(file = snakemake@output[["dim_plot"]])
DimPlot(seurat_obj[[sample]], label = TRUE) +
  labs(x = levels((seurat_obj[[sample]]@meta.data$orig.ident)))

dev.off()
