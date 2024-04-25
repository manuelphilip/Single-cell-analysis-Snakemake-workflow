log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(scCustomize)
library(ggplot2)
library(ggrepel)


seurat_obj <- readRDS(snakemake@input[["intergrated_seurat_object"]])

column_name <- snakemake@params[["column_name"]]
base_level <- snakemake@params[["base_level"]]
comparison <- snakemake@params[["comparison_variable"]]

print(base_level)
print(comparison)

accessor_string <- paste0("seurat_obj$", column_name)
Idents(seurat_obj) <- column_name


markers <- FindMarkers(seurat_obj, ident.1 = base_level, ident.2 = comparison)



write.csv(markers,
          file = snakemake@output[["all_markers"]], quote = FALSE)
pos_markers <- FindMarkers(seurat_obj, ident.1 = base_level,
                           ident.2 = comparison, only.pos = TRUE)
pos_markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_markers
write.csv(top10_markers,
          file = snakemake@output[["top_10_markers"]], quote = FALSE)


heatmap <-
  DoHeatmap(seurat_obj, features = top10_markers$gene) + NoLegend()

pdf(file = snakemake@output[["heatmap"]])
heatmap
dev.off()

#volcano plot

markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
markers$diffexpressed[markers$avg_log2FC > 1 & markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
markers$diffexpressed[markers$avg_log2FC < -1 & markers$p_val < 0.05] <- "DOWN"



markers$gene <- row.names(markers)
markers$delabel <- NA
markers$delabel[markers$diffexpressed != "NO"] <-
  markers$gene[markers$diffexpressed != "NO"]


volcano_plot <- ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val),
                                        col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  labs(title = column_name)

pdf(file = snakemake@output[["volcano_plot"]])
volcano_plot
dev.off()