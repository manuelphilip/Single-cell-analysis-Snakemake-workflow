log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")



library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(forcats)
library(cowplot)
library(enrichplot)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)

expression_file <-
  read.csv(snakemake@input[["diff_exp_genes"]], sep = ",", header = TRUE)

expression_file$gene <- expression_file$X


expression_file$ENSEMBL <- mapIds(org.Mm.eg.db,
                                  keys = expression_file$gene,
                                  column = "ENSEMBL",
                                  keytype = "SYMBOL",
                                  multiVals = "first")

expression_file$ENTREZID <- mapIds(org.Mm.eg.db,
                                   keys = expression_file$ENSEMBL,
                                   column = "ENTREZID",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")

ensembl_list <-
  expression_file %>%
  dplyr::filter(avg_log2FC > 0.1 | avg_log2FC < -0.1) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::select(ENSEMBL)

ego2_cc <- enrichGO(gene       = ensembl_list$ENSEMBL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = "ENSEMBL",
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

ego2_bp <- enrichGO(gene         = ensembl_list$ENSEMBL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = "ENSEMBL",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

ego2_mf <- enrichGO(gene         = ensembl_list$ENSEMBL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = "ENSEMBL",
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

go_cc_table <- data.frame(ego2_cc)
go_bp_table <- data.frame(ego2_bp)
go_mf_table <- data.frame(ego2_mf)



gene_entrezid_fc <-
  expression_file %>%
  filter(p_val < 0.05) %>%
  dplyr::select(ENTREZID, avg_log2FC)

genelist <- gene_entrezid_fc[, 2]
names(genelist) <- as.character(gene_entrezid_fc[, 1])
genelist <- sort(genelist, decreasing = TRUE)


length(genelist)

pathway <- gsePathway(genelist,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE, organism = "mouse")

pathway_table <- data.frame(pathway)



print(head(pathway_table))

write.csv(go_cc_table,
          file = snakemake@output[["GO_CC"]], quote = FALSE)
write.csv(go_bp_table,
          file = snakemake@output[["GO_BP"]], quote = FALSE)
write.csv(go_mf_table,
          file = snakemake@output[["GO_MF"]], quote = FALSE)

write.csv(pathway_table,
          file = snakemake@output[["pathway"]], quote = FALSE)


#Barplot

bar_plot_cc <- mutate(ego2_cc, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore")
bar_plot_bp <- mutate(ego2_bp, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore")
bar_plot_mf <- mutate(ego2_mf, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore")
pdf(file = snakemake@output[["bar_plot_cc"]])
bar_plot_cc
dev.off()
pdf(file = snakemake@output[["bar_plot_bp"]])
bar_plot_bp
dev.off()
pdf(file = snakemake@output[["bar_plot_mf"]])
bar_plot_mf
dev.off()

#Doptplot
if (length(ego2_cc) > 30) {
  dotpot_enrichment_cc <-
    dotplot(ego2_cc, showCategory = 30) +
    ggtitle("dotplot for GSEA Cellular component")

  dotpot_enrichment_bp <-
    dotplot(ego2_bp, showCategory = 30) +
    ggtitle("dotplot for GSEA Biological process")
  dotpot_enrichment_mf <-
    dotplot(ego2_mf, showCategory = 30) +
    ggtitle("dotplot for GSEA Molecular function")
}else {
  dotpot_enrichment_cc <-
    dotplot(ego2_cc, showCategory = length(ego2_cc)) +
    ggtitle("dotplot for GSEA Cellular component")

  dotpot_enrichment_bp <-
    dotplot(ego2_bp, showCategory = length(ego2_bp)) +
    ggtitle("dotplot for GSEA Biological process")
  dotpot_enrichment_mf <-
    dotplot(ego2_mf, showCategory = length(ego2_mf)) +
    ggtitle("dotplot for GSEA Molecular function")
}


pdf(file = snakemake@output[["dot_plot_cc"]])
dotpot_enrichment_cc
dev.off()
pdf(file = snakemake@output[["dot_plot_bp"]])
dotpot_enrichment_bp
dev.off()
pdf(file = snakemake@output[["dot_plot_mf"]])
dotpot_enrichment_mf
dev.off()

#Treeplot
ego2_cc <- pairwise_termsim(ego2_cc)
p1 <- treeplot(ego2_cc)
p2 <- treeplot(ego2_cc, hclust_method = "average")
tree_plot_cc <- aplot::plot_list(p1, p2, tag_levels = "A")

ego2_bp <- pairwise_termsim(ego2_bp)
p1 <- treeplot(ego2_bp)
p2 <- treeplot(ego2_bp, hclust_method = "average")
tree_plot_bp <- aplot::plot_list(p1, p2, tag_levels = "A")

ego2_mf <- pairwise_termsim(ego2_mf)
p1 <- treeplot(ego2_mf)
p2 <- treeplot(ego2_mf, hclust_method = "average")
tree_plot_mf <- aplot::plot_list(p1, p2, tag_levels = "A")

pdf(file = snakemake@output[["tree_plot_cc"]], height = 20, width = 30)
tree_plot_cc
dev.off()
pdf(file = snakemake@output[["tree_plot_bp"]],height = 20, width = 30)
tree_plot_bp
dev.off()
pdf(file = snakemake@output[["tree_plot_mf"]], height = 20, width = 30)
tree_plot_mf
dev.off()

#Enrichment Map
p1 <- emapplot(ego2_cc, cex_category=1.5)
p2 <- emapplot(ego2_cc, cex_category=1.5,layout="kk") 
enrichment_map_cc <-
  cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:4])



p1 <- emapplot(ego2_bp, cex_category=1.5)
p2 <- emapplot(ego2_bp, cex_category=1.5,layout="kk") 
enrichment_map_bp <-
  cowplot::plot_grid( p1, p2, ncol=2, labels=LETTERS[1:4])



p1 <- emapplot(ego2_mf, cex_category=1.5)
p2 <- emapplot(ego2_mf, cex_category=1.5,layout="kk") 
enrichment_map_mf <- 
  cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:4])


pdf(file = snakemake@output[["enrichment_map_cc"]], height = 20, width = 30)
enrichment_map_cc
dev.off()
pdf(file = snakemake@output[["enrichment_map_bp"]], height = 20, width = 30)
enrichment_map_bp
dev.off()
pdf(file = snakemake@output[["enrichment_map_mf"]], height = 20, width = 30)
enrichment_map_mf
dev.off()


#UpSet Plot
upset_plot_cc <- upsetplot(ego2_cc)
upset_plot_bp <- upsetplot(ego2_bp)
upset_plot_mf <- upsetplot(ego2_mf)

pdf(file = snakemake@output[["upset_plot_cc"]])
upset_plot_cc
dev.off()
pdf(file = snakemake@output[["upset_plot_bp"]])
upset_plot_bp
dev.off()
pdf(file = snakemake@output[["upset_plot_mf"]])
upset_plot_mf
dev.off()

#Pathway
if (nrow(pathway_table) > 30) {
  dotplot_pathway <- dotplot(pathway, showCategory = 30) +
    ggtitle("dotplot for Reactome Pathway")
  pdf(file = snakemake@output[["dotplot_pathway"]], height = 20, width = 20)
  dotplot_pathway
}else if (nrow(pathway_table) > 1) {
  dotplot_pathway <- dotplot(pathway, showCategory = nrow(pathway_table)) +
    ggtitle("dotplot for Reactome Pathway")
  pdf(file = snakemake@output[["dotplot_pathway"]], height = 20, width = 20)
  dotplot_pathway

}else {
  pdf(file = snakemake@output[["dotplot_pathway"]], height = 20, width = 20)
  dotplot_pathway <- "no significant term enriched"
  plot.new()
  text(.5, .5, dotplot_pathway, font=2, cex=1.5)
}


if (nrow(pathway_table) > 1) {
  arrange_pathway <- dplyr::arrange(pathway, abs(NES)) %>%
  dplyr::group_by(sign(NES))
  pdf(file = snakemake@output[["bar_plot_pathway"]])
  ggplot(arrange_pathway, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=10)+
    geom_col(orientation='y') +
    scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
    theme_minimal() + ylab(NULL)
}else {
  pdf(file = snakemake@output[["bar_plot_pathway"]])
  txt <- "no significant term enriched"
  plot.new()
  text(.5, .5, txt, font=2, cex=1.5)
}
