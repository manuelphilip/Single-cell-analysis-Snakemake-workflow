rule per_sam_diff_expr_and_marker_ident:
    input:
        sleuth_object="results/seurat/clustering/{model}.seurat_objt.rds",
    output:
        #all_markers="results/tables/diffexp/{model}.{sample}.diff-exp-genes.tsv",
        top_10_markers="results/tables/diffexp/{model}.{sample}.top-10-markers.tsv",
        #feature_plot="results/plots/diffexp/{model}.{sample}.Features-plot.pdf",
        heatmap="results/plots/diffexp/{model}.{sample}.Heatmap-plot.pdf",
        #Vln_plot="results/plots/diffexp/{model}.{sample}.Top-features-Vln-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000
    conda:
        "../envs/diff_exp.yaml"
    log:
        "logs/seurat/diffexp/{model}.{sample}.per_sam_diff_expr_and_marker_ident.log",
    script:
        "../scripts/diff_exp_marker_genes.R"

"""
rule diff_exp_results:
    input:
        sleuth_object="results/seurat/diffexp/{model}.seurat_objt.rds",
    output:
        #all_markers="results/tables/diffexp/{model}.{sample}.diff-exp-genes.tsv",
        top_10_markers="results/tables/diffexp/{model}.{sample}.top-10-markers.tsv",
        #feature_plot="results/plots/diffexp/{model}.{sample}.Features-plot.pdf",
        heatmap="results/plots/diffexp/{model}.{sample}.Heatmap-plot.pdf",
        #Vln_plot="results/plots/diffexp/{model}.{sample}.Top-features-Vln-plot.pdf",
        sleuth_object=sleuth_object="results/seurat/diffexp/{model}.seurat_objt.rds",
    threads:20
    conda:
        "../envs/diff_exp.yaml"
    log:
        "logs/seurat/diffexp_results/{model}.{sample}.diff_exp_results.log",
    script:
        "../scripts/diff_exp_results.R"
"""