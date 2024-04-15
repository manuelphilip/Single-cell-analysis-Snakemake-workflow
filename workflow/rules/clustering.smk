rule clustering:
    input:
       sleuth_object="results/seurat/preprocessing/{model}.seurat_objt.rds",
    output:
        sleuth_object="results/seurat/clustering/{model}.seurat_objt.rds",
        dim_plot="results/plots/clustering/{model}.Dim-plot.pdf",
    params:
        dims=config["clustering"]["dim"],
        resolution=config["clustering"]["resolution"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/clustering/{model}.clustering.log",
    script:
        "../scripts/clustering.R"


rule per_sam_diff_expr_and_marker_ident:
    input:
        sleuth_object="results/seurat/clustering/{model}.seurat_objt.rds",
    output:
        all_markers="results/tables/diffexp/{model}.diff-exp-genes.tsv",
        top_10_markers ="results/tables/diffexp/{model}.top-10-markers.tsv",
        feature_plot="results/plots/diffexp/{model}.Features-plot.pdf",
        heatmap="results/plots/diffexp/{model}.Heatmap-plot.pdf",
        Vln_plot="results/plots/diffexp/{model}.Top-features-Vln-plot.pdf",
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/diffexp/{model}.per_sam_diff_expr_and_marker_ident.log",
    script:
        "../scripts/Normalization_dim_red.R"
