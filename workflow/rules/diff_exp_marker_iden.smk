if config["visualize_marker_expression"]["activate"]:
    rule per_sample_diff_expression_and_marker_identification:
        input:
            seurat_object="results/seurat/clustering/all.seurat_objt.rds",
        output:
            all_markers="results/tables/diffexp/{sample}.diff-exp-genes.tsv",
            top_10_markers="results/tables/diffexp/{sample}.top-10-markers.tsv",
            heatmap="results/plots/diffexp/{sample}.Heatmap-plot.pdf",
            feature_plot="results/plots/diffexp/{sample}.Features-plot.pdf",
            Vln_plot="results/plots/diffexp/{sample}.Top-features-Vln-plot.pdf",
        resources:
            cpus_per_task=20,
            mem_mb=94000,
            nodes=10
        params:
            genes_of_interest=config["visualize_marker_expression"]["genes_of_interest"],
        conda:
            "../envs/diff_exp.yaml"
        log:
            "logs/seurat/diffexp/{sample}.per_sam_diff_expr_and_marker_ident.log",
        script:
            "../scripts/diff_exp_marker_genes.R"
else:
    rule per_sample_diff_expression_and_marker_identification:
        input:
            seurat_object="results/seurat/clustering/all.seurat_objt.rds",
        output:
            all_markers="results/tables/diffexp/{sample}.diff-exp-genes.tsv",
            top_10_markers="results/tables/diffexp/{sample}.top-10-markers.tsv",
            heatmap="results/plots/diffexp/{sample}.Heatmap-plot.pdf",
        resources:
            cpus_per_task=20,
            mem_mb=94000,
            nodes=10
        conda:
            "../envs/diff_exp.yaml"
        log:
            "logs/seurat/diffexp/{sample}.per_sam_diff_expr_and_marker_ident.log",
        script:
            "../scripts/diff_exp_marker_genes.R"

rule assign_cell_type_identity_to_cluster:
    input:
        seurat_object="results/seurat/clustering/all.seurat_objt.rds",
        file_celltype=config["celltype_annotation"]["path"],
        top_10_markers="results/tables/diffexp/{sample}.top-10-markers.tsv",
    output:
        dim_plot="results/plots/celltype/{sample}.Dim-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/plots/celltype/{sample}.dimplot.log",
    script:
        "../scripts/celltype_ident.R"