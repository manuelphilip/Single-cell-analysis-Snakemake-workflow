rule sample_diffexp:
    input:
        intergrated_seurat_object="results/seurat/integration/{model}.seurat_objt_integration.rds",
    output:
        all_markers="results/tables/seurat/integration/{model}.diff-exp-genes.tsv",
        top_10_markers="results/tables/seurat/integration/{model}.top-10-markers.tsv",
        heatmap="results/plots/seurat/integration/{model}.Heatmap-plot.pdf",
        volcano_plot="results/plots/seurat/integration/volcano_plots/{model}.volcano_plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10
    params:
        column_name=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
        base_level=lambda wc: config["diffexp"]["models"][wc.model][
            "base_level"
        ],
        comparison_variable=lambda wc: config["diffexp"]["models"][wc.model][
            "comparison_variable"
        ],
    conda:
        "../envs/diff_exp.yaml"

    log:
        "logs/seurat/integration/sample_diffexp/{model}.sample_diffexp.log",
    script:
        "../scripts/sample-diffexp.R"

