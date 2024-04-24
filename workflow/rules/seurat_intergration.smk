rule perform_merge_seurat:
    input:
        intergrated_seurat_object="results/seurat/preprocessing/all.seurat_intergrated_objt.rds",
        samples="config/samples.tsv",
    output:
        merge_seurat_object="results/seurat_intergration/merge/{model}.seurat_objt_before_intergration.rds",
    resources:
        runtime = 20,
    params:
        model=get_model_samples,
        groups=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
        base_level=lambda wc: config["diffexp"]["models"][wc.model][
            "base_level"
        ],
        dims=config["merged_sample_clustering"]["dim"],
        resolution=config["merged_sample_clustering"]["resolution"],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat_intergration/merge/{model}.seurat_object_before_intergration.log",
    script:
        "../scripts/seurat_before_intergration.R"



rule plot_dim_plots_for_merged_seurat:
    input:
        seurat_object="results/seurat_intergration/merge/{model}.seurat_objt_before_intergration.rds",
    output:
        merge_dim_plot="results/plots/seurat_intergration/merge/clustering/{model}.Dim-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5
    params:
        column_name=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/plots/seurat_intergration/merge/clustering/{model}.Dim-plot.log",
    script:
        "../scripts/plot-merge-dimplot.R"

