rule perform_clustering_and_dim_reduction:
    input:
       sleuth_object="results/seurat/{model}.seurat_objt.rds",
    output:
        sleuth_object="results/seurat/clustering/{model}.seurat_objt.rds",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10
    params:
        dims=config["clustering"]["dim"],
        resolution=config["clustering"]["resolution"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/{model}.seurat_object.log",
    script:
        "../scripts/clustr-dim-red.R"


rule plot_clustering_dim_reduction_plots:
    input:
       sleuth_object="results/seurat/clustering/{model}.seurat_objt.rds",
    output:
        dim_plot="results/plots/clustering/{model}.Dim-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5
    params:
        dims=config["clustering"]["dim"],
        resolution=config["clustering"]["resolution"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/clustering/{model}.clustering.log",
    script:
        "../scripts/plot-dimplot.R"


