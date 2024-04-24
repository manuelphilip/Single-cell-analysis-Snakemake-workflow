rule perform_clustering_and_dim_reduction:
    input:
       seurat_object="results/seurat/preprocessing/all.seurat_objt.rds",
    output:
        seurat_object="results/seurat/clustering/all.seurat_objt.rds",
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
        "logs/seurat/clustering/all.seurat_object.log",
    script:
        "../scripts/clustr-dim-red.R"


rule plot_clustering_dim_reduction_plots:
    input:
       seurat_object="results/seurat/clustering/all.seurat_objt.rds",
    output:
        dim_plot="results/plots/clustering/all.Dim-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5
    params:
        dims=config["clustering"]["dim"],
        resolution=config["clustering"]["resolution"],
        
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/clustering/all.clustering.log",
    script:
        "../scripts/plot-dimplot.R"


