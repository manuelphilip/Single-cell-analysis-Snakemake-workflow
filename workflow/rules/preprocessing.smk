rule perform_preprocessing:
    input:
        samples="config/samples.tsv",
    output:
        sleuth_object="results/seurat/{model}.seurat_objt.rds",
    resources:
        cpus_per_task=20,
        mem_mb=94000
    params:
        path=config["resources"]["path"],
        clustering=config["clustering"]["activate"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/{model}.seurat_object.log",
    script:
        "../scripts/preprocessing.R"


rule plot_preprocessing_plots:
    input:
        sleuth_object="results/seurat/{model}.seurat_objt.rds",
    output:
        QC_vln_plot="results/plots/preprocessing/{model}.QC-Vln-plot.pdf",
        variable_features="results/plots/preprocessing/{model}.Highly_variable_features-plot.pdf",
        elbow_plot="results/plots/preprocessing/{model}.Elbow-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/{model}.QC-Vln-plot.log",
    script:
        "../scripts/Plot_preprocessing_plots.R"


