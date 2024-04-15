rule setup_seurat_object:
    input:
        samples="config/samples.tsv",
    output:
        sleuth_object="results/seurat/{model}.seurat_objt.rds",
        QC_vln_plot="results/plots/preprocessing/{model}.QC-Vln-plot.pdf",
    params:
        path=config["resources"]["path"],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/{model}.seurat_object.log",
    script:
        "../scripts/Preprocessing.R"


rule normalization_and_dim_reduction:
    input:
        sleuth_object="results/seurat/{model}.seurat_objt.rds",
    output:
        variable_features="results/plots/preprocessing/{model}.Highly_variable_features-plot.pdf",
        elbow_plot="results/plots/preprocessing/{model}.Elbow-plot.pdf",
        sleuth_object="results/seurat/preprocessing/{model}.seurat_objt.rds",
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/preprocessing/{model}.highly_variable_features.log",
    script:
        "../scripts/Normalization_dim_red.R"
