rule perform_preprocessing:
    input:
        samples="config/samples.tsv",
    output:
        seurat_object= "results/seurat/preprocessing/all.seurat_objt.rds",
        intergrated_seurat_object="results/seurat/preprocessing/all.seurat_integration_objt.rds",
    resources:
        cpus_per_task=20,
        runtime = 20,
        mem_mb=94000,
        nodes=10
    params:
        path=config["resources"]["path"],
        clustering=config["clustering"]["activate"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/preprocessing/all.seurat_object.log",
    script:
        "../scripts/preprocessing.R"


rule plot_preprocessing_plots:
    input:
        seurat_object="results/seurat/preprocessing/all.seurat_objt.rds",
    output:
        QC_vln_plot=report("results/plots/preprocessing/all.QC-Vln-plot.pdf",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "pre-processing Vln plot",
                },
        ),
        variable_features=report("results/plots/preprocessing/all.Highly_variable_features-plot.pdf",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "Highly variable features plot",
                },
        ),
        elbow_plot=report("results/plots/preprocessing/all.Elbow-plot.pdf",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "Elbow plot",
                },        
        )
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/plots/seurat-preprocessing/all.QC-Vln-plot.log",
    script:
        "../scripts/Plot_preprocessing_plots.R"


