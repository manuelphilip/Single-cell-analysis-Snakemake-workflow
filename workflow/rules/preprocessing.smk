rule setup_seurat_object:
    input:
        samples="config/samples.tsv",
    output:
        sleuth_object="results/sleuth/seurat_object.rds",
    params:
        path=config["resources"]["path"],
    conda:
        "../envs/seurat.yaml",
    log:
        "logs/setup_seurat_object/seurat_object/seurat_obj.log",
    script:
        "../scripts/Preprocessing.R" 



