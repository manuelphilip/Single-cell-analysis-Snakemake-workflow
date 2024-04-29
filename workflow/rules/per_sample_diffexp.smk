rule sample_diffexp:
    input:
        intergrated_seurat_object="results/seurat/integration/{model}.seurat_objt_integration.rds",
    output:
        all_markers=report("results/tables/seurat/integration/{model}.diff-exp-genes.tsv",
            caption="../report/model_all_markers.rst",
            category="Per model differential expression",
            subcategory="differential expression tables",
            labels={
                   "model": "{model}", "table": "differentially expressed genes"
                },  
        
        ),
        top_10_markers=report("results/tables/seurat/integration/{model}.top-10-markers.tsv",
            caption="../report/model_top_10_markers.rst",
            category="Per model differential expression",
            subcategory="differential expression tables",
            labels={
                   "model": "{model}", "table": "top 10 markers"
                },  
        ),
        volcano_plot=report("results/plots/seurat/integration/volcano_plots/{model}.volcano_plot.pdf",
            caption="../report/model_volcano_plot.rst",
            category="Per model differential expression",
            subcategory="Volcano plot",
            labels={
                   "model": "{model}", "plot": "Volcano plot",
                },  
        )
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

