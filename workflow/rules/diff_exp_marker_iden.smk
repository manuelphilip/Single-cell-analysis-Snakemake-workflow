if config["visualize_marker_expression"]["activate"]:
    rule per_sample_diff_expression_and_marker_identification:
        input:
            seurat_object="results/seurat/clustering/all.seurat_objt.rds",
        output:
            all_markers=report("results/tables/diffexp/{sample}.diff-exp-genes.tsv",
                    caption="../report/sample_all_markers.rst",
                    category="Per sample marker identification",  
                    subcategory="Differential expression tables",
                    labels={
                        "sample":"{sample}", "table": "differentially expressed genes",
                    },
            ),
            top_10_markers=report("results/tables/diffexp/{sample}.top-10-markers.tsv",
                    caption="../report/sample_top_10_markers.rst",
                    category="Per sample marker identification",  
                    subcategory="Differential expression tables",
                    labels={
                        "sample":"{sample}", "table": "top-10-markers",
                    },
            ),
            heatmap=report("results/plots/diffexp/{sample}.Heatmap-plot.pdf",
                    caption="../report/heatmap.rst",
                    category="Per sample marker identification",  
                    subcategory="per sample Heatmap",
                    labels={
                        "sample":"{sample}", "plot":"Heatmap",
                    },            
            
            ),
            feature_plot=report("results/plots/diffexp/{sample}.Features-plot.pdf",
                    caption="../report/sample_features_plot.rst",
                    category="Per sample marker identification",  
                    subcategory="per sample Top features plot",
                    labels={
                        "sample":"{sample}", "plot":"Top features plot",
                    },            
            ),
            Vln_plot=report("results/plots/diffexp/{sample}.Top-features-Vln-plot.pdf",
                    caption="../report/sample_volin_plot.rst",
                    category="Per sample marker identification",  
                    subcategory="per sample Top features Vln plot",
                    labels={
                        "sample":"{sample}", "plot":"Top features Vln plot",
                    },            
                ),
        resources:
            cpus_per_task=20,
            mem_mb=94000,
            nodes=10
        params:
            genes_of_interest=config["visualize_marker_expression"]["genes_of_interest"],
        conda:
            "../envs/diff_exp.yaml",
        log:
            "logs/seurat/diffexp/{sample}.per_sam_diff_expr_and_marker_ident.log",
        script:
            "../scripts/diff_exp_marker_genes.R"
else:
    rule per_sample_diff_expression_and_marker_identification:
        input:
            seurat_object="results/seurat/clustering/all.seurat_objt.rds",
        output:
            all_markers=report("results/tables/diffexp/{sample}.diff-exp-genes.tsv",
                caption="../report/sample_all_markers.rst",
                category="Per sample marker identification",
                subcategory="Differential expression tables",
                labels={
                    "sample":"{sample}", "table":"differentially expressed genes",
                },  
            
            ),
            top_10_markers=report("results/tables/diffexp/{sample}.top-10-markers.tsv",
                caption="../report/sample_top_10_markers.rst",
                category="Per sample marker identification",
                subcategory="Differential expression tables",
                labels={
                    "sample":"{sample}", "table":"top 10 markers",
                },  
            ),
            heatmap=report("results/plots/diffexp/{sample}.Heatmap-plot.pdf",
                caption="../report/heatmap.rst",
                category="Per sample marker identification",
                subcategory="per sample heatmap",
                labels={
                    "sample":"{sample}", "plot":"Heatmap",
                },             
            
            )
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
        dim_plot=report("results/plots/celltype/{sample}.Dim-plot.pdf",
            caption="../report/sample_dimplot.rst",
            category="Per sample assigned celltype to Dimensions",
            subcategory="celltype Dimensions plot",
            labels={
                "sample":"{sample}", "plot":"Dim plot",
            }, 
        
        
        )
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