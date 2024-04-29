rule per_sample_enrichment_analysis:
    input:
        diff_exp_genes="results/tables/seurat/integration/{model}.diff-exp-genes.tsv",
    output:
        GO_CC=report("results/tables/seurat/integration/enrichment/groups/{model}.top_cellular_components.tsv",
                caption="../report/sample_enrichment_table_cc.rst",
                category="Per model enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "table": "Significantly enriched cellular components",
                },  
        ),
        GO_BP=report("results/tables/seurat/integration/enrichment/groups/{model}.top_biological_process.tsv",
                caption="../report/sample_enrichment_table_bp.rst",
                category="Per model enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "table": "Significantly enriched Biological process",
                },  
        
        ),
        GO_MF=report("results/tables/seurat/integration/enrichment/groups/{model}.top_molecular_function.tsv",
                caption="../report/sample_enrichment_table_mf.rst",
                category="Per model enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "table": "Significantly enriched Molecular functions",
                },         
        
        ),
        pathway=report("results/tables/seurat/integration/pathway/groups/{model}.pathway_results.tsv",
                caption="../report/sample_pathway_table.rst",
                category="Per model pathway analysis",
                subcategory="tables",
                labels={
                   "model": "{model}", "table": "Significantly enriched pathways",
                },          
        
        ),
        bar_plot_cc= report("results/plots/seurat/integration/enrichment/groups/{model}.bar_plot_cc.pdf",
                category="Per model enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "plot": "Cellular component bar plot",
                },  
        ),
        bar_plot_bp=report("results/plots/seurat/integration/enrichment/groups/{model}.bar_plot_bp.pdf",
                category="Per model enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "plot": "Biological process bar plot",
                },          
        ),

        bar_plot_mf=report("results/plots/seurat/integration/enrichment/groups/{model}.bar_plot_mf.pdf",
                category="Per model enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "plot": "Molecular function bar plot",
                },        
        ),

        dot_plot_cc=report("results/plots/seurat/integration/enrichment/groups/{model}.dot_plot_cc.pdf",
                category="Per model enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "plot": "Cellular component dot plot",
                },         
        ),

        dot_plot_bp=report("results/plots/seurat/integration/enrichment/groups/{model}.dot_plot_bp.pdf",
                category="Per model enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "plot": "Biological process dot plot",
                },         
        ),

        dot_plot_mf=report("results/plots/seurat/integration/enrichment/groups/{model}.dot_plot_mf.pdf",
                category="Per model enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "plot": "Molecular function dot plot",
                },           
        ),

        tree_plot_cc=report("results/plots/seurat/integration/enrichment/groups/{model}.tree_plot_cc.pdf",
                caption="../report/sample_tree_plot.rst",
                category="Per model enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "plot": "Cellular component tree plot",
                },          
        ),
        tree_plot_bp=report("results/plots/seurat/integration/enrichment/groups/{model}.tree_plot_bp.pdf",
               caption="../report/sample_tree_plot.rst",
               category="Per model enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "plot": "Biological process tree plot",
                },         
        ),
        tree_plot_mf=report("results/plots/seurat/integration/enrichment/groups/{model}.tree_plot_mf.pdf",
                caption="../report/sample_tree_plot.rst",
                category="Per model enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "plot": "Molecular function tree plot",
                },          
        ),
        enrichment_map_cc=report("results/plots/seurat/integration/enrichment/groups/{model}.enrichment_map_cc.pdf",
                caption="../report/sample_enrichment_map.rst",
                category="Per model enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "plot": "Cellular component enrichment map",
                },          
        
        ),
        enrichment_map_bp=report("results/plots/seurat/integration/enrichment/groups/{model}.enrichment_map_bp.pdf",
               caption="../report/sample_enrichment_map.rst",
               category="Per model enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "plot": "Biological process enrichment map",
                }, 
        ),
        enrichment_map_mf=report("results/plots/seurat/integration/enrichment/groups/{model}.enrichment_map_mf.pdf",
                caption="../report/sample_enrichment_map.rst",
                category="Per model enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "plot": "Molecular function enrichemt map",
                },         
        ),
        upset_plot_cc=report("results/plots/seurat/integration/enrichment/groups/{model}.upset_plot_cc.pdf",
                category="Per model enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "plot": "Cellular component upset plot",
                }, 
        
        ),
        upset_plot_bp=report("results/plots/seurat/integration/enrichment/groups/{model}.upset_plot_bp.pdf",
               category="Per model enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "plot": "Biological process upset plot",
                },
        
        ),
        upset_plot_mf=report("results/plots/seurat/integration/enrichment/groups/{model}.upset_plot_mf.pdf",
                category="Per model enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "plot": "Molecular function upset plot",
                },  
        
        ),
        dotplot_pathway=report("results/plots/seurat/integration/pathway/groups/{model}.dotplot_pathway.pdf",
                category="Per model pathway analysis",
                subcategory="dotplot",
                labels={
                   "model": "{model}", "plot": "dotplot",
                },  

        ),
        bar_plot_pathway=report("results/plots/seurat/integration/pathway/groups/{model}.bar_plot_pathway.pdf",
                caption="../report/sample_pathway_barplot.rst",
                category="Per model pathway analysis",
                subcategory="barplot",
                labels={
                   "model": "{model}", "plot": "barplot",
                },  
        ),
    params:
      gene_ontology_p_val=config["enrichment"]["sig-level"]["gene_ontology_p_val"],
      gene_ontology_q_val=config["enrichment"]["sig-level"]["gene_ontology_q_val"],
      pathway_significant=config["enrichment"]["sig-level"]["pathway"],
    resources:
        runtime = 20,
    conda:
        "../envs/enrichment.yaml",
    log:
        "logs/seurat/integration/per_sample_enrichment_analysis/{model}.per_sample_enrichment_analysis.log",
    script:
        "../scripts/enrichment.R"

use rule per_sample_enrichment_analysis as per_celltype_enrichement_analysis with:
    input:
        diff_exp_genes="results/tables/seurat/integration/per_celltype/{model}.{celltype}.celltype-diff-exp-genes.tsv",
    output:
        GO_CC=report("results/tables/seurat/integration/enrichment/celltype/{model}.{celltype}.top_cellular_components.tsv",
                caption="../report/sample_enrichment_table_cc.rst", 
                category="Per model celltype enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "celltype":"{celltype}","table": "Significantly enriched Cellular components",
                },        
        ),
        GO_BP=report("results/tables/seurat/integration/enrichment/celltype/{model}.{celltype}.top_biological_process.tsv",
                caption="../report/sample_enrichment_table_bp.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "celltype":"{celltype}","table": "Significantly enriched Biological process",
                },          
        ),
        GO_MF=report("results/tables/seurat/integration/enrichment/celltype/{model}.{celltype}.top_molecular_function.tsv",
                caption="../report/sample_enrichment_table_mf.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "celltype":"{celltype}","table": "Significantly enriched Molecular functions",
                },         
        ),
        pathway=report("results/tables/seurat/integration/pathway/celltype/{model}.{celltype}.pathway_results.tsv",
                caption="../report/sample_pathway_table.rst",
                category="Per model celltype pathway analysis",
                subcategory="tables",
                labels={
                   "model": "{model}", "celltype":"{celltype}","table": "Significantly enriched pathways",
                },        
        
        ),
        
        bar_plot_cc=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.bar_plot_cc.pdf", 
                category="Per model celltype enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Cellular component bar plot",
                },         
        ),
        bar_plot_bp=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.bar_plot_bp.pdf", 
                category="Per model celltype enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Biological process bar plot",
                },         
        ),
        bar_plot_mf=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.bar_plot_mf.pdf", 
                category="Per model celltype enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Molecular function bar plot",
                },          
        ),
        dot_plot_cc=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.dot_plot_cc.pdf", 
                category="Per model celltype enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Cellular component dot plot",
                },         
        ),
        dot_plot_bp=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.dot_plot_bp.pdf", 
                category="Per model celltype enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Biological process dot plot",
                },        
        ),
        dot_plot_mf=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.dot_plot_mf.pdf", 
                category="Per model celltype enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Molecular function dot plot",
                },        
        ),
        tree_plot_cc=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.tree_plot_cc.pdf",
                caption="../report/sample_tree_plot.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Cellular component tree plot",
                },
        
        ),
        tree_plot_bp=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.tree_plot_bp.pdf",
                caption="../report/sample_tree_plot.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Biological process tree plot",
                }, 
        ),
        tree_plot_mf=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.tree_plot_mf.pdf",
                caption="../report/sample_tree_plot.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Molecular function tree plot",
                }, 
        ),
        enrichment_map_cc=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.enrichment_map_cc.pdf",
                
                category="Per model celltype enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Cellular component enrichment map",
                },        
        
        ),
        enrichment_map_bp=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.enrichment_map_bp.pdf",
                 caption="../report/sample_enrichment_map.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Biological process enrichment map",
                },
        
        
        ),
        enrichment_map_mf=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.enrichment_map_mf.pdf",
                caption="../report/sample_enrichment_map.rst",
                category="Per model celltype enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Molecular function enrichment map",
                }, 
        
        ),
        upset_plot_cc=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.upset_plot_cc.pdf",
                category="Per model celltype enrichment analysis",
                subcategory="Cellular component",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Cellular component upset plot",
                },         
        
        ),
        upset_plot_bp=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.upset_plot_bp.pdf",
                category="Per model celltype enrichment analysis",
                subcategory="Biological process",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Biological process upset plot",
                },

        ),
        upset_plot_mf=report("results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.upset_plot_mf.pdf",
                category="Per model celltype enrichment analysis",
                subcategory="Molecular function",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "Molecular function upset plot",
                },         
        
        ),
        dotplot_pathway=report("results/plots/seurat/integration/pathway/celltype/{model}.{celltype}.dotplot_pathway.pdf",
                category="Per model celltype pathway analysis",
                subcategory="dotplot",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "dot plot",
                },
        ),
        bar_plot_pathway=report("results/plots/seurat/integration/pathway/celltype/{model}.{celltype}.bar_plot_pathway.pdf",
                caption="../report/sample_pathway_barplot.rst",
                category="Per model celltype pathway analysis",
                subcategory="barplot",
                labels={
                   "model": "{model}", "celltype":"{celltype}","plot": "bar plot",
                },
        
        )
    log:
        "logs/seurat/integration/per_celltype_enrichement_analysis/{model}.{celltype}.per_celltype_enrichement_analysis.log",