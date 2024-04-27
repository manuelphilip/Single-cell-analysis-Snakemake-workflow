from snakemake.utils import validate
import pandas as pd
import yaml
from pathlib import Path

##### load config and sample sheets #####


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index(
    "sample", drop=False
)
samples.index.names = ["sample_id"]


def drop_unique_cols(df):
    singular_cols = df.nunique().loc[(df.nunique().values <= 1)].index
    return df.drop(singular_cols, axis=1)


samples = drop_unique_cols(samples)
validate(samples, schema="../schemas/samples.schema.yaml")

#units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(
#    ["sample", "path"], drop=False
#)
#validate(units, schema="../schemas/units.schema.yaml")
# print(units)


report: "../report/workflow.rst"


##### wildcard constraints #####

wildcard_constraints:
    sample="|".join(samples.index),
    models="|".join(list(config["diffexp"].get("models", [])) + ["all"]),

####### helpers ###########

#def get_fastqs(wildcards):
#    """Determine whether unit is single-end."""
#    fastq_path = units.loc[(sample), "path"]

#    return fastq_path


def get_transcriptome(wildcards):
    transcriptome = expand(
        ["resources/refdata-gex-{build}-{release}-A.tar.gz"],
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
    )
    return transcriptome


def get_model_samples(wildcards):
    samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    gps = config["diffexp"]["models"][wildcards.model]["column_name"]
    sample_groups = samples.loc[samples[gps].notnull(), ["sample"]]
    condition_groups = samples.loc[samples[gps].notnull(), [gps]]
    condition = condition_groups[gps].values
    samples = sample_groups["sample"].values
    return samples, condition

def get_model(wildcards):
    if wildcards.model == "all":
        return {"full": None}
    return config["diffexp"]["models"][wildcards.model]

def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """

    wanted_input = []

    wanted_input.extend(
        expand(
            ["results/seurat/preprocessing/all.seurat_objt.rds",
            "results/plots/preprocessing/all.QC-Vln-plot.pdf",
            "results/plots/preprocessing/all.Highly_variable_features-plot.pdf",
            "results/plots/preprocessing/all.Elbow-plot.pdf"
        ],
        model=config["diffexp"]["models"],
        )
    )
    if config["clustering"]["activate"]:
        wanted_input.extend(
            expand(
            "results/plots/clustering/all.Dim-plot.pdf"),
        )

    if config["visualize_marker_expression"]["activate"]:
        wanted_input.extend(
            expand(
                [
                "results/tables/diffexp/{unit.sample}.diff-exp-genes.tsv",
                "results/tables/diffexp/{unit.sample}.top-10-markers.tsv",
                "results/plots/diffexp/{unit.sample}.Heatmap-plot.pdf",    
                "results/plots/diffexp/{unit.sample}.Top-features-Vln-plot.pdf",
                "results/plots/diffexp/{unit.sample}.Features-plot.pdf",
                "results/plots/celltype/{unit.sample}.Dim-plot.pdf"
                ],
                unit=samples[["sample"]].itertuples(),

            )
        
        )
    if config["celltype_annotation"]["activate"]:
        wanted_input.extend(
            expand(
                [
                "results/plots/celltype/{unit.sample}.Dim-plot.pdf"
                ],
                unit=samples[["sample"]].itertuples(),

            )
        
        )
    if config["diffexp"]["activate"]:
        wanted_input.extend(
            expand(
                [
                "results/seurat/integration/merge/{model}.seurat_objt_before_integration.rds",
                ],
                model=config["diffexp"]["models"],
            )
        
        )
    if config["merged_sample_clustering"]["activate"]:
        wanted_input.extend(
            expand(
                [
                "results/plots/seurat/integration/merge/clustering/{model}.Dim-plot.pdf",
                #Need to fix the rule issue plot_integration_plots
                #"results/plots/seurat/integration/integration/{model}.Dim-plot.pdf",
                "results/tables/seurat/integration/{model}.diff-exp-genes.tsv",
                "results/plots/seurat/integration/per_celltype/{model}.Dim-plot.pdf",
                "results/tables/seurat/integration/{model}.top-10-markers.tsv",
                "results/plots/seurat/integration/{model}.Heatmap-plot.pdf",
                "results/plots/seurat/integration/volcano_plots/{model}.volcano_plot.pdf",
                ],
                model=config["diffexp"]["models"],
            )
        
        )
    if config["cell_type_diff_exp"]["activate"]:
        wanted_input.extend(
            expand(
                [
                "results/tables/seurat/integration/per_celltype/{model}.{celltype}.celltype-diff-exp-genes.tsv",
                "results/tables/seurat/integration/per_celltype/{model}.{celltype}.top-10-celltype-markers.tsv",
                "results/plots/seurat/integration/per_celltype/{model}.{celltype}.celltype-Heatmap-plot.pdf",
                "results/plots/seurat/integration/per_celltype/volcano_plots/{model}.{celltype}.celltype-volcano_plot.pdf",
                ],
                model=config["diffexp"]["models"],
                celltype=config["cell_type_diff_exp"]["cell_type"],
            )
        
        )
        if config["enrichment"]["activate"]:
            wanted_input.extend(
                expand(
                    [
                    "results/tables/seurat/integration/enrichment/celltype/{model}.{celltype}.top_cellular_components.tsv",
                    "results/tables/seurat/integration/enrichment/celltype/{model}.{celltype}.top_biological_process.tsv",
                    "results/tables/seurat/integration/enrichment/celltype/{model}.{celltype}.top_molecular_function.tsv",
                    "results/tables/seurat/integration/pathway/celltype/{model}.{celltype}.pathway_results.tsv",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.bar_plot_cc.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.bar_plot_bp.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.bar_plot_mf.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.dot_plot_cc.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.dot_plot_bp.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.dot_plot_mf.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.tree_plot_cc.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.tree_plot_bp.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.tree_plot_mf.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.enrichment_map_cc.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.enrichment_map_bp.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.enrichment_map_mf.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.upset_plot_cc.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.upset_plot_bp.pdf",
                    "results/plots/seurat/integration/enrichment/celltype/{model}.{celltype}.upset_plot_mf.pdf",
                    "results/plots/seurat/integration/pathway/celltype/{model}.{celltype}.dotplot_pathway.pdf",
                    "results/plots/seurat/integration/pathway/celltype/{model}.{celltype}.bar_plot_pathway.pdf",
                ],
                model=config["diffexp"]["models"],
                celltype=config["cell_type_diff_exp"]["cell_type"],
            )
        )
    if config["enrichment"]["activate"]:
        wanted_input.extend(
            expand(
                [
                "results/tables/seurat/integration/enrichment/groups/{model}.top_cellular_components.tsv",
                "results/tables/seurat/integration/enrichment/groups/{model}.top_biological_process.tsv",
                "results/tables/seurat/integration/enrichment/groups/{model}.top_molecular_function.tsv",
                "results/tables/seurat/integration/pathway/groups/{model}.pathway_results.tsv",
                "results/plots/seurat/integration/enrichment/groups/{model}.bar_plot_cc.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.bar_plot_bp.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.bar_plot_mf.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.dot_plot_cc.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.dot_plot_bp.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.dot_plot_mf.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.tree_plot_cc.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.tree_plot_bp.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.tree_plot_mf.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.enrichment_map_cc.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.enrichment_map_bp.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.enrichment_map_mf.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.upset_plot_cc.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.upset_plot_bp.pdf",
                "results/plots/seurat/integration/enrichment/groups/{model}.upset_plot_mf.pdf",
                "results/plots/seurat/integration/pathway/groups/{model}.bar_plot_pathway.pdf",
                "results/plots/seurat/integration/pathway/groups/{model}.dotplot_pathway.pdf",
                ],
                model=config["diffexp"]["models"],
            )
        
        )
    return wanted_input
