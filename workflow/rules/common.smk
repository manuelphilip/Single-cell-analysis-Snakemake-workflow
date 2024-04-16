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
    sample = "|".join(samples.index),
    model="|".join(list(config["diffexp"].get("models", [])) + ["all"])
# model="|".join(list(config["diffexp"].get("models", [])) + ["all"]),

####### helpers ###########
"""

def check_config():
    representative_transcripts_keywords = ["canonical", "mostsignificant"]
    representative_transcripts = config["resources"]["ref"][
        "representative_transcripts"
    ]
    if representative_transcripts not in representative_transcripts_keywords:
        if not os.path.exists(representative_transcripts):
            raise ValueError(
                f"Invalid value given for resources/ref/representative_transcripts in "
                "configuration. Must be 'canonical', 'mostsignificant' or valid path, "
                "but {representative_transcripts} does not exist or is not readable."
            )


check_config()

"""


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

def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """

    wanted_input = []

    wanted_input.extend(
        expand(
            ["results/seurat/{model}.seurat_objt.rds",
            "results/plots/preprocessing/{model}.QC-Vln-plot.pdf",
            "results/plots/preprocessing/{model}.Highly_variable_features-plot.pdf",
            "results/plots/preprocessing/{model}.Elbow-plot.pdf"
        ],
        model=config["diffexp"]["models"],
        )
    )
    if config["clustering"]["activate"]:
        wanted_input.extend(
            expand(
                "results/plots/clustering/{model}.Dim-plot.pdf",
                model=config["diffexp"]["models"],    
            )
        
        )
    wanted_input.extend(
        expand(
            [#"results/tables/diffexp/{model}.diff-exp-genes.tsv",
            #"results/tables/diffexp/{model}.top-10-markers.tsv",
            #"results/plots/diffexp/{model}.Features-plot.pdf",
            "results/plots/diffexp/{model}.Heatmap-plot.pdf",
            #"results/plots/diffexp/{model}.Top-features-Vln-plot.pdf",
            "results/seurat/diffexp/{model}.seurat_objt.rds",

        ],
        model=config["diffexp"]["models"],
        )
    )
    return wanted_input
