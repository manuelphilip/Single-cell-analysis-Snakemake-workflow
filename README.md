# Snakemake workflow: SingleCell-10X-Genomics-RNA-seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/manuelphilip/Single-cell-analysis-Snakemake-workflow/workflows/Tests/badge.svg?branch=main)](https://github.com/manuelphilip/Single-cell-analysis-Snakemake-workflow/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for single-cell data analysis of 10X genomics data including cell type annotation, differential expression (marker gene identification), scRNA-seq integration


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<manuelphilip>%2F<Single-cell-analysis-Snakemake-workflow>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <Single-cell-analysis-Snakemake-workflow>sitory and its DOI (see above).

# Description
 A Snakemake workflow designed to analyse single cell data from cell ranger output. The workflow expects the cell ranger output which contains per sample bc_matrix (raw and filtered) under the `~/sample_name/outs` folder.

