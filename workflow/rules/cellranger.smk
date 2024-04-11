rule cellranger_count:
    input:
        fastq=get_fastqs,
        transcriptome=expand(
            ["resources/refdata/refdata-gex-{build}-{release}-A/"],
            build=config["resources"]["ref"]["build"],
            release=config["resources"]["ref"]["release"],
        ),
    output:
        cellranger_count=directory("results/cellranger_count_out/{samples}/"),
    log:
        "logs/cellranger_count/{samples}.log",
    params:
        samples=expand(
            ["{unit.sample}"],
            unit=units[["sample"]].itertuples(),
        ),
    threads: 55
    shell:
        # more info https://github.com/10XGenomics/cellranger/issues/82
        "mkdir -p {output.cellranger_count} &&"
        "cellranger count --id={params.samples} --fastqs {input.fastq} --sample={params.samples} --localcores={threads} --transcriptome {input.transcriptome} &>{log} && mv {params.samples}/outs/*  {output.cellranger_count} && rm -r {params.samples}"
