rule download_transcriptome:
    output:
        "resources/refdata-gex-{build}-{release}-A.tar.gz",
    log:
        "logs/download_transcriptome/{build}-{release}.log",
    params:
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
    shell:
        "wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-{params.build}-{params.release}-A.tar.gz -O {output} &> {log}"


rule extract_transcriptome_data:
    input: 
        transcriptome_zipped=get_transcriptome,
    output:
        transcriptome=directory("resources/refdata/refdata-gex-{build}-{release}-A"),
    params:
        # Somehow only one directory before the last directory is created in the ' snakemake output',
        outdir=lambda w, output: output[0][:-26],
    log:
        "logs/extract_transcriptome_data/refdata-gex-{build}-{release}-A.log",
    shell:
        "mkdir -p {params.outdir} &&"
        "tar -zxvf {input.transcriptome_zipped} -C {params.outdir}  &> {log}"
