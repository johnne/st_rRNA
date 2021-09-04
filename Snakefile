import pandas as pd
import os
from snakemake.utils import validate

container: "docker://continuumio/miniconda3:4.9.2"

validate(config,schema="config/config.schema.yml",set_default=True)


def parse_samples(f):
    df = pd.read_csv(f, index_col=0, sep="\t", header=0)
    return df.to_dict(orient="index")

samples = parse_samples(config["sample_list"])

rule all:
    input:
        "results/dada2/taxonomy.tsv"

rule download_sortmerna_db:
    output:
        "resources/sortmerna/silva-euk-18s-id95.fasta"
    params:
        url_base = "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/"
    log:
        "resources/sortmerna/silva-euk-18s-id95.log"
    shell:
        """
        curl -L -o {output[0]} {params.url_base}/silva-euk-18s-id95.fasta > {log} 2>&1
        """

rule sortmerna:
    input:
        R1 = lambda wildcards: samples[wildcards.sample]["R1"],
        R2 = lambda wildcards: samples[wildcards.sample]["R2"],
        db = "resources/sortmerna/silva-euk-18s-id95.fasta"
    output:
        R1 = "results/rRNA/{sample}_R1.rRNA.fastq.gz",
        R2 = "results/rRNA/{sample}_R2.rRNA.fastq.gz"
    conda: "envs/sortmerna.yml"
    log:
        runlog="results/logs/sortmerna/{sample}.log",
        reportlog="results/sortmerna/{sample}.log"
    params:
        workdir = "$TMPDIR/sortmerna/{sample}.wd",
        outdir = lambda wildcards, output: os.path.dirname(output.R1)
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 48
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        sortmerna --threads {threads} --workdir {params.workdir} --fastx \
            --reads {input.R1} --reads {input.R2} --ref {input.db} --paired_in \
            --out2 --aligned {params.workdir}/{wildcards.sample}.rRNA > {log.runlog} 2>&1
        gzip {params.workdir}/*.fastq
        mv {params.workdir}/*.gz {params.outdir}
        mv {params.workdir}/{wildcards.sample}.rRNA.log {log.reportlog}
        rm -rf {params.workdir}
        """