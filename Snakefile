import pandas as pd
import os
from snakemake.utils import validate
#TODO: Also identify 16S rRNA
rRNA = {"16S": ["silva-arc-16s-id95.fasta", "silva-bac-16s-id90.fasta"],
        "18S": ["silva-euk-18s-id95.fasta"]}
#TODO: Figure out how to link taxonomy/counts to spots
container: "docker://continuumio/miniconda3:4.9.2"

validate(config,schema="config/config.schema.yml",set_default=True)


def parse_samples(f):
    df = pd.read_csv(f, index_col=0, sep="\t", header=0)
    return df.to_dict(orient="index")

samples = parse_samples(config["sample_list"])

rule all:
    input:
        expand("results/taxonomy/{sample}.taxonomy.tsv", sample=samples.keys())

rule download_sortmerna_db:
    output:
        "resources/sortmerna/{f}"
    params:
        url_base = "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/"
    log:
        "resources/sortmerna/download-{f}.log"
    shell:
        """
        curl -L -o $f {params.url_base}/{wildcards.f} > {log} 2>&1
        """

def rRNA_fasta(wildcards):
    input = []
    for f in rRNA[wildcards.subunit]:
        input.append(f"resources/sortmerna/{f}")
    return input

rule sortmerna:
    input:
        R1 = lambda wildcards: samples[wildcards.sample]["R1"],
        R2 = lambda wildcards: samples[wildcards.sample]["R2"],
        db = rRNA_fasta
    output:
        R1 = "results/rRNA/{subunit}/{sample}.rRNA_fwd.fastq.gz",
        R2 = "results/rRNA/{subunit}/{sample}.rRNA_rev.fastq.gz"
    conda: "envs/sortmerna.yml"
    log:
        runlog="results/logs/rRNA/{subunit}/{sample}.log",
        reportlog="results/rRNA/{subunit}/{sample}.log"
    params:
        workdir = "$TMPDIR/rRNA/{subunit}.{sample}.wd",
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

rule download_pr2:
    output:
        "resources/pr2/pr2.dada2.fasta"
    log:
        "results/logs/PR2/download.log"
    params:
        url=config["PR2"]["url"]
    shell:
        """
        curl -L -o {output[0]}.gz {params.url} > {log} 2>&1
        gunzip {output[0]}.gz
        """

rule assign_taxonomy:
    input:
        seqs = "results/rRNA/{sample}.rRNA_rev.fastq.gz",
        refFasta = "resources/pr2/pr2.dada2.fasta"
    output:
        taxdf = "results/taxonomy/{sample}.taxonomy.tsv",
        bootdf = "results/taxonomy/{sample}.bootstrap.tsv"
    log:
        "results/logs/assignTaxonomy.{sample}.log"
    params:
        minBoot = config["assignTaxonomy"]["minBoot"],
        outputBootstraps = config["assignTaxonomy"]["outputBootstraps"]
    conda: "envs/dada2.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    threads: 10
    script:
        "scripts/assignTaxa.R"
