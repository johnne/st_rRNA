import os
from snakemake.utils import validate
include: "scripts/common.py"
localrules: multiqc, report
rRNA = {"16S": ["silva-arc-16s-id95.fasta", "silva-bac-16s-id90.fasta"],
        "18S": ["silva-euk-18s-id95.fasta"]}
#TODO: Figure out how to link taxonomy/counts to spots
container: "docker://continuumio/miniconda3:4.9.2"

validate(config,schema="config/config.schema.yml",set_default=True)
samples = parse_samples(config["sample_list"])

rule report:
    input:
        expand("results/taxonomy/{sample}.{subunit}.taxonomy.tsv",
            sample=samples.keys(),
            subunit = config["subunits"]),
        "results/report/sample_report.html"
    log: "results/log/report.log"
    shell:
        """
        exec &> {log}
        mkdir -p results/report
        snakemake --unlock 
        snakemake --report results/report/snakemake_report.html
        """

rule download_metaxa_db:
    output:
        expand("resources/metaxa2/blast.{s}",
            s = ["nhr","nin","nsd","nsi","nsq","taxonomy.txt","cutoffs.txt"]),
        expand("resources/metaxa2/HMMs/{L}.hmm{suffix}",
            L = ["A","B","C","E","M","N"],
            suffix = ["",".h3f",".h3i",".h3m",".h3p"])
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    conda: "envs/metaxa.yml"
    shell:
        """
        metaxa2_install_database -g SSU -d {params.outdir}
        """

rule metaxa:
    conda: "envs/metaxa.yml"
    shell:
        """
        metaxa2_x 
        """

rule download_sortmerna_db:
    """
    Downloads fasta files for 16S/18S for use with sortmerna
    """
    output:
        "resources/sortmerna/{f}"
    params:
        url_base = "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/"
    log:
        "resources/sortmerna/download-{f}.log"
    shell:
        """
        curl -L -o {output[0]} {params.url_base}/{wildcards.f} > {log} 2>&1
        """

rule sortmerna:
    """
    Runs sortmerna for each sample/subunit combination. 
    
    To run with both 16S and 18S subunits, both should be specified in config.
    """
    input:
        R2 = lambda wildcards: samples[wildcards.sample]["R2"],
        db = rRNA_fasta
    output:
        "results/rRNA/{subunit}/{sample}.rRNA.fastq.gz"
    conda: "envs/sortmerna.yml"
    log:
        runlog="results/logs/rRNA/{subunit}/{sample}.log",
        reportlog="results/rRNA/{subunit}/{sample}.log"
    params:
        workdir = "$TMPDIR/rRNA/{subunit}.{sample}.wd",
        R2 = "$TMPDIR/rRNA/{subunit}.{sample}.wd/R2.fastq",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        ref_string = lambda wildcards, input: " ".join([f"--ref {x}" for x in input.db])
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 120
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        gunzip -c {input.R2} > {params.R2}
        sortmerna --threads {threads} --workdir {params.workdir} --fastx \
            --reads {params.R2} {params.ref_string}
            --aligned {params.workdir}/{wildcards.sample}.rRNA > {log.runlog} 2>&1
        rm {params.R2}
        gzip {params.workdir}/*.fastq
        mv {params.workdir}/*.gz {params.outdir}
        mv {params.workdir}/{wildcards.sample}.rRNA.log {log.reportlog}
        #rm -rf {params.workdir}
        """

rule multiqc:
    """
    Runs multiqc on sortmerna logs to summarize counts
    """
    input:
        expand("results/rRNA/{subunit}/{sample}.log",
            sample = samples.keys(), subunit = config["subunits"])
    output:
        report("results/report/sample_report.html")
    log:
        "results/logs/multiqc/multiqc.log"
    params:
        file_list = "results/report/multiqc_file_list.txt",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        name = lambda wildcards, output: os.path.basename(output[0])
    conda: "envs/multiqc.yml"
    shell:
        """
        echo {input} | tr " " "\n" > {params.file_list}
        multiqc -f --file-list {params.file_list} -o {params.outdir} -n {params.name} > {log} 2>&1
        rm {params.file_list}
        """

rule download_train_set:
    """
    Downloads a DADA2-ready train set for use with the assignTaxonomy function
    """
    output:
        "resources/DADA2/{subunit}.train.fasta"
    log:
        "results/logs/DADA2/download.{subunit}.train.log"
    params:
        url=lambda wildcards: config["assignTaxonomy"]["database"][wildcards.subunit]
    shell:
        """
        curl -L -o {output[0]}.gz {params.url} > {log} 2>&1
        gunzip {output[0]}.gz
        """

rule assign_taxonomy:
    """
    Assigns taxonomy to the R2 reads identified as rRNA by sortmerna 
    """
    input:
        seqs = "results/rRNA/{subunit}/{sample}.rRNA.fastq.gz",
        refFasta = "resources/DADA2/{subunit}.train.fasta"
    output:
        taxdf = report("results/taxonomy/{sample}.{subunit}.taxonomy.tsv"),
        bootdf = "results/taxonomy/{sample}.{subunit}.bootstrap.tsv"
    log:
        "results/logs/assignTaxonomy.{sample}.{subunit}.log"
    params:
        minBoot = config["assignTaxonomy"]["minBoot"],
        outputBootstraps = config["assignTaxonomy"]["outputBootstraps"]
    conda: "envs/dada2.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    threads: 10
    script:
        "scripts/assignTaxa.R"