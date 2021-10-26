import os
from snakemake.utils import validate
include: "scripts/common.py"
localrules:
    multiqc,
    report,
    download_sortmerna_db,
    download_train_set,
    extract_R1,
    sample,
    gather_spots

configfile: "config/config.yml"
rRNA = {"16S": ["silva-arc-16s-id95.fasta", "silva-bac-16s-id90.fasta"],
        "18S": ["silva-euk-18s-id95.fasta"]}
container: "docker://continuumio/miniconda3:4.9.2"

validate(config,schema="config/config.schema.yml",set_default=True)
samples = parse_samples(config["sample_list"])
barcodes = read_barcodes(config["barcode_list"])

wildcard_constraints:
    barcode = "({})".format("|".join(barcodes.keys())),

def taxinput(config, samples):
    input = []
    if config["taxtool"] == "vsearch":
        input += expand("results/taxonomy/{sample}.{subunit}.vsearch.tsv",
                        sample = samples.keys(), subunit = config["subunits"])
    else:
        input += expand("results/taxonomy/{sample}.{subunit}.assignTaxonomy.tsv",
                        sample = samples.keys(), subunit = config["subunits"])
    return input

rule report:
    input:
        taxinput(config, samples),
        "results/report/sample_report.html"

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
        "results/rRNA/{subunit}/{sample}.R2.rRNA.fastq.gz"
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
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        gunzip -c {input.R2} > {params.R2}
        sortmerna --threads {threads} --workdir {params.workdir} --fastx \
            --reads {params.R2} {params.ref_string} \
            --aligned {params.workdir}/{wildcards.sample}.R2.rRNA > {log.runlog} 2>&1
        rm {params.R2}
        gzip {params.workdir}/*.fastq
        mv {params.workdir}/*.gz {params.outdir}
        mv {params.workdir}/{wildcards.sample}.R2.rRNA.log {log.reportlog}
        rm -rf {params.workdir}
        """

rule extract_R1:
    """
    Extracts R1 sequences for reads identified as rRNA
    """
    input:
        R1 = lambda wildcards: samples[wildcards.sample]["R1"],
        R2 = "results/rRNA/{subunit}/{sample}.R2.rRNA.fastq.gz"
    output:
        R1 = "results/rRNA/{subunit}/{sample}.R1.rRNA.fastq.gz"
    log:
        "results/logs/rRNA/{subunit}.{sample}.seqtk.log"
    params:
        ids = "$TMPDIR/{subunit}.{sample}.ids",
        R1 = "$TMPDIR/{subunit}.{sample}.R1.rRNA.fastq.gz"
    conda: "envs/seqkit.yml"
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        gunzip -c {input.R2} | egrep "^@" | cut -f1 -d ' ' | sed 's/@//g' > {params.ids}
        seqkit grep -f {params.ids} {input.R1} | gzip -c > {params.R1}
        mv {params.R1} {output.R1}
        """

rule sample_spots:
    """
    1. Splits fastq file into spots
    2. Clusters R1 sequences by UMIs
    3. Samples X reads per UMI cluster
    """
    input:
        R1 = "results/rRNA/{subunit}/{sample}.R1.rRNA.fastq.gz",
        R2 = "results/rRNA/{subunit}/{sample}.R2.rRNA.fastq.gz"
    output:
        R2 = temp("results/rRNA/{subunit}/spots/{barcode}/{sample}.R2.rRNA.subsampled.fastq.gz"),
        f = temp("results/rRNA/{subunit}/spots/{barcode}/{sample}.map.tsv")
    log:
        "results/logs/{subunit}/{sample}.sample_spots.{barcode}.log"
    params:
        X = config["seqs_per_spot"],
        umi_start = config["UMI"]["start"],
        umi_stop = config["UMI"]["stop"],
        ids = "$TMPDIR/{subunit}.{sample}.{barcode}/{sample}",
        tmpdir = "$TMPDIR/{subunit}.{sample}.{barcode}",
        R1 = "$TMPDIR/{subunit}.{sample}.{barcode}/{sample}.{barcode}.R1.rRNA.fastq",
        R2 = "$TMPDIR/{subunit}.{sample}.{barcode}/{sample}.{barcode}.R2.rRNA.subsampled.fastq.gz"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 2
    conda: "envs/seqkit.yml"
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        mkdir -p {params.tmpdir}
        seqkit grep -m 0 -P -r -s -p "^{wildcards.barcode}" {input.R1} > {params.R1}
        python scripts/sample_spot.py {params.R1} --num_reads {params.X} \
            --start {params.umi_start} --stop {params.umi_stop} \
            --barcode {wildcards.barcode} --mapfile {output.f} > {params.ids}
        seqkit grep -f {params.ids} {input.R2} | gzip -c > {params.R2}
        mv {params.R2} {output.R2}
        """

rule gather_spots:
    input:
        fq = expand("results/rRNA/{{subunit}}/spots/{barcode}/{{sample}}.R2.rRNA.subsampled.fastq.gz",
            barcode = barcodes.keys()),
        tsv = expand("results/rRNA/{{subunit}}/spots/{barcode}/{{sample}}.map.tsv",
            barcode = barcodes.keys())
    output:
        "results/rRNA/{subunit}/{sample}.sampled.R2.rRNA.fastq.gz",
        "results/rRNA/{subunit}/{sample}.map.tsv"
    params:
        tmp1 = "$TMPDIR/{subunit}.{sample}.sampled.R2.rRNA.fastq.gz",
        tmp2 = "$TMPDIR/{subunit}.{sample}.mapfile.tsv"
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        cat {input.fq} > {params.tmp1}
        mv {params.tmp1} {output[0]}
        cat {input.tsv} > {params.tmp2}
        mv {params.tmp2} {output[1]}
        """

## Sample target rule
rule sample:
    input:
        expand("results/rRNA/{subunit}/{sample}.sampled.R2.rRNA.fastq.gz",
            subunit = config["subunits"], sample = samples.keys()),
        expand("results/rRNA/{subunit}/{sample}.map.tsv",
            subunit = config["subunits"], sample = samples.keys())

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
        url=lambda wildcards: config["assignTaxonomy"]["taxdb"][wildcards.subunit]
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
        seqs = "results/rRNA/{subunit}/{sample}.sampled.R2.rRNA.fastq.gz",
        refFasta = "resources/DADA2/{subunit}.train.fasta"
    output:
        taxdf = report("results/taxonomy/{sample}.{subunit}.assignTaxonomy.tsv"),
        bootdf = "results/taxonomy/{sample}.{subunit}.assignTaxonomy.bootstrap.tsv"
    log:
        "results/logs/assignTaxonomy.{sample}.{subunit}.log"
    params:
        minBoot = config["assignTaxonomy"]["minBoot"],
        outputBootstraps = config["assignTaxonomy"]["outputBootstraps"],
        tryRC = config["assignTaxonomy"]["tryRC"]
    conda: "envs/dada2.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    threads: 20
    script:
        "scripts/assignTaxa.R"

rule download_vsearch:
    output:
        "resources/sintax/{subunit}.fasta"
    log:
        "results/logs/vsearch/{subunit}.download.log"
    params:
        url = lambda wildcards: config["vsearch"][wildcards.subunit]
    shell:
        """
        curl -L -o {output[0]}.gz {params.url} > {log} 2>&1
        gunzip {output[0]}.gz
        """

rule vsearch:
    input:
        seqs = "results/rRNA/{subunit}/{sample}.sampled.R2.rRNA.fastq.gz",
        db = "resources/sintax/{subunit}.fasta"
    output:
        report("results/taxonomy/{sample}.{subunit}.vsearch.tsv")
    log:
        "results/logs/vsearch/{sample}.{subunit}.log"
    params:
        cutoff = config["vsearch"]["cutoff"]
    conda: "envs/vsearch.yml"
    threads: 10
    shell:
        """
        vsearch --threads {threads} --sintax {input.seqs} --db {input.db} \
            --sintax_cutoff {params.cutoff} --tabbedout {output[0]} > {log} 2>&1 
        """