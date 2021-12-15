import os
from snakemake.utils import validate

include: "scripts/common.py"
localrules:
    multiqc,
    report,
    extract_R1,
    sample,
    gather_spots,
    sample_spots,
    concat_fasta

configfile: "config/config.yml"
rRNA = {"16S": ["silva-arc-16s-id95.fasta", "silva-bac-16s-id90.fasta"],
        "18S": ["silva-euk-18s-id95.fasta"]}
container: "docker://continuumio/miniconda3:4.9.2"

validate(config,schema="config/config.schema.yml",set_default=True)
samples = parse_samples(config["sample_list"])
barcodes = read_barcodes(config["barcode_list"])

for sample, d in samples.items():
    try:
        d["filter_rank"]
    except KeyError:
        samples[sample]["filter_rank"] = "genus"

wildcard_constraints:
    barcode = "({})".format("|".join(barcodes.keys())),
    subunit = "16S|18S",
    sample = "({})".format("|".join(samples.keys())),

rule report:
    input:
        expand("results/metaxa2/{sample}.taxonomy.txt",
            sample = samples.keys()),
        expand("results/metaxa2/{sample}.{domain}.fasta",
            sample = samples.keys(), domain = ["archaea","bacteria","chloroplast","eukaryota"]),
        "results/report/multiqc.html"

rule cutadapt:
    input:
        R1 = lambda wildcards: samples[wildcards.sample]["R1"],
        R2 = lambda wildcards: samples[wildcards.sample]["R2"]
    output:
        R1 = "results/cutadapt/{sample}.cutadapt.R1.fastq.gz",
        R2 = "results/cutadapt/{sample}.cutadapt.R2.fastq.gz"
    conda:
        "envs/cutadapt.yml"
    log:
        "results/logs/cutadapt/{sample}.log"
    threads: 4
    params:
        R1_adapter = config["cutadapt"]["R1_adapter"],
        R2_adapter=config["cutadapt"]["R2_adapter"],
        minlen=config["cutadapt"]["minlen"],
        extra_params=config["cutadapt"]["extra_params"]
    resources:
        runtime = lambda wildcards,attempt: attempt ** 2 * 60 * 2
    shell:
        """
        cutadapt {params.extra_params} -m {params.minlen} -j {threads} \
            -a {params.R1_adapter} -A {params.R2_adapter} -o {output.R1} \
            -p {output.R2} {input.R1} {input.R2} > {log} 2>&1
        """

rule fastp:
    input:
        R2 = rules.cutadapt.output.R2
    output:
        R2 = "results/fastp/{sample}_R2.fastq.gz"
    log:
        "results/logs/fastp/{sample}.log"
    conda: "envs/fastp.yml"
    params:
        complexity_cutoff = config["fastp"]["complexity_cutoff"],
        extra_params = config["fastp"]["extra_params"]
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 2
    shell:
        """
        fastp {params.extra_params} -y -Y {params.complexity_cutoff} \
            -i {input.R2} -o {output.R2} --thread 4 > {log} 2>&1
        """

rule extract_R1:
    """
    Extracts R1 sequences
    """
    input:
        R1 = lambda wildcards: samples[wildcards.sample]["R1"],
        R2 = rules.fastp.output.R2
    output:
        R1 = "results/fastp/{sample}_R1.fastq.gz"
    log:
        "results/logs/seqtk/{sample}.log"
    params:
        ids = "$TMPDIR/{sample}.ids",
        R1 = "$TMPDIR/{sample}.R1.fastp.fastq.gz"
    conda: "envs/seqkit.yml"
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        gunzip -c {input.R2} | egrep "^@" | cut -f1 -d ' ' | sed 's/@//g' > {params.ids}
        seqkit grep -f {params.ids} {input.R1} | gzip -c > {params.R1}
        mv {params.R1} {output.R1}
        """

rule fastqc:
    input:
        R1 = rules.extract_R1.output.R1,
        R2 = rules.fastp.output.R2
    output:
        R1 = "results/fastqc/{sample}_R1_fastqc.zip",
        R2 = "results/fastqc/{sample}_R2_fastqc.zip"
    log:
        "results/logs/fastqc/{sample}.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.R1)
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 2
    conda:
        "envs/fastqc.yml"
    shell:
        """
        fastqc --noextract -o {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """

rule multiqc:
    input:
        expand("results/logs/fastp/{sample}.log",
            sample=samples.keys()),
        expand("results/logs/cutadapt/{sample}.log",
            sample=samples.keys()),
        expand("results/fastqc/{sample}_{R}_fastqc.zip",
            R=["R1","R2"], sample=samples.keys())
    output:
        "results/report/multiqc.html"
    log:
        "results/logs/multiqc.log"
    conda:
        "envs/multiqc.yml"
    params:
        filelist = lambda wildcards, output: os.path.dirname(output[0])+"/filelist",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        echo {input} | tr ' ' '\n' > {params.filelist}
        multiqc -o {params.outdir} -n multiqc.html -l {params.filelist} -f > {log} 2>&1
        """

rule sample_spots:
    """
    1. Splits fastq file into spots
    2. Clusters R1 sequences by UMIs
    3. Samples X reads per UMI cluster
    """
    input:
        R1 = rules.extract_R1.output.R1,
        R2 = rules.fastp.output.R2
    output:
        R2 = temp("results/spots/{barcode}/{sample}.R2.fastq.gz"),
        f = temp("results/spots/{barcode}/{sample}.map.tsv"),
        stats = temp("results/spots/{barcode}/{sample}.stats")
    log:
        "results/logs/spots/{sample}.{barcode}.log"
    params:
        X = config["seqs_per_spot"],
        umi_start = config["UMI"]["start"],
        umi_stop = config["UMI"]["stop"],
        ids = "$TMPDIR/{sample}.{barcode}/{sample}",
        tmpdir = "$TMPDIR/{sample}.{barcode}",
        R1 = "$TMPDIR/{sample}.{barcode}/{sample}.{barcode}.R1.rRNA.fastq",
        R2 = "$TMPDIR/{sample}.{barcode}/{sample}.{barcode}.R2.rRNA.subsampled.fastq.gz"
    conda: "envs/seqkit.yml"
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        mkdir -p {params.tmpdir}
        seqkit grep -m 0 -P -r -s -p "^{wildcards.barcode}" {input.R1} > {params.R1}
        python scripts/sample_spot.py {params.R1} --stats {output.stats} \
            --num_reads {params.X} --start {params.umi_start} --stop {params.umi_stop} \
            --barcode {wildcards.barcode} --mapfile {output.f} > {params.ids}
        seqkit grep -f {params.ids} {input.R2} | gzip -c > {params.R2}
        mv {params.R2} {output.R2}
        """

## Run Metaxa2
rule metaxa:
    input:
        R2 = rules.sample_spots.output.R2
    output:
        expand("results/metaxa2/spots/{{barcode}}/{{sample}}.{domain}.fasta",
            domain=["archaea","bacteria","chloroplast","eukaryota"]),
        expand("results/metaxa2/spots/{{barcode}}/{{sample}}.{f}.txt",
            f = ["taxonomy","summary"])
    log:
        "results/logs/metaxa2/{sample}.{barcode}.log"
    params:
        out = "results/metaxa2/spots/{barcode}/{sample}",
        qual_percent = 10,
        qual = 20,
        relscore = 0,
        R2 = "$TMPDIR/{sample}.{barcode}.R2.fastq"
    conda: "envs/metaxa.yml"
    threads: 2
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 2
    shell:
        """
        gunzip -c {input.R2} > {params.R2}
        metaxa2 -R 0 --cpu {threads} -i {params.R2} --plus -g ssu -f fastq \
            --quality_trim T -q {params.qual} \
            --quality_percent {params.qual_percent} -o {params.out} > {log} 2>&1
        rm {params.R2}
        """

def concat_files(files, outfile):
    import shutil
    with open(outfile,'wb') as fhout:
        for f in files:
            with open(f,'rb') as fhin:
                shutil.copyfileobj(fhin,fhout)

def concat_df(files, outfile):
    import pandas
    df = pandas.DataFrame()
    for f in files:
        try:
            _df = pandas.read_csv(f, sep="\t", index_col=0, header=0)
        except pandas.errors.EmptyDataError:
            continue
        df = pandas.concat([df, _df])
    df.to_csv(outfile, sep="\t")

rule gather_spots:
    """
    Concatenates spot files per subunit/sample
    """
    input:
        tax = expand("results/metaxa2/spots/{barcode}/{{sample}}.taxonomy.txt",
            barcode = barcodes.keys()),
        tsv = expand("results/spots/{barcode}/{{sample}}.map.tsv",
            barcode = barcodes.keys()),
        stats = expand("results/spots/{barcode}/{{sample}}.stats",
            barcode = barcodes.keys())
    output:
        tax = "results/metaxa2/{sample}.taxonomy.txt",
        tsv = "results/metaxa2/{sample}.map.tsv",
        stats = "results/metaxa2/{sample}.stats.tsv"
    run:
        # Concatenate tsv and stats files
        concat_files(input.tsv, output.tsv)
        concat_files(input.stats, output.stats)
        concat_files(input.tax, output.tax)

## Sample target rule
rule sample:
    input:
        tax=expand("results/metaxa2/spots/{barcode}/{sample}.taxonomy.txt",
            barcode=barcodes.keys(), sample=samples.keys()),
        tsv=expand("results/spots/{barcode}/{sample}.map.tsv",
            barcode=barcodes.keys(), sample=samples.keys()),
        stats=expand("results/spots/{barcode}/{sample}.stats",
            barcode=barcodes.keys(), sample=samples.keys())

rule concat_fasta:
    input:
        expand("results/metaxa2/spots/{barcode}/{{sample}}.{{domain}}.fasta",
            barcode = barcodes.keys())
    output:
        "results/metaxa2/{sample}.{domain}.fasta"
    run:
        concat_files(input, output[0])

rule spot_taxonomy:
    input:
        tax = rules.gather_spots.output.tax,
        mapfile = rules.gather_spots.output.tsv,
        barcodes = config["barcode_list"]
    output:
        "results/taxonomy/{sample}.metaxa2.spot_taxonomy.tsv"
    params:
        score = config["metaxa2"]["score_cutoff"]
    run:
        fhlog = open(log[0], 'w')
        # Read taxonomy table
        taxdf = pd.read_csv(input.tax, sep="\t", index_col=0, header=None,
            names=["ID", "Classification", "Identity", "Length", "Score"])
        # Filter taxonomy to assignments above threshold
        taxdf = taxdf.loc[taxdf.Score>=params.score]

        # Read barcodes
        barcodes = pd.read_csv(input.barcodes, sep="\t",index_col=0,header=0,
            dtype=str,names=["barcode", "x", "y"])
        # Make coord column and create renaming dictionary
        barcodes["coord"] = barcodes[["x", "y"]].agg("x".join,axis=1)
        barcodes = barcodes.drop(["x", "y"],axis=1).to_dict()["coord"]

        # Read file with mapping of read to UMI and barcode
        mapdf = pd.read_csv(input.mapfile,sep="\t",index_col=0,header=0,
            names=["read_id", "umi", "spot"])

        # Merge taxonomy and mapdf
        spot_taxdf = pd.merge(taxdf,mapdf.drop("umi",axis=1),how="inner",left_index=True,
            right_index=True)
        spot_taxdf.fillna("Unassigned",inplace=True)
        # groupby and count assignments per spot
        spot_taxonomy_counts = spot_taxdf.groupby(["spot", "Classification"]).count().loc[:, "Score"]
        spot_taxonomy_counts = pd.DataFrame(spot_taxonomy_counts.reset_index())
        spot_taxonomy_counts = pd.pivot_table(pd.DataFrame(spot_taxonomy_counts),columns="spot",index="Classification")
        spot_taxonomy_counts = spot_taxonomy_counts["Score"].fillna(0)
        index = list(spot_taxonomy_counts.sum(axis=1).sort_values(ascending=False).index)
        # Transpose dataframe and change index names
        spot_taxonomy_counts = spot_taxonomy_counts.rename(columns=barcodes)
        spot_taxonomy_counts.loc[index].T.to_csv(output[0], sep="\t")