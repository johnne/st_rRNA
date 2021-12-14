import os
from snakemake.utils import validate

include: "scripts/common.py"
localrules:
    multiqc,
    report,
    download_sortmerna_db,
    download_taxref,
    download_vsearch,
    download_speciesref,
    extract_R1,
    sample,
    gather_spots,
    filter_species_fasta,
    filter_seqs,
    reformat_taxref,
    spot_taxonomy

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
        expand("results/taxonomy/{sample}.{subunit}.{taxtool}.spot_taxonomy.tsv",
            sample = samples.keys(), subunit = config["subunits"],
            taxtool = config["taxtool"]),
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
        "results/rRNA/{subunit}/{sample}.R2.fastq.gz"
    conda: "envs/sortmerna.yml"
    log:
        runlog="results/logs/rRNA/{subunit}/{sample}.log",
        reportlog="results/rRNA/{subunit}/{sample}.log"
    params:
        evalue = config["sortmerna"]["evalue"],
        workdir = "$TMPDIR/rRNA/{subunit}.{sample}.wd",
        aligned = "{sample}.R2",
        R2 = "$TMPDIR/rRNA/{subunit}.{sample}.wd/R2.fastq",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        ref_string = lambda wildcards, input: " ".join([f"--ref {x}" for x in input.db])
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 4
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=temp; fi
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        gunzip -c {input.R2} > {params.R2}
        sortmerna -e {params.evalue} --threads {threads} --workdir {params.workdir} --fastx \
            --reads {params.R2} {params.ref_string} \
            --aligned {params.workdir}/{params.aligned} > {log.runlog} 2>&1
        rm {params.R2}
        gzip {params.workdir}/*.fastq
        mv {params.workdir}/*.gz {params.outdir}
        mv {params.workdir}/{wildcards.sample}.R2.log {log.reportlog}
        rm -rf {params.workdir}
        """

rule extract_R1:
    """
    Extracts R1 sequences for reads identified as rRNA
    """
    input:
        R1 = lambda wildcards: samples[wildcards.sample]["R1"],
        R2 = rules.sortmerna.output
    output:
        R1 = "results/rRNA/{subunit}/{sample}.R1.fastq.gz"
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
        R1 = rules.extract_R1.output,
        R2 = rules.sortmerna.output
    output:
        R2 = temp("results/rRNA/{subunit}/spots/{barcode}/{sample}.R2.subsampled.fastq.gz"),
        f = temp("results/rRNA/{subunit}/spots/{barcode}/{sample}.map.tsv"),
        stats = temp("results/rRNA/{subunit}/spots/{barcode}/{sample}.stats")
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
        runtime = lambda wildcards, attempt: attempt ** 2 * 10
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

rule filter_seqs:
    """
    Removes sequences with characters other than ["A","C","G","T"]
    """
    input:
        rules.sample_spots.output.R2
    output:
        "results/rRNA/{subunit}/spots/{barcode}/{sample}.R2.subsampled.filtered.fastq.gz"
    log:
        "results/logs/{subunit}/{sample}.{barcode}.filter_nonDNA.log"
    conda: "envs/seqkit.yml"
    params:
        seq_type = "fastq"
    script: "scripts/filter_seqs.py"

rule assignTaxonomySpot:
    """
    Runs assignTaxonomy on subsampled and filtered sequences for a spot
    """
    input:
        refFasta = "resources/DADA2/{subunit}.assignTaxonomy.reformat.fasta",
        spFasta = "resources/DADA2/{subunit}.addSpecies.filtered.fasta",
        seqs = rules.filter_seqs.output[0]
    output:
        taxdf = "results/rRNA/{subunit}/spots/{barcode}/{sample}.assignTaxonomy.tsv",
        bootdf = "results/rRNA/{subunit}/spots/{barcode}/{sample}.assignTaxonomy.bootstrap.tsv"
    log:
        "results/rRNA/{subunit}/spots/{barcode}/{sample}.assignTaxonomy.log"
    params:
        minBoot = config["assignTaxonomy"]["minBoot"],
        outputBootstraps=config["assignTaxonomy"]["outputBootstraps"],
        tryRC=config["assignTaxonomy"]["tryRC"],
        taxLevels=lambda wildcards: config["assignTaxonomy"]["ranks"][
            wildcards.subunit]
    conda: "envs/dada2.yml"
    resources:
        runtime = lambda wildcards,attempt: attempt ** 2 * 60 * 10
    threads: config["assignTaxonomy"]["threads"]
    script:
        "scripts/assignTaxa.R"


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
        tax = expand("results/rRNA/{{subunit}}/spots/{barcode}/{{sample}}.assignTaxonomy.tsv",
            barcode = barcodes.keys()),
        boot = expand("results/rRNA/{{subunit}}/spots/{barcode}/{{sample}}.assignTaxonomy.bootstrap.tsv",
            barcode = barcodes.keys()),
        tsv = expand("results/rRNA/{{subunit}}/spots/{barcode}/{{sample}}.map.tsv",
            barcode = barcodes.keys()),
        stats = expand("results/rRNA/{{subunit}}/spots/{barcode}/{{sample}}.stats",
            barcode = barcodes.keys())
    output:
        tax = "results/rRNA/{subunit}/{sample}.assignTaxonomy.tsv",
        boot = "results/rRNA/{subunit}/{sample}.assignTaxonomy.bootstrap.tsv",
        tsv = "results/rRNA/{subunit}/{sample}.map.tsv",
        stats = "results/rRNA/{subunit}/{sample}.stats.tsv"
    params:
        tmp1 = "$TMPDIR/{subunit}.{sample}.sampled.R2.rRNA.fastq.gz",
        tmp2 = "$TMPDIR/{subunit}.{sample}.mapfile.tsv",
        tmp3 = "$TMPDIR/{subunit}.{sample}.stats.tsv"
    run:
        # Concatenate tsv and stats files
        concat_files(input.tsv, output.tsv)
        concat_files(input.stats, output.stats)
        concat_df(input.tax, output.tax)
        concat_df(input.boot, output.boot)

## Sample target rule
rule sample:
    input:
        expand("results/rRNA/{subunit}/{sample}.assignTaxonomy.tsv",
            subunit = config["subunits"], sample = samples.keys()),
        expand("results/rRNA/{subunit}/{sample}.map.tsv",
            subunit = config["subunits"], sample = samples.keys()),
        expand("results/rRNA/{subunit}/{sample}.stats.tsv",
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

rule download_taxref:
    """
    Downloads a DADA2-ready train set for use with the assignTaxonomy function
    """
    output:
        "resources/DADA2/{subunit}.assignTaxonomy.fasta"
    log:
        "results/logs/DADA2/download.{subunit}.assignTaxonomy.log"
    params:
        url=lambda wildcards: config["assignTaxonomy"]["taxdb"][wildcards.subunit]
    shell:
        """
        curl -L -o {output[0]}.gz {params.url} > {log} 2>&1
        gunzip {output[0]}.gz
        """

rule reformat_taxref:
    input:
        "resources/DADA2/{subunit}.assignTaxonomy.fasta"
    output:
        "resources/DADA2/{subunit}.assignTaxonomy.reformat.fasta"
    log:
        "results/logs/{subunit}.reformat_taxref.log"
    params:
        num_ranks = lambda wildcards: len(config["assignTaxonomy"]["ranks"][wildcards.subunit])
    conda: "envs/seqkit.yml"
    script:
        "scripts/reformat_taxref.py"


def postprocess(wildcards):
    if wildcards.subunit == "18S":
        s = " | sed '/^>/s/>\([^;]*\);.*,s:\(.*\)/>\\1 \\2/' | sed 's/_/ /g' | sed 's/ \([A-Z]\) /_\\1 /' > "
    else:
        s = " > "
    return s

rule download_speciesref:
    """
    Downloads a DADA2-ready train set for use with addSpecies function
    """
    output:
        "resources/DADA2/{subunit}.addSpecies.fasta"
    log:
        "results/DADA2/download.{subunit}.addSpecies.log"
    params:
        url = lambda wildcards: config["assignTaxonomy"]["speciesdb"][wildcards.subunit],
        postprocess = lambda wildcards: postprocess(wildcards)
    shell:
        """
        curl -L -o {output[0]}.gz {params.url} > {log} 2>&1
        gunzip -c {output[0]}.gz {params.postprocess} {output[0]}
        rm {output[0]}.gz
        """

rule filter_species_fasta:
    """
    Removes sequences containing non DNA characters
    """
    input:
        "resources/DADA2/{subunit}.addSpecies.fasta"
    output:
        "resources/DADA2/{subunit}.addSpecies.filtered.fasta"
    log:
        "resources/DADA2/{subunit}.filter_species.log"
    conda: "envs/seqkit.yml"
    params:
        seq_type = "fasta"
    script:
        "scripts/filter_seqs.py"

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
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60
    conda: "envs/vsearch.yml"
    threads: 10
    shell:
        """
        vsearch --threads {threads} --sintax {input.seqs} --db {input.db} \
            --sintax_cutoff {params.cutoff} --tabbedout {output[0]} > {log} 2>&1 
        """

def split_common_reads(taxf, bootf):
    tax16f, tax18f = taxf
    boot16f, boot18f = bootf
    tax16 = pd.read_csv(tax16f,header=0,sep="\t",index_col=0)
    tax18 = pd.read_csv(tax18f,header=0,sep="\t",index_col=0)
    boot16 = pd.read_csv(boot16f,header=0,sep="\t",index_col=0)
    boot18 = pd.read_csv(boot18f,header=0,sep="\t",index_col=0)
    common_reads = set(tax18.index).intersection(tax16.index)
    common_ranks = set(tax18.columns).intersection(tax16.columns)
    common_ranks = [rank for rank in tax18.columns if rank in common_ranks]
    drop = {"16S": [], "18S": []}
    for r in common_reads:
        tax18s = tax18.loc[r]
        tax16s = tax16.loc[r]
        boot_18s = boot18.loc[r]
        boot_16s = boot16.loc[r]
        # Start from most specific rank
        for rank in reversed(common_ranks):
            if tax18s[rank]==tax18s[rank] and tax16s[rank]!=tax16s[rank]:
                drop["16S"].append(r)
                break
            elif tax18s[rank]!=tax18s[rank] and tax16s[rank]==tax16s[rank]:
                drop["18S"].append(r)
                break
            elif tax18s[rank]!=tax18s[rank] and tax16s[rank]!=tax16s[rank]:
                continue
            if boot_18s[rank] >= boot_16s[rank]:
                drop["16S"].append(r)
                break
            else:
                drop["18S"].append(r)
                break
    return drop

rule spot_taxonomy:
    input:
        tax = expand("results/rRNA/{subunit}/{{sample}}.assignTaxonomy.tsv",
            subunit = config["subunits"]),
        boot = expand("results/rRNA/{subunit}/{{sample}}.assignTaxonomy.bootstrap.tsv",
            subunit = config["subunits"]),
        mapfile = expand("results/rRNA/{subunit}/{{sample}}.map.tsv",
            subunit = config["subunits"]),
        barcodes = config["barcode_list"]
    output:
        "results/taxonomy/{sample}.16S.assignTaxonomy.spot_taxonomy.tsv",
        "results/taxonomy/{sample}.18S.assignTaxonomy.spot_taxonomy.tsv"
    log:
        "results/logs/{sample}.spot_taxonomy.log"
    params:
        filter_rank = lambda wildcards: samples[wildcards.sample]["filter_rank"],
        minBoot = 75,
        subunits = config["subunits"]
    run:
        logfh = open(log[0], 'w')
        # Read barcodes
        logfh.write(f"Reading barcodes from {input.barcodes}\n")
        barcodes = pd.read_csv(input.barcodes, sep="\t", index_col=0, header=0,
        dtype=str, names=["barcode", "x", "y"])
        # Make coord column and create renaming dictionary
        barcodes["coord"] = barcodes[["x", "y"]].agg("x".join,axis=1)
        barcodes = barcodes.drop(["x","y"], axis=1).to_dict()["coord"]

        # If running with more than one subunit, identify common reads and
        # split them depending on bootstrap to the different subunits
        drop = {"16S": [], "18S": []}
        if len(config["subunits"]) > 1:
            drop = split_common_reads(sorted(input.tax), sorted(input.boot))
            common_reads = len(drop['16S'])+len(drop['18S'])
            logfh.write(f"Found a total of {common_reads} reads with matches to both subunits\n")
        for i, taxfile in enumerate(sorted(input.tax)):
            bootfile = sorted(input.boot)[i]
            mapfile = sorted(input.mapfile)[i]
            outfile = sorted(output)[i]
            subunit = taxfile.split("/")[-2]
            logfh.write(f"Taxfile: {taxfile}\nBootstraps: {bootfile}\nMapfile: {mapfile}\nOutfile: {outfile}\nSubunit: {subunit}\n")
            taxdf = pd.read_csv(taxfile, sep="\t", header=0, index_col=0)
            taxdf.drop(drop[subunit])
            if len(drop[subunit])>0:
                logfh.write(f"Dropped {len(drop[subunit])} reads that had better assignment for other subunit\n")
            bootdf = pd.read_csv(bootfile, sep="\t", header=0, index_col=0)
            bootdf.drop(drop[subunit])
            bootdf = bootdf.loc[bootdf[params.filter_rank] >= params.minBoot]
            reads = set(bootdf.index).intersection(taxdf.index)
            taxdf = taxdf.loc[reads]
            ranks = list(taxdf.columns)
            mapdf = pd.read_csv(mapfile, sep="\t", index_col=0, header=0,
                names=["read_id", "umi", "spot"])
            spot_taxdf = pd.merge(taxdf, mapdf.drop("umi", axis=1), how="inner", left_index=True,
                right_index=True)
            spot_taxdf.fillna("Unassigned", inplace=True)
            # add taxonomy column
            spot_taxdf["taxonomy"] = spot_taxdf[ranks].agg(";".join,axis=1)
            # groupby and count
            spot_taxonomy_counts = spot_taxdf.groupby(["spot", "taxonomy"]).count().loc[:, ranks[0]]
            spot_taxonomy_counts = pd.DataFrame(spot_taxonomy_counts.reset_index().rename(columns={ranks[0]: 'n'}))
            spot_taxonomy_counts = pd.pivot_table(spot_taxonomy_counts, columns="spot", index="taxonomy")
            spot_taxonomy_counts = spot_taxonomy_counts["n"].fillna(0)
            index = list(spot_taxonomy_counts.sum(axis=1).sort_values(ascending=False).index)
            # Transpose dataframe and change index names
            spot_taxonomy_counts = spot_taxonomy_counts.rename(columns=barcodes)
            spot_taxonomy_counts.loc[index].T.to_csv(outfile, sep="\t")
        logfh.close()
