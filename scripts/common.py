import pandas as pd


def rRNA_fasta(wildcards):
    input = []
    for f in rRNA[wildcards.subunit]:
        input.append(f"resources/sortmerna/{f}")
    return input


def parse_samples(f):
    df = pd.read_csv(f, index_col=0, sep="\t", header=0)
    return df.to_dict(orient="index")


def read_barcodes(f):
    df = pd.read_csv(f, index_col=0, sep="\t", header=None,
                     names=["barcode","x","y"])
    return df.to_dict(orient="index")