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


def spot_taxonomy(sm):
    import pandas as pd
    # Read taxonomy table
    taxdf = pd.read_csv(sm.input.tax, sep="\t", index_col=0, header=None,
                        names=["ID", "Classification", "Identity", "Length",
                               "Score"])
    # Filter taxonomy to assignments above threshold
    taxdf = taxdf.loc[taxdf.Score >= sm.params.score]

    # Read barcodes
    barcodes = pd.read_csv(sm.input.barcodes, sep="\t", index_col=0, header=0,
                           dtype=str, names=["barcode", "x", "y"])
    # Make coord column and create renaming dictionary
    barcodes["coord"] = barcodes[["x", "y"]].agg("x".join, axis=1)
    barcodes = barcodes.drop(["x", "y"], axis=1).to_dict()["coord"]

    # Read file with mapping of read to UMI and barcode
    mapdf = pd.read_csv(sm.input.mapfile, sep="\t", index_col=0, header=0,
                        names=["read_id", "umi", "spot"])

    # Merge taxonomy and mapdf
    spot_taxdf = pd.merge(taxdf, mapdf.drop("umi", axis=1), how="inner",
                          left_index=True, right_index=True)
    spot_taxdf.fillna("Unassigned", inplace=True)
    # groupby and count assignments per spot
    spot_taxonomy_counts = spot_taxdf.groupby(["spot", "Classification"]).count().loc[:, "Score"]
    spot_taxonomy_counts = pd.DataFrame(spot_taxonomy_counts.reset_index())
    spot_taxonomy_counts = pd.pivot_table(pd.DataFrame(spot_taxonomy_counts),
                                          columns="spot",
                                          index="Classification")
    spot_taxonomy_counts = spot_taxonomy_counts["Score"].fillna(0)
    index = list(spot_taxonomy_counts.sum(axis=1).sort_values(ascending=False).index)
    # Transpose dataframe and change index names
    spot_taxonomy_counts = spot_taxonomy_counts.rename(columns=barcodes)
    spot_taxonomy_counts.loc[index].T.to_csv(sm.output[0], sep="\t")


def main(sm):
    toolbox = {'spot_taxonomy': spot_taxonomy}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)