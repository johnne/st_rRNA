#!/usr/bin/env python

from Bio.SeqIO import parse


def main(sm):
    with open(sm.input[0], 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for record in parse(fhin, "fasta"):
            items = (record.id).split(";")
            if len(items) == 7 and items[0] == items[1]:
                recid = ";".join(items[1:])
                fhout.write(f">{recid}\n{record.seq}\n")


if __name__ == "__main__":
    main(snakemake)