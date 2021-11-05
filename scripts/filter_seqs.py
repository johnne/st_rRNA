#!/usr/bin/env python

from Bio import SeqIO
from gzip import open as gzopen


def check_seq(seq):
    test_string = "".join(sorted(set(seq)))
    for nt in list(test_string):
        if nt.upper() not in ["A","C","G","T"]:
            return False
    return True


def main(sm):
    openfn = open
    if sm.input[0].endswith(".gz"):
        openfn = gzopen
    records = []
    with openfn(sm.input[0], 'rt') as fhin, openfn(sm.output[0], 'wt') as fhout:
        for record in SeqIO.parse(fhin, sm.params.seq_type):
            if check_seq(record.seq):
                records.append(record)
    with openfn(sm.output[0], 'wt') as fhout:
        SeqIO.write(records, fhout, sm.params.seq_type)


if __name__ == "__main__":
    main(snakemake)