#!/usr/bin/env python

from Bio.SeqIO import parse


def check_seq(seq):
    test_string = "".join(sorted(set(seq)))
    for nt in list(test_string):
        if nt.upper() not in ["A","C","G","T"]:
            return False
    return True


def main(sm):
    with open(sm.input[0], 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for record in parse(fhin, "fasta"):
            if check_seq(record.seq):
                fhout.write(f">{record.description}\n{str(record.seq)}\n")


if __name__ == "__main__":
    main(snakemake)