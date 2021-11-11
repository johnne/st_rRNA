#!/usr/bin/env python

from Bio.SeqIO import parse


def get_len(f_in, f_out, l=8):
    lengths = {}
    with open(f_in, 'r') as fhin, open(f_out, 'w') as fhout:
        for record in parse(fhin, "fasta"):
            items = (record.id).split(";")
            try:
                lengths[len(items)]+=1
            except KeyError:
                lengths[len(items)]=1
            if len(items) == l:
                fhout.write(f">{record.description}\n{record.seq}\n")
    return lengths




def main(sm):
    num_ranks = sm.params.num_ranks
    with open(sm.input[0], 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for record in parse(fhin, "fasta"):
            record.id = (record.id).rstrip(";")
            items = (record.id).split(";")
            if items[0] == items[1]:
                items = items[1:]
            recid = ";".join(items)
            if len(items) == num_ranks:
                fhout.write(f">{recid}\n{record.seq}\n")


if __name__ == "__main__":
    main(snakemake)