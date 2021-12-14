#!/usr/bin/env python

import sys
from Bio import SeqIO
import random
random.seed(99)
from argparse import ArgumentParser


def read_fastq(f, umi_start=18, umi_stop=25):
    """
    This function performs a very naive clustering of reads by UMI
    by simply adding UMIs as keys in a dictionary and the corresponding
    reads to a list

    :param f: fastq file
    :param umi_start: start position of UMI in each read
    :param umi_stop: stop position of UMI in each read
    :return: dictionary with UMIs as keys and a list of reads as value
    """
    umi_dict = {}
    for read in SeqIO.parse(f, "fastq"):
        # Extract the umi based on umi start and stop
        umi = read.seq[umi_start:umi_stop]
        # Attempt to add read to a list with this umi as key
        try:
            umi_dict[umi].append(read)
        # If umi has not been added, create a new list and add the read
        except KeyError:
            umi_dict[umi] = [read]
    return umi_dict


def get_umi_rep(umi_dict):
    """
    Get one random representative for each UMI

    :param umi_dict: UMI clustered reads
    :return: dictionary with one read per UMI
    """
    umi_reps = {}
    # Iterate the key->value (UMI -> list of reads) pairs
    # and select 1 read at random from the list corresponding to the UMI
    for umi, reads in umi_dict.items():
        umi_reps[umi] = random.sample(reads, k=1)[0]
    return umi_reps


def sample_spot(umi_reps, X=100):
    """
    Gathers all UMI clustered reads into a list, then samples X reads
    from that list

    :param umi_dict: UMI clustered reads
    :param X: number of reads to sample from the spot
    :return: tuple of sampled reads and a dictionary mapping reads to UMIs
    """
    reads = []
    read_to_umi = {}
    # Iterate umi->representative read key/value pairs
    for umi, read in umi_reps.items():
        # Store a link of read ID -> UMI map
        read_to_umi[read.id] = umi
        # Add representative read to list
        reads.append(read)
    # Finally, sample X reads from the list of representative reads
    sampled_reads = random.sample(reads, k=min(X, len(reads)))
    return sampled_reads, read_to_umi


def write_mapfile(f, barcode, sampled_reads, read_to_umi):
    """
    Writes a mapfile of sampled reads to their UMI and (optionally)
    their barcode.

    :param f: Output file
    :param barcode: optional barcode to append
    :param read_to_umi: dictionary of reads to UMIs
    :return:
    """

    if barcode is None:
        b = ""
    else:
        b = f"\t{barcode}"
    with open(f, 'w') as fhout:
        for read in sampled_reads:
            umi = read_to_umi[read.id]
            fhout.write(f"{read.id}\t{umi}{b}\n")


def write_read_ids(sampled_reads):
    """
    Writes read ids to file

    :param sampled_reads: list of sampled reads
    :return:
    """
    with sys.stdout as fhout:
        for read in sampled_reads:
            fhout.write(f"{read.id}\n")


def write_stats(umi_dict, statsfile, barcode):
    """
    Write a table with statistics of UMIs and number of reads per UMI found
    in each spot (barcode)

    :param umi_dict: Dictionary of UMIs and their reads
    :param statsfile: Output file
    :param barcode: Barcode (spot) being sampled
    :return:
    """
    if barcode:
        barcode_string = f"\t{barcode}\n"
    else:
        barcode_string = "\n"
    with open(statsfile, 'w') as fhout:
        for umi, reads in umi_dict.items():
            # Write UMI<tab>number of reads<barcode>
            fhout.write(f"{umi}\t{len(reads)}{barcode_string}")


def main(args):
    # Get a dictionary of unique UMIs (keys) and a list of reads with
    # those UMIs (values)
    umi_dict = read_fastq(args.R1, args.start, args.stop)
    # Write to file
    if args.stats:
        write_stats(umi_dict, args.stats, args.barcode)
    # Extract 1 representative sequence for each UMI
    umi_reps = get_umi_rep(umi_dict)
    # Sample args.num_reads from representative reads
    sampled_reads, read_to_umi = sample_spot(umi_reps, args.num_reads)

    write_read_ids(sampled_reads)
    if args.mapfile is not None:
        write_mapfile(args.mapfile, args.barcode, sampled_reads, read_to_umi)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("R1", help="Read 1 file")
    parser.add_argument("--stats",
                        help="Write number of reads per UMI for spot")
    parser.add_argument("--num_reads", type=int, default=100,
                        help="Number of reads to sample from a spot")
    parser.add_argument("--start", type=int, default=18,
                        help="Start position of UMI")
    parser.add_argument("--stop", type=int, default=25,
                        help="Stop position of UMI")
    parser.add_argument("--barcode", type=str,
                        help="Barcode string to add to mapfile")
    parser.add_argument("--mapfile", type=str,
                        help="Tab-separated file mapping read ids to "
                             "barcode and UMI")
    args = parser.parse_args()
    main(args)