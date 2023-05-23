#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Output basic coordinates of each seq in a FASTA file (all starting from 0 and going to the length of the sequence). The scaffold is taken to be the FASTA headerline. If --subset is set then only the sequence ids specified in that file will be retained.",

    epilog='''Usage example:

    python fasta_to_basic_bed.py --input FASTA --subset TABLE --output BED

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FASTA", type=str,
                        help="Path to input FASTA", required=True)

    parser.add_argument("-s", "--subset", metavar="IDS", type=str,
                        help="Path to optional input file with ids to keep - one per line.",
                        required=False, default = None)

    parser.add_argument("-o", "--output", metavar="BED", type=str,
                        help="Path to output BED", required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.input, cut_header = False)

    if args.subset is not None:
        subset_ids = set()
        with open(args.subset, "r") as in_subset:
            for subset_line in in_subset:
                subset_ids.add(subset_line.rstrip())

        for subset_id in subset_ids:
            if subset_id not in seqs.keys():
                sys.exit("Error specified subset id " + subset_id + " not in FASTA")

    with open(args.output, "w") as output_bed:
        for seqid, sequence in seqs.items():
            if args.subset is not None and seqid not in subset_ids:
                continue
            print(seqid + "\t" + str(0) + "\t" + str(len(sequence)), file = output_bed)


if __name__ == '__main__':
    main()
