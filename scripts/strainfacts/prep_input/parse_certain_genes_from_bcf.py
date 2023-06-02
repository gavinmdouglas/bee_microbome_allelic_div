#!/usr/bin/python3

import argparse
import os
import sys
import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
from functions.io_utils import read_ids

class bed_coor:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


def main():

    parser = argparse.ArgumentParser(

            description='Parse subset of lines from bcf. This was useful for testing.',

    epilog='''Usage example:

    python parse_certain_genes_from_bcf.py --input INPUT --gene_set GENES

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='INPUT', type=str, required=True,
                        help='Path to input BCF of all samples.')

    parser.add_argument('-c', '--gene_set', metavar='GENES', type=str, required=True,
                        help='Print out all lines that match genes in this file.')

    parser.add_argument('-b', '--bed', metavar='BED', type=str, required=True,
                        help='Path to reference bed file with gene names and coordinates of interest.')


    args = parser.parse_args()

    bcf_in = pysam.VariantFile(args.input)

    gene_set = read_ids(args.gene_set)

    gene_info = dict()

    with open(args.bed, 'r') as bed_filehandle:

        for bed_line in bed_filehandle:

            bed_line_split = bed_line.split()

            in_gene = bed_line_split[0]

            if in_gene in gene_set:

                gene_info[in_gene] = bed_coor(int(bed_line_split[1]),
                                              int(bed_line_split[2]))


    print(str(bcf_in.header), end = '')

    for gene in gene_set:

        for rec in bcf_in.fetch(contig = gene,
                                start = gene_info[gene].start,
                                stop = gene_info[gene].stop):

            rec_split = str(rec).split()

            # Skip line if it is not polymorphic, or if it is multi-allelic.
            alt_allele = rec_split[4]
            if alt_allele == '.' or ',' in alt_allele:
                continue

            print(str(rec), end = '')


if __name__ == '__main__':
    main()
