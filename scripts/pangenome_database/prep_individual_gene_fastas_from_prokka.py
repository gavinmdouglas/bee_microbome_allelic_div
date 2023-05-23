#!/usr/bin/python3

import argparse
import os
import sys
import textwrap
import gzip

from functions.io_utils import read_fasta, write_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Create FASTA files for each separate gene in prokka output.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="GENE_FASTA", type=str,
                        help="Path to FASTA file containing nucleotide sequences of genes.",
                        required=True)

    parser.add_argument("-p", "--prefix", metavar="PREFIX", type=str,
                        help="Prefix to add to beginning of file names.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="FOLDER", type=str,
                        help="Path to output folder, which will contain a FASTA for each gene.",
                        required=True)

    args = parser.parse_args()

    gene_seqs = read_fasta(args.fasta, cut_header = True)

    os.mkdir(args.output)        

    for gene in gene_seqs:

        gene_subset = { gene : gene_seqs[gene] }

        outfile = args.output + '/' + args.prefix + '_' + gene + '.fa'

        write_fasta(gene_subset, outfile)


if __name__ == '__main__':
    main()
