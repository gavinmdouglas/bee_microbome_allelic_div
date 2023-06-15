#!/usr/bin/python3

import argparse
import sys
import os
from evolution_metric_functions import (codon_to_aa,
                                        pairwise_dnds,
                                        read_fasta) 
import itertools
import pandas as pd
import numpy as np

def main():

    parser = argparse.ArgumentParser(

        description='Read in FASTA file of CDS sequences and compute dN, dS, and dN/dS for '
                    'each pairwise comparison. Will optionally ignore the first and last X bases. '
                    'NOTE: This approach does not correct for mutational biases (e.g., different '
                    'rates of transition vs transversion mutations), as implemented in Li 1993. '
                    'Note that this script is a modified version of \'mean_pairwise_dnds.py\' (v1.1.0) in '
                    'Gavin Douglas\' \'handy_pop_gen\' GitHub repository (under GPL 3 license).',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='Path to codon-aligned FASTA.', required=True)

    parser.add_argument('-o', '--output', metavar='OUTPUT', type=str,
                        help='Path to output table of all pairwise comparisons.', required=True)

    parser.add_argument('--bases_to_trim', metavar='INT', type=int,
                        help='Number of positions to trim from the start and end of the DNA alignment.',
                        required=False, default = 0)

    args = parser.parse_args()

    # Read in input sequences.
    seqs = read_fasta(args.input, convert_upper=True)

    if args.bases_to_trim > 0:
        trimmed_seqs = dict()
        for seq_id, seq_val in seqs.items():
            trimmed_seqs[seq_id] = seq_val[args.bases_to_trim: -args.bases_to_trim]

        seqs = trimmed_seqs

    if len(seqs.keys()) == 0:
        sys.exit("Stopping: no sequences were found in the input file: " + args.input)
    elif len(seqs.keys()) == 1:
        sys.exit("Stopping: there is only one input sequence present: " + args.input)

    # Get mean dn, ds, and dn/ds based on all pairwise comparisons.
    pairwise_combos = list(itertools.combinations(seqs, 2))

    dnds = pd.DataFrame(columns=['n_subs', 'n_sites',
                                 's_subs', 's_sites',
                                 'dn', 'ds', 'dnds'],
                        index=pd.MultiIndex.from_tuples(pairwise_combos,
                                                        names=('seq1', 'seq2')))
    for combo in pairwise_combos:
        dnds.loc[combo, :] = list(pairwise_dnds(seqs[combo[0]],
                                                seqs[combo[1]]))

    dnds.to_csv(args.output,
                sep='\t', header=True, na_rep='NA')


if __name__ == '__main__':
    main()
