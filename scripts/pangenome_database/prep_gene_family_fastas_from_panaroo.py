#!/usr/bin/python3

import argparse
import gzip
from collections import defaultdict
import os
from functions.io_utils import read_fasta, write_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Create FASTA per all genes in gene family (based on specified set) by parsing Panaroo output. Optionally will trim bases from the start and end of each gene.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-d", "--data", metavar="GENE_DATA", type=str,
                        help="Path to file containing info per Panaroo gene (gene_data.csv.gz).",
                        required=True)
    
    parser.add_argument("-p", "--presence", metavar="PANAROO_PRESENCE", type=str,
                        help="Path to Panaroo presence/absence table (gene_presence_absence.csv.gz).",
                        required=True)

    parser.add_argument("-g", "--genes", metavar="GENE_LIST", type=str,
                        help="Path to file with gene families to parse (one per line).",
                        required=True)

    parser.add_argument("-s", "--species_prefix", metavar="SPECIES_NAME", type=str,
                        help="Species name to prefix to gene family in output filenames (and to remove from input gene family IDs in --genes).",
                        required=True)

    parser.add_argument("-t", "--trim_length", metavar="INT", type=int, default=0,
                        help="Number of bases to remove from the start and end of each gene.")

    parser.add_argument("-m", "--min_length", metavar="INT", type=int, default=100,
                        help="Min. length of all genes in a gene family for it to be included in output.")

    # Add flag to specify whether genes with multi copies in the same genome should be included.
    parser.add_argument("--singlecopy_core_only", action="store_true",
                        help="Exclude genes that are not single copy and found in every genome.")

    parser.add_argument("-o", "--output", metavar="FOLDER", type=str,
                        help="Path to output folder, which will contain a FASTA for each gene family.",
                        required=True)

    args = parser.parse_args()

    # Read in gene family IDs from specified file.
    gene_families = set()
    with open(args.genes, 'r') as gene_fh:
        for gene_line in gene_fh:
            raw_gf = gene_line.rstrip()
            raw_gf = raw_gf.replace(args.species_prefix + '_', '')
            gene_families.add(raw_gf)

    # Parse Panaroo presence/absence table and prep mapfile of gene IDs to gene families and genomes.
    gene_family_map = {}
    genome_map = {}
    if args.singlecopy_core_only:
        multi_copies_to_exclude = set()
    with gzip.open(args.presence, 'rt') as presence_fh:
        headerline = presence_fh.readline()
        genomes = headerline.rstrip().split(',')[3:]
        for line in presence_fh:
            line = line.rstrip()
            line = line.split(',')
            gene_family = line[0]

            if gene_family not in gene_families:
                continue

            for i in range(3, len(line)):
                if line[i] != '':
                    if ';' in line[i] and args.singlecopy_core_only:
                        multi_copies_to_exclude.add(gene_family)
                    for g in line[i].split(';'):
                        genome_map[g] = genomes[i - 3]
                        gene_family_map[g] = gene_family

    # Parse Panaroo gene data and write FASTA per gene family.
    gene_family_seqs = defaultdict(dict)
    gene_family_length_fails = set()
    with gzip.open(args.data, 'rt') as data_fh:
        for line in data_fh:
            line = line.rstrip()
            line = line.split(',')

            gene_id = line[3]

            if gene_id not in gene_family_map.keys():
                continue

            gene_family = gene_family_map[gene_id]
            genome = genome_map[gene_id]
            seq = line[5]

            if args.trim_length > 0:
                seq = seq[args.trim_length:-args.trim_length]

            if len(seq) < args.min_length:
                gene_family_length_fails.add(gene_family)

            if gene_family not in gene_family_seqs:
                gene_family_seqs[gene_family] = dict()

            gene_family_seqs[gene_family][genome + '||' + gene_id] = seq

    # Create output folder if it does not exist.
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Write out FASTAs.
    for gf in gene_family_seqs.keys():
        if gf in gene_family_length_fails or gf in multi_copies_to_exclude:
            continue
            
        if args.singlecopy_core_only and len(gene_family_seqs[gf]) != len(genomes):
            continue
 
        write_fasta(gene_family_seqs[gf],
                    args.output + '/' + args.species_prefix + '_' + gf + '.fa')

if __name__ == '__main__':
    main()
