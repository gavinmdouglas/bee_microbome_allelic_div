#!/usr/bin/python3

import argparse
import sys
import textwrap
import gzip
from collections import defaultdict
import os


def main():

    parser = argparse.ArgumentParser(
            description="Parse StrainFacts output files, and supplementary files, specifically for concatenated "
                        "core genes/core genome, which will produce a core genome FASTA file, with one sequence per inferred strain. "
                        "Note that INDELs will be ignored in this output. "
                        "Also, the leading and trailing 102 bp from each core gene will be trimmed off.",
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--species',
                        metavar='SPECIES', type=str,
                        help="Species name to parse.",
                        required=True)
    
    parser.add_argument('--dataset',
                        metavar='DATASET', type=str,
                        help="Datset name to parse.",
                        required=True)
    
    parser.add_argument('-c', '--core_genes',
                        metavar='GENE_LIST', type=str,
                        help="Path to file containing core genes to parse: one per line.",
                        required=True)

    parser.add_argument('-r', '--ref_fasta_file',
                        metavar='REF_FASTA', type=str,
                        help="Path to gzipped FASTA file containing all representative gene sequences (separately and untrimmed).",
                        required=True)

    parser.add_argument('--strain_presence_file',
                        metavar='STRAIN_PRESENCE', type=str,
                        help="Path to gzipped table with mapping of species names to strains numbers called as present across samples.",
                        required=True)

    parser.add_argument('-s', '--sites_file',
                        metavar='SITES_FILE', type=str,
                        help="Path to gzipped StrainFacts-output sites file (note that this includes the gene names as well).",
                        required=True)

    parser.add_argument('-g', '--geno_file',
                        metavar='GENO_FILE', type=str,
                        help="Path to gzipped StrainFacts-output genotype file.",
                        required=True)

    parser.add_argument('-p', '--output_seq_prefix',
                        metavar='SEQNAME_PREFIX', type=str,
                        help="String to affix to start of output sequence number (delimited by '_')",
                        required=False, default='strain')

    parser.add_argument('-f', '--output_fasta',
                        metavar='OUTPUT_FASTA', type=str,
                        help="Path to output concatenated core genome FASTA.",
                        required=True)

    parser.add_argument('--output_individual_fastas',
                        metavar='PATH', type=str,
                        help="Path to folder for individual core gene sequences for all strains, one per strain.",
                        required=True)

    parser.add_argument('-m', '--output_map',
                        metavar='OUTPUT_MAP', type=str,
                        help="Path to output mapfile (positions of each core gene in concatenated core genome).",
                        required=True)

    args = parser.parse_args()

    species = args.species
    dataset = args.dataset
    ref_fasta_file = args.ref_fasta_file
    strain_presence_file = args.strain_presence_file
    sites_file = args.sites_file
    geno_file = args.geno_file
    output_seq_prefix = args.output_seq_prefix
    output_fasta = args.output_fasta
    output_map = args.output_map

    core_genes = set()
    with open(args.core_genes, 'r') as gene_fh:
        for gene_line in gene_fh:
            core_genes.add(gene_line.rstrip())

    with gzip.open(strain_presence_file, 'rt') as strain_presence_file_fh:
        strain_presence_columns = strain_presence_file_fh.readline().rstrip().split()

        col_i = strain_presence_columns.index(dataset)
        for strain_line in strain_presence_file_fh:
            strain_line = strain_line.rstrip()
            strain_split = strain_line.split()
            if strain_split[0] == species:
                if strain_split[col_i] == 'NA':
                    sys.exit('Stopping - no strains called as present')
                else:
                    present_strains = set(strain_split[col_i].split(','))

    sites = {}
    with gzip.open(sites_file, 'rt') as sites_file_fh:
        for sites_line in sites_file_fh:
            sites_line = sites_line.rstrip()
            sites_split = sites_line.split()

            site_index = int(sites_split[0])
            info = sites_split[1].split('|')
            sites[site_index] = {}

            sites[site_index]['gene'] = info[0]
            sites[site_index]['pos'] = int(info[1])
            sites[site_index]['ref'] = info[2]
            sites[site_index]['alt'] = info[3]

    gene_refs = defaultdict(str)
    current_name = None
    with gzip.open(ref_fasta_file, 'rt') as ref_fasta_file_fh:
        for ref_fasta_line in ref_fasta_file_fh:
            ref_fasta_line = ref_fasta_line.rstrip()
            if len(ref_fasta_line) == 0:
                continue
            if ref_fasta_line[0] == '>':
                seq_name = ref_fasta_line.split()[0][1:]
                if seq_name in core_genes:
                    current_sequence = True
                    current_name = seq_name
                else:
                    current_sequence = False
                    current_name = None
            elif current_sequence:
                gene_refs[current_name] += ref_fasta_line

    # Run length sanity check.
    for gene_seq in gene_refs.values():
        if len(gene_seq) % 3 != 0:
            sys.exit('Error - gene sequence length not divisible by three.')

    # Initialize strain sequences.
    inferred_seqs = {}
    for strain_i in present_strains:
        strain_identifier = output_seq_prefix + '_' + strain_i
        inferred_seqs[strain_identifier] = {}
        for core_gene_id, core_gene_seq in gene_refs.items():
            inferred_seqs[strain_identifier][core_gene_id] = core_gene_seq

    # Parse genotype output table. Note that the genotype is coded as a float, where
    # (as stated in the tutorial): '0.0 means entirely reference and 1.0 means entirely
    # alternative allele.'
    strainfacts_parsed_genes = set()

    with gzip.open(geno_file, 'rt') as geno_file_fh:
        geno_file_fh.readline()
        for geno_line in geno_file_fh:
            geno_line = geno_line.rstrip()
            geno_split = geno_line.split()

            # Strain not called as present, so ignore.
            if geno_split[0] not in present_strains:
                continue

            strain_id = output_seq_prefix + '_' + geno_split[0]

            seq_index = int(geno_split[1])
            geno_prob = float(geno_split[2])

            if geno_prob > 0.5:
                ref_geno = sites[seq_index]['ref']
                ref_length = len(ref_geno)

                # Ignore INDELs.
                if ref_length != len(sites[seq_index]['alt']):
                    continue

                gene = sites[seq_index]['gene']

                strainfacts_parsed_genes.add(gene)

                seq_slice1 = inferred_seqs[strain_id][gene][:sites[seq_index]['pos'] - 1]
                seq_slice2 = inferred_seqs[strain_id][gene][sites[seq_index]['pos'] + ref_length - 1:]

                exp_ref_seq = inferred_seqs[strain_id][gene][sites[seq_index]['pos'] - 1: sites[seq_index]['pos'] + ref_length - 1]

                if exp_ref_seq != ref_geno:
                    print('Error - expected reference genotype not found for this genotype line and site line:\n' + geno_line,
                          file=sys.stderr)
                    print('Expected ref: ' + exp_ref_seq, file=sys.stderr)
                    print('Observed ref: ' + ref_geno, file=sys.stderr)
                    print('Seq. position: ' + str(sites[seq_index]['pos']),
                          file=sys.stderr)
                    sys.exit()

                inferred_seqs[strain_id][gene] = seq_slice1 + sites[seq_index]['alt'] + seq_slice2

    # Loop through all core genes present in StrainFacts output in alphabetical order and produce a combined sequence for each strain.
    # Note that 102 bases are removed from the beginning and end of each core gene.
    # (This is to help correct for the problem of reads not mapping very well to the edges).
    full_core = {}
    sorted_genes = sorted(list(strainfacts_parsed_genes))
    out_map = open(output_map, 'w')
    individual_genes = defaultdict(dict)
    print('\t'.join(['gene', 'start', 'stop']), file = out_map)
    for strain_i in present_strains:
        strain_identifier = output_seq_prefix + '_' + strain_i
        full_core[strain_identifier] = ''
        for core_gene_id in sorted_genes:
            gene_start = len(full_core[strain_identifier]) + 1
            gene_stop = gene_start + len(inferred_seqs[strain_identifier][core_gene_id][102:-102]) - 1
            print('\t'.join([core_gene_id, str(gene_start), str(gene_stop)]), file = out_map)
            full_core[strain_identifier] += inferred_seqs[strain_identifier][core_gene_id][102:-102]
            individual_genes[core_gene_id][strain_identifier] = inferred_seqs[strain_identifier][core_gene_id][102:-102]
    out_map.close()

    # Final sanity check on core genomes.
    exp_length = len(full_core[strain_identifier])
    for strain_identifier in full_core.keys():
        if len(full_core[strain_identifier]) % 3 != 0:
            sys.exit('Error - final core genome not divisible by three.')

        if len(full_core[strain_identifier]) != exp_length:
            sys.exit('Error - not all final core genome sizes are equal.')

    # Write out core genome FASTA.
    with open(output_fasta, 'w') as out_fasta:
        for s in sorted(full_core.keys()):
            out_fasta.write('>' + dataset + '.' + s + '\n')
            out_fasta.write(textwrap.fill(full_core[s], width=70) + '\n')

    # Write out individual core gene FASTAs to separate files as well.
    os.makedirs(args.output_individual_fastas, exist_ok=True)
    for core_gene_id in sorted_genes:
        individual_gene_outfile = args.output_individual_fastas + '/' + core_gene_id + '.fna'
        with open(individual_gene_outfile, 'w') as individual_gene_out:
            for strain_id in sorted(individual_genes[core_gene_id].keys()):
                individual_gene_out.write('>' + dataset + '.' + strain_id + '\n')
                individual_gene_out.write(textwrap.fill(individual_genes[core_gene_id][strain_id], width=70) + '\n')

if __name__ == '__main__':
    main()
