#!/usr/bin/python3

import argparse
import sys
import textwrap
import gzip


def main():

    parser = argparse.ArgumentParser(

            description='''

            Parse StrainFacts output files, and supplementary files, to produce FASTA of all inferred sequences.

            Note that this script is designed to be run on a single gene at a time.

            Will create two output FASTAs: one containing the full inferred genes (untrimmed), and one of only the sites that had sufficient depth to be called.

            Note that INDELs will be ignored.
            ''',

            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-n', '--gene_name',
                        metavar='GENE_NAME', type=str,
                        help="Gene name to parse.",
                        required=True)

    parser.add_argument('-r', '--ref_fasta_file',
                        metavar='REF_FASTA', type=str,
                        help="Path to gzipped FASTA file containing all representative gene sequences",
                        required=True)

    parser.add_argument('-a', '--allele_presence_file',
                        metavar='ALLELE_PRESENCE', type=str,
                        help="Path to gzipped table with mapping of gene IDs to allele numbers called as present across samples.",
                        required=True)

    parser.add_argument('-s', '--sites_file',
                        metavar='SITES_FILE', type=str,
                        help="Path gzipped StrainFacts-output sites file.",
                        required=True)

    parser.add_argument('-g', '--geno_file',
                        metavar='GENO_FILE', type=str,
                        help="Path gzipped StrainFacts-output genotype file.",
                        required=True)

    parser.add_argument('-i', '--invariant_file',
                        metavar='INVARIANT_FILE', type=str,
                        help="Path to gzipped invariant sites file created when prepping StrainFacts input.",
                        required=True)

    parser.add_argument('-p', '--output_seq_prefix',
                        metavar='SEQNAME_PREFIX', type=str,
                        help="String to affix to start of output sequence number (delimited by '_')",
                        required=False, default='allele')

    parser.add_argument('-o1', '--output_file_full',
                        metavar='OUTPUT_FILE', type=str,
                        help="Path to output file.",
                        required=True)

    parser.add_argument('-o2', '--output_file_informative',
                        metavar='OUTPUT_FILE', type=str,
                        help="Path to output file.",
                        required=True)

    args = parser.parse_args()

    gene_name = args.gene_name
    ref_fasta_file = args.ref_fasta_file
    allele_presence_file = args.allele_presence_file
    sites_file = args.sites_file
    geno_file = args.geno_file
    invariant_file = args.invariant_file
    output_seq_prefix = args.output_seq_prefix
    output_file_full = args.output_file_full
    output_file_informative = args.output_file_informative

    with gzip.open(allele_presence_file, 'rt') as allele_presence_file_fh:
        allele_presence_file_fh.readline()
        for allele_line in allele_presence_file_fh:
            allele_line = allele_line.rstrip()
            allele_split = allele_line.split()
            if allele_split[0] == gene_name:
                present_alleles = set(allele_split[1].split(','))

    gene_ref = ''
    with gzip.open(ref_fasta_file, 'rt') as ref_fasta_file_fh:
        for ref_fasta_line in ref_fasta_file_fh:
            ref_fasta_line = ref_fasta_line.rstrip()
            if len(ref_fasta_line) == 0:
                continue
            if ref_fasta_line[0] == '>':
                seq_name = ref_fasta_line.split()[0][1:]
                if seq_name == gene_name:
                    current_sequence = True
                else:
                    current_sequence = False
            elif current_sequence:
                gene_ref += ref_fasta_line

    sites = {}
    with gzip.open(sites_file, 'rt') as sites_file_fh:
        for sites_line in sites_file_fh:
            sites_line = sites_line.rstrip()
            sites_split = sites_line.split()

            site_index = int(sites_split[0])
            info = sites_split[1].split('_')
            sites[site_index] = {}
            sites[site_index]['pos'] = int(info[0])
            sites[site_index]['ref'] = info[1]
            sites[site_index]['alt'] = info[2]

    # Keep track of invariant sites as informative, which
    # the variable sites (with sufficient coverage)
    # will be added to.
    invariant_sites = set()
    with gzip.open(invariant_file, 'rt') as invariant_file_fh:
        for invariant_line in invariant_file_fh:
            invariant_line = invariant_line.rstrip()
            invariant_split = invariant_line.split()

            # Sanity check:
            if int(invariant_split[1]) < 101:
                print()
                print('Error - invariant site position is less than 101, which should have been trimmed off. In this file: ' + invariant_file,
                      file=sys.stderr)
                sys.exit()
            invariant_sites.add(int(invariant_split[1]))

    # Parse genotype output table. Note that the genotype is coded as a float, where
    # (as stated in the tutorial): '0.0 means entirely reference and 1.0 means entirely
    # alternative allele.'
    full_untrimmed_seqs = {}

    informative_sites = invariant_sites.copy()

    with gzip.open(geno_file, 'rt') as geno_file_fh:
        geno_file_fh.readline()
        for geno_line in geno_file_fh:
            geno_line = geno_line.rstrip()
            geno_split = geno_line.split()

            # Strain not called as present, so ignore.
            if geno_split[0] not in present_alleles:
                continue

            seq_index = int(geno_split[1])
            geno_prob = float(geno_split[2])

            ref_geno = sites[seq_index]['ref']
            ref_length = len(ref_geno)

            # Ignore INDELs.
            if ref_length != len(sites[seq_index]['alt']):
                continue

            strain_id = output_seq_prefix + '_' + geno_split[0]

            if strain_id not in full_untrimmed_seqs.keys():
                full_untrimmed_seqs[strain_id] = gene_ref

            # If site is informative, add to informative_sites.
            seq_pos = sites[seq_index]['pos']
            if seq_pos in invariant_sites:
                print(geno_file)
                print(seq_index)
                print(seq_pos)
                print('Error - variant site is already in informative_sites, but is being added again.',
                      file=sys.stderr)
                sys.exit()
            informative_sites.add(seq_pos)

            if geno_prob > 0.5:

                seq_slice1 = full_untrimmed_seqs[strain_id][:seq_pos - 1]
                seq_slice2 = full_untrimmed_seqs[strain_id][seq_pos + ref_length - 1:]

                exp_ref_seq = full_untrimmed_seqs[strain_id][seq_pos - 1: seq_pos + ref_length - 1]

                if exp_ref_seq != ref_geno:
                    print('Error - expected reference genotype not found for this genotype line and site line:\n' + geno_line,
                          file=sys.stderr)
                    print('Expected ref: ' + exp_ref_seq, file=sys.stderr)
                    print('Observed ref: ' + ref_geno, file=sys.stderr)
                    print('Seq. position: ' + str(seq_pos),
                          file=sys.stderr)
                    sys.exit()

                full_untrimmed_seqs[strain_id] = seq_slice1 + sites[seq_index]['alt'] + seq_slice2

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).

    # Also get FASTA of informative sites.
    informative_fasta = {}
    informative_sites = sorted(list(informative_sites))

    with open(output_file_full, "w") as out_fasta:
        for s in sorted(full_untrimmed_seqs.keys()):
            informative_fasta[s] = ''
            for i in informative_sites:
                informative_fasta[s] += full_untrimmed_seqs[s][i - 1]
            out_fasta.write(">" + s + "\n")
            out_fasta.write(textwrap.fill(full_untrimmed_seqs[s], width=70) + "\n")

    # Then write out informative FASTA.
    with open(output_file_informative, "w") as out_informative:
        for s in sorted(informative_fasta.keys()):
            out_informative.write(">" + s + "\n")
            out_informative.write(textwrap.fill(informative_fasta[s], width=70) + "\n")


if __name__ == '__main__':
    main()
