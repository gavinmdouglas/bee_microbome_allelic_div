#!/usr/bin/python3

import argparse
import os
import sys
import pysam
import pandas as pd
import numpy as np
import random
from collections import defaultdict


def read_ids(filename):
    '''Read ids from a file into a list (one id per line)'''
    ids = list()
    with open(filename, 'r') as id_file:
        for line in id_file:
            ids.append(line.rstrip())

    return ids


class bed_coor:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


def main():

    parser = argparse.ArgumentParser(

            description='''
            Parse merged BCF to output depth info needed for StrainFacts (note that INDELs are ignored).
            Specfically used for producing the input for inferring overall strains,
            so a file for all core genes overall is produced rather than for individual genes.
            Also used for writing all sites with sufficient depth (including invariant sites),
            so that the core genome sequence can be reconstructed later.
            *Note that this version subsamples the reads at each site per sample to the specified number of reads. Any sites with fewer reads are ignored.*''',

            epilog='''Usage example:

        python bcfs_to_core_genome_strainfacts_input_subsampled.py --input INPUT
                                                                   --core_genes CORE_GENES
                                                                   --pos_samples SAMPLES_LIST
                                                                   --bed BED
                                                                   --outdir OUTDIR
                                                                   --subsample 50
                                                                   --random_seed 1

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='INPUT', type=str, required=True,
                        help='Path to input BCF of all samples.')

    parser.add_argument('-c', '--core_genes', metavar='CORE_GENES', type=str, required=True,
                        help='Path to text file with id of core genes to use - one per line.')

    parser.add_argument('-p', '--pos_samples', metavar='POS_SAMPLES', type=str, required=True,
                        help='Path to text file with id of samples where species was called as present - one per line. All other samples will be ignored.')

    parser.add_argument('-b', '--bed', metavar='BED', type=str, required=True,
                        help='Path to reference bed file with gene names and coordinates of interest.')

    parser.add_argument('-o', '--output', metavar='OUTPUT', type=str, required=True,
                        help='Path to directory for output files.')

    parser.add_argument('-m', '--min_prev', metavar='PREVALENCE', type=float, required=False, default = 0.05,
                        help='Minimum sample prevalence (as proportion) a site must have reads match alternative allele to be included.')

    parser.add_argument('-s', '--subsample', metavar='SUBSAMPLE_DEPTH', type=int, required=True,
                        help='Minimum depth for each site per sample to be included, and number of reads that will be subsampled per site/sample.')

    parser.add_argument('-r', '--random_seed', metavar='SEED', type=int, required=True,
                        help='Random seed to use (which will affect the subsampling of reads at each site per sample).')

    args = parser.parse_args()

    # Set random seed.
    np.random.seed(args.random_seed)

    min_depth = args.subsample

    bcf_in = pysam.VariantFile(args.input)

    gene_info = dict()

    pos_samples = read_ids(args.pos_samples)
    core_genes = read_ids(args.core_genes)

    core_genes_observed = 0

    with open(args.bed, 'r') as bed_filehandle:

        for bed_line in bed_filehandle:

            bed_line_split = bed_line.split()

            in_gene = bed_line_split[0]

            if in_gene in core_genes:
                core_genes_observed += 1

                gene_info[in_gene] = bed_coor(int(bed_line_split[1]),
                                              int(bed_line_split[2]))

    if core_genes_observed != len(core_genes):
        sys.exit('Error - number of observed genes does not match expected.')

    # Get sample ids, assuming that they are the first field after splitting by '.'
    # (and removing and preceding path if present).
    raw_samples = list(bcf_in.header.samples)
    samples = []
    pos_samples_observed = 0
    for raw_sample in raw_samples:
        sample_basename = os.path.basename(raw_sample)
        sample_id = sample_basename.split('.')[0]
        samples.append(sample_id)

        if sample_id in pos_samples:
            pos_samples_observed += 1

    num_all_samples = len(samples)

    if pos_samples_observed != len(pos_samples):
        sys.exit('Error - number of observed samples where species is present does not match expected.')

    site_info = []

    sample_alt_depth = defaultdict(list)
    sample_ref_depth = defaultdict(list)
    sample_total_depth = defaultdict(list)

    # Also keep track of depth at invariant sites, so these
    # can be used to reconstruct the final strain core genome
    # sequences later.
    invariant_site_info = []
    sample_invariant_total_depth = defaultdict(list)

    for gene in core_genes:

        for rec in bcf_in.fetch(contig = gene,
                                start = gene_info[gene].start,
                                stop = gene_info[gene].stop):

            rec_split = str(rec).split()

            site_pos = rec_split[1]
            ref_allele = rec_split[3]
            alt_allele = rec_split[4]

            # Skip INDELs.
            if len(ref_allele) > 1 or len(alt_allele) > 1:
                continue

            # Skip multi-allelic sites.
            if ',' in alt_allele:
                continue

            sample_genotypes = rec_split[-num_all_samples:]
            format_fields = rec_split[len(rec_split) - num_all_samples - 1].split(':')
            AD_index = format_fields.index('AD')

            # Keep track of invariant sites (but treat separately).
            if alt_allele == '.':
                invariant_site_info.append('|'.join([gene, site_pos]))

                for idx, geno_info in enumerate(sample_genotypes):

                    if samples[idx] not in pos_samples:
                        continue

                    invariant_depth_raw = geno_info.split(':')[AD_index]
                    invariant_ref_depth = int(invariant_depth_raw.replace('.', '0'))

                    if invariant_ref_depth > min_depth:
                        invariant_ref_depth = min_depth

                    sample_invariant_total_depth[idx].append(invariant_ref_depth)

            else:
                # Keep track of biallelic sites separately, as these are the main focus.
                site_info.append('|'.join([gene, site_pos, ref_allele, alt_allele]))

                for idx, geno_info in enumerate(sample_genotypes):

                    if samples[idx] not in pos_samples:
                        continue

                    allele_depth_raw = geno_info.split(':')[AD_index]
                    allele_depth_raw = allele_depth_raw.replace('.', '0')
                    allele_depth = allele_depth_raw.split(',')

                    ref_depth = int(allele_depth[0])

                    if len(allele_depth) == 2:
                        alt_depth = int(allele_depth[1])
                    elif len(allele_depth) == 1:
                        alt_depth = 0
                    else:
                        print(geno_info, file = sys.stderr)
                        sys.exit('Problem with allele depth?')

                    site_total_dp = ref_depth + alt_depth

                    # If the site depth >= the minimum depth, then subsample the reads to this cut-off.
                    if site_total_dp >= args.subsample:
                        if alt_depth == 0:
                            ref_depth = args.subsample
                        elif ref_depth == 0:
                            alt_depth = args.subsample
                        else:
                            allele_population = np.array(['ref'] * ref_depth + ['alt'] * alt_depth)
                            sampled_indices = np.random.choice(len(allele_population), args.subsample, replace=False)
                            sampled_elements = allele_population[sampled_indices]
                            ref_depth = np.sum(sampled_elements == 'ref')
                            alt_depth = np.sum(sampled_elements == 'alt')

                        site_total_dp = ref_depth + alt_depth

                    sample_ref_depth[idx].append(ref_depth)
                    sample_alt_depth[idx].append(alt_depth)
                    sample_total_depth[idx].append(site_total_dp)

    total_depth_df = pd.DataFrame.from_dict(sample_total_depth)

    # When parsing out core genes - filter by site first, since these genes are expected to be found in these samples already.
    # I.e., the limiting data here is # samples rather than # sites. The # sites is definitely more of a problem for the
    # individual genes rather than overall strain.
    site_nonzero_sites_prop = total_depth_df[total_depth_df >= min_depth].count(axis = 1) / total_depth_df.shape[1]
    total_depth_df_filt = total_depth_df.loc[site_nonzero_sites_prop > 0.9, :]

    # Filter by minimum prevalence of samples with at least 1 alt read, if specified.
    if args.min_prev > 0:
        alt_depth_df = pd.DataFrame.from_dict(sample_alt_depth)
        alt_depth_df_filt = alt_depth_df.loc[site_nonzero_sites_prop > 0.9, :]
        site_alt_sites_prop_filt = alt_depth_df_filt[alt_depth_df_filt > 0].count(axis = 1) / alt_depth_df_filt.shape[1]
        total_depth_df_filt = total_depth_df_filt.loc[site_alt_sites_prop_filt >= args.min_prev, :]

    sample_nonzero_sites_prop = total_depth_df_filt[total_depth_df_filt >= min_depth].count(axis = 0) / total_depth_df_filt.shape[0]
    total_depth_df_filt = total_depth_df_filt.loc[:, sample_nonzero_sites_prop > 0.9]

    sample_subset = [samples[i] for i in list(total_depth_df_filt.columns)]
    site_info_subset = [site_info[i] for i in list(total_depth_df_filt.index)]

    # If this approach resulted in too few datapoints, then try doing the reverse and see if that makes a difference.
    if len(sample_subset) < 2 or len(site_info_subset) < 10:

        print('Filtering by sites first resulted in too many datapoints being thrown out. Will try filtering out samples first instead.', file = sys.stderr)

        sample_nonzero_sites_prop = total_depth_df[total_depth_df >= min_depth].count(axis = 0) / total_depth_df.shape[0]
        total_depth_df_filt = total_depth_df.loc[:, sample_nonzero_sites_prop > 0.9]

        site_nonzero_sites_prop = total_depth_df_filt[total_depth_df_filt >= min_depth].count(axis = 1) / total_depth_df_filt.shape[1]
        total_depth_df_filt = total_depth_df_filt.loc[site_nonzero_sites_prop > 0.9, :]

        # Filter by minimum prevalence of samples with at least 1 alt read, if specified.
        if args.min_prev > 0:
            alt_depth_df_filt = alt_depth_df.loc[site_nonzero_sites_prop > 0.9, :]
            site_alt_sites_prop_filt = alt_depth_df_filt[alt_depth_df_filt > 0].count(axis = 1) / alt_depth_df_filt.shape[1]
            total_depth_df_filt = total_depth_df_filt.loc[site_alt_sites_prop_filt >= args.min_prev, :]

        sample_subset = [samples[i] for i in list(total_depth_df_filt.columns)]
        site_info_subset = [site_info[i] for i in list(total_depth_df_filt.index)]

        if len(sample_subset) < 2 or len(site_info_subset) < 10:
            sys.exit('Stopping as at least the number of samples or number of sites are below very lenient cut-offs of min of 2 samples and 10 sites.')

    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    else:
        print('Stopping - output folder already exists.', file = sys.stderr)
        sys.exit()

    ref_depth_df = pd.DataFrame.from_dict(sample_ref_depth)
    alt_depth_df = pd.DataFrame.from_dict(sample_alt_depth)

    sample_ref_depth_filt = ref_depth_df.loc[total_depth_df_filt.index, total_depth_df_filt.columns]
    sample_alt_depth_filt = alt_depth_df.loc[total_depth_df_filt.index, total_depth_df_filt.columns]

    sample_ref_depth_filt.reset_index(inplace = True, drop = True)
    sample_ref_depth_filt.columns = range(sample_ref_depth_filt.shape[1])

    sample_ref_depth_filt_melt = pd.melt(sample_ref_depth_filt, ignore_index = False,
                                             value_name = 'metagenotype', var_name = 'sample')
    sample_ref_depth_filt_melt['position'] = sample_ref_depth_filt_melt.index
    sample_ref_depth_filt_melt['allele'] = 'ref'
    sample_ref_depth_filt_melt.reset_index(inplace = True, drop = True)

    sample_alt_depth_filt.reset_index(inplace = True, drop = True)
    sample_alt_depth_filt.columns = range(sample_alt_depth_filt.shape[1])

    sample_alt_depth_filt_melt = pd.melt(sample_alt_depth_filt, ignore_index = False,
                                             value_name = 'metagenotype', var_name = 'sample')
    sample_alt_depth_filt_melt['position'] = sample_alt_depth_filt_melt.index
    sample_alt_depth_filt_melt['allele'] = 'alt'
    sample_alt_depth_filt_melt.reset_index(inplace = True, drop = True)

    combined_gene_data = pd.concat([sample_ref_depth_filt_melt, sample_alt_depth_filt_melt],
                                    ignore_index = True, sort = False)

    combined_gene_data = combined_gene_data[['sample', 'position', 'allele', 'metagenotype']]

    outfile = args.output + '/metagenotypes.tsv'
    combined_gene_data.to_csv(path_or_buf = outfile, sep = '\t', index = False)

    sample_outfile = args.output + '/samples.tsv'
    pd.Series(sample_subset).to_csv(path_or_buf = sample_outfile, sep = '\t', index = True, header = False)

    site_outfile = args.output + '/sites.tsv'
    pd.Series(site_info_subset).to_csv(path_or_buf = site_outfile, sep = '\t', index = True, header = False)

    # Finally, also determine invariant sites that pass the minimum depth cut-off (for the subset of samples selected).
    total_invariant_depth_df = pd.DataFrame.from_dict(sample_invariant_total_depth)
    total_invariant_depth_df = total_invariant_depth_df.loc[:, total_depth_df_filt.columns]
    invariant_depth_passing_prop = total_invariant_depth_df[total_invariant_depth_df >= min_depth].count(axis = 1) / total_depth_df_filt.shape[1]

    total_invariant_depth_df_filt = total_invariant_depth_df.loc[invariant_depth_passing_prop > 0.9, :]
    passable_invariant_sites = [invariant_site_info[i] for i in list(total_invariant_depth_df_filt.index)]

    passable_invariant_site_outfile = args.output + '/passable_invariant_sites.tsv'
    with open(passable_invariant_site_outfile, 'w') as invariant_out_fh:
        for passable_invariant_site in passable_invariant_sites:
            passable_invariant_site = passable_invariant_site.split('|')
            print('\t'.join(passable_invariant_site), file = invariant_out_fh)


if __name__ == '__main__':
    main()
