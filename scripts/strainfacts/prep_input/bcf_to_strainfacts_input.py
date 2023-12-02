#!/usr/bin/python3

import argparse
import os
import sys
import pysam
import pandas as pd
from collections import defaultdict


class bed_coor:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


def main():

    parser = argparse.ArgumentParser(

            description='''
            Parse merged BCF to output depth info needed for StrainFacts.
            Assumes that each separate contig is a separate gene that should be input separately into StrainFacts.
            Will also output invariant sites that nonetheless have sufficient depth, for later downstream use.
            Note that INDELs are ignored.
            ''',

            epilog='''Usage example:

    python bcfs_to_strainfacts_input.py --input INPUT --bed BED --outdir OUTDIR

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='INPUT', type=str, required=True,
                        help='Path to input BCF of all samples.')

    parser.add_argument('-b', '--bed', metavar='INPUT', type=str, required=True,
                        help='Path to reference bed file with gene names and coordinates of interest.')

    parser.add_argument('-o', '--output', metavar='OUTPUT', type=str,
                        help='Path to directory for output files.', required=True)

    parser.add_argument('--min_sites', metavar='MIN_SITES', type=int,
                        help='Min number of sites. Note that I believe sfacts cannot work without at least two sites', required=False, default = 10)

    parser.add_argument('--min_samples', metavar='MIN_SAMPLES', type=int,
                        help='Min number of samples. Note that I believe sfacts cannot work without at least two samples', required=False, default = 2)

    parser.add_argument('-m', '--min_prev', metavar='PREVALENCE', type=float, required=False, default = 0.05,
                        help='Minimum sample prevalence (as proportion) a site must have reads match alternative allele to be included.')

    args = parser.parse_args()

    bcf_in = pysam.VariantFile(args.input)

    gene_info = dict()

    with open(args.bed, 'r') as bed_filehandle:

        for bed_line in bed_filehandle:

            bed_line_split = bed_line.split()

            in_gene = bed_line_split[0]

            gene_info[in_gene] = bed_coor(int(bed_line_split[1]),
                                          int(bed_line_split[2]))

    if not os.path.isdir(args.output):
        os.makedirs(args.output)
        os.makedirs(args.output + '/metagenotypes')
        os.makedirs(args.output + '/samples')
        os.makedirs(args.output + '/sites')
        os.makedirs(args.output + '/passable_invariant_sites')
    else:
        print('Stopping - output folder already exists.', file = sys.stderr)
        sys.exit()

    for gene in sorted(gene_info.keys()):

        # Get sample ids, assuming that they are the first field after splitting by '.'
        # (and removing and preceding path if present).
        raw_samples = list(bcf_in.header.samples)
        samples = []
        for raw_sample in raw_samples:
            sample_basename = os.path.basename(raw_sample)
            samples.append(sample_basename.split('.')[0])

        num_samples = len(samples)

        site_info = []

        sample_alt_depth = defaultdict(list)
        sample_ref_depth = defaultdict(list)
        sample_total_depth = defaultdict(list)

        # Also keep track of depth at invariant sites, so these
        # can be used to reconstruct the final strain core genome
        # sequences later.
        invariant_site_info = []
        sample_invariant_total_depth = defaultdict(list)

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

            sample_genotypes = rec_split[-num_samples:]
            format_fields = rec_split[len(rec_split) - num_samples - 1].split(':')
            AD_index = format_fields.index('AD')

            # Keep track of invariant sites (but treat separately).
            if alt_allele == '.':
                invariant_site_info.append('|'.join([gene, site_pos]))

                for idx, geno_info in enumerate(sample_genotypes):
                    invariant_depth_raw = geno_info.split(':')[AD_index]
                    invariant_ref_depth = int(invariant_depth_raw.replace('.', '0'))
                    sample_invariant_total_depth[idx].append(invariant_ref_depth)

            else:
                # Keep track of biallelic sites separately, as these are the main focus.
                site_info.append('_'.join([site_pos, ref_allele, alt_allele]))

                sample_genotypes = rec_split[-num_samples:]

                format_fields = rec_split[len(rec_split) - num_samples - 1].split(':')
                AD_index = format_fields.index('AD')

                for idx, geno_info in enumerate(sample_genotypes):

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

                    sample_ref_depth[idx].append(ref_depth)
                    sample_alt_depth[idx].append(alt_depth)
                    sample_total_depth[idx].append(ref_depth + alt_depth)

        total_depth_df = pd.DataFrame.from_dict(sample_total_depth)

        # First exclude all samples without at least 10 reads at > 90% of sites for this gene.
        sample_nonzero_sites_prop = total_depth_df[total_depth_df >= 10.0].count(axis = 0) / total_depth_df.shape[0]
        total_depth_df_filt = total_depth_df.loc[:, sample_nonzero_sites_prop > 0.9]

        # Then exclude all sites that don't have at least 10 reads across > 90% of these samples.
        site_nonzero_sites_prop = total_depth_df_filt[total_depth_df_filt >= 10.0].count(axis = 1) / total_depth_df_filt.shape[1]
        total_depth_df_filt = total_depth_df_filt.loc[site_nonzero_sites_prop > 0.9, :]

        # Filter by minimum prevalence of samples with at least 1 alt read, if specified.
        if args.min_prev > 0:
            alt_depth_df = pd.DataFrame.from_dict(sample_alt_depth)
            alt_depth_df_filt = alt_depth_df.loc[site_nonzero_sites_prop > 0.9, :]
            site_alt_sites_prop_filt = alt_depth_df_filt[alt_depth_df_filt > 0].count(axis = 1) / alt_depth_df_filt.shape[1]
            total_depth_df_filt = total_depth_df_filt.loc[site_alt_sites_prop_filt >= args.min_prev, :]

        sample_subset = [samples[i] for i in list(total_depth_df_filt.columns)]
        site_info_subset = [site_info[i] for i in list(total_depth_df_filt.index)]

        if len(sample_subset) < args.min_samples or len(site_info_subset) < args.min_sites:
            print('Skipping gene ' + gene + ' as at least the number of samples or number of sites are below the set cut-offs.', file = sys.stderr)
            continue

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

        outfile = args.output + '/metagenotypes/' + gene + '_metagenotype.tsv'
        combined_gene_data.to_csv(path_or_buf = outfile, sep = '\t', index = False)

        sample_outfile = args.output + '/samples/' + gene + '_samples.tsv'
        pd.Series(sample_subset).to_csv(path_or_buf = sample_outfile, sep = '\t', index = True, header = False)

        site_outfile = args.output + '/sites/' + gene + '_sites.tsv'
        pd.Series(site_info_subset).to_csv(path_or_buf = site_outfile, sep = '\t', index = True, header = False)

        # Finally, also determine invariant sites that pass the minimum depth cut-off (for the subset of samples selected).
        total_invariant_depth_df = pd.DataFrame.from_dict(sample_invariant_total_depth)

        # Return empty list if dataframe if dimensions are 0.
        # (This is the case if there are no invariant sites.)
        if total_invariant_depth_df.shape[0] == 0:
            passable_invariant_sites = []
        else: 
            total_invariant_depth_df = total_invariant_depth_df.loc[:, total_depth_df_filt.columns]
            invariant_depth_passing_prop = total_invariant_depth_df[total_invariant_depth_df >= 10.0].count(axis = 1) / total_depth_df_filt.shape[1]

            total_invariant_depth_df_filt = total_invariant_depth_df.loc[invariant_depth_passing_prop > 0.9, :]
            passable_invariant_sites = [invariant_site_info[i] for i in list(total_invariant_depth_df_filt.index)]

        passable_invariant_site_outfile = args.output + '/passable_invariant_sites/' + gene + '.tsv'
        with open(passable_invariant_site_outfile, 'w') as invariant_out_fh:
            for passable_invariant_site in passable_invariant_sites:
                passable_invariant_site = passable_invariant_site.split('|')
                print('\t'.join(passable_invariant_site), file = invariant_out_fh)


if __name__ == '__main__':
    main()
