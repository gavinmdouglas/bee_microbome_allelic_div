#!/usr/bin/python3

import sys
import os
from collections import defaultdict
import textwrap

# Parse core gene alignments per species and get merged alignments per species (split by each dataset separately).
# Replace all non-passable sites with N's.

def read_fasta(filename, cut_header=False, convert_upper=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    def parse_fasta_lines(file_handle):
        for line in file_handle:

            line = line.rstrip()

            if len(line) == 0:
                continue

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                name = name.rstrip("\r\n")

                # Make sure that sequence id is not already in dictionary.
                if name in seq:
                    sys.stderr("Stopping due to duplicated id in file: " + name)

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        with gzip.open(filename, "rt") as fasta_in:
            parse_fasta_lines(fasta_in)
    else:
        with open(filename, "r") as fasta_in:
            parse_fasta_lines(fasta_in)

    if convert_upper:
        for seq_name in seq.keys():
            seq[seq_name] = seq[seq_name].upper()

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).
    for s in sorted(seq.keys()):
        out_fasta.write(">" + s + "\n")
        out_fasta.write(textwrap.fill(seq[s], width=70) + "\n")

    out_fasta.close()

all_species = []
with open('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/species_w_strains.txt', 'r') as species_fh:
    for species_line in species_fh:
        all_species.append(species_line.rstrip())

for species in all_species:

    # Get path to all files ending in '.fna' in alignment_dir.
    alignment_dir = '/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/trimmed_individual_core_gene_fastas_w_ref_align/' + species
    fasta_paths = [alignment_dir + '/' + f for f in os.listdir(alignment_dir) if f.endswith('.fna')]

    gene_alignments = dict()
    for fasta_path in fasta_paths:
        gene_id = os.path.basename(fasta_path).split('.')[0]
        gene_alignments[gene_id] = read_fasta(fasta_path)
    
    # Loop through individual datasets and merge alignments (excluding sequences for all other datasets).
    # Also, restrict to sites subset as comparable within this dataset specifically (and do not consider all other sites).
    sp_datasets = []
    with open('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/' + species + '.datasets.txt', 'r') as datasets_fh:
        for dataset_line in datasets_fh:
            sp_datasets.append(dataset_line.rstrip())

    for dataset in sp_datasets:
        dataset_sites = defaultdict(set)
        dataset_genes = set()

        with open('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/' + species + '.' + dataset + '.sites.txt', 'r') as sites_fh:
            for site_line in sites_fh:
                site_line = site_line.rstrip().split('\t')
                # Assuming that 102 bases were trimmed from the start and end of the genes.
                dataset_sites[site_line[0]].add(int(site_line[1]) - 103)
                dataset_genes.add(site_line[0])

        dataset_merged_alignment = defaultdict(str)
        first_gene = True
        for gf in sorted(list(dataset_genes)):
            raw_seqs = gene_alignments[gf]

            uniq_seq_ids = list(raw_seqs.keys())
            seq_length = len(raw_seqs[uniq_seq_ids[0]])
            for uniq_seq_id in uniq_seq_ids[1:]:
                if seq_length != len(raw_seqs[uniq_seq_id]):
                    sys.exit('Error: sequences for this gene are of different length, but they should be aligned: ' + gf)

            # Read through and run sanity check that all expected genes and genomes are represented.
            # First read through sequences is used to identify all indices that should be excluded.
            represented_genomes = set()
            gene_past_genomes = set()
            indices_to_ignore = set()

            raw_to_genome_id = dict()

            for raw_seq_id, raw_seq in raw_seqs.items():

                if '||' in raw_seq_id:
                    # Ref. genome - parse out gene ID specifically.
                    genome_id = raw_seq_id.split('||')[0]
                else:
                    # Non-ref. genome - use whole ID (as long as it's the right dataset - skip otherwise).
                    if dataset not in raw_seq_id:
                        continue

                    genome_id = raw_seq_id

                raw_to_genome_id[raw_seq_id] = genome_id

                if genome_id in gene_past_genomes:
                    sys.exit('Error: duplicate genome ID found in this gene family: ' + gf)
                else:
                    gene_past_genomes.add(genome_id)

                raw_seq = raw_seq.upper()

                dna_i = 0
                for seq_i in range(len(raw_seq)):

                    # Skip gaps.
                    if raw_seq[seq_i] == '.' or raw_seq[seq_i] == '-':
                        continue

                    if dna_i not in dataset_sites[gf]:
                        indices_to_ignore.add(seq_i)
                    
                    dna_i += 1

                if genome_id not in represented_genomes:
                    represented_genomes.add(genome_id)
                else:
                    sys.exit('Error: multiple gene ids matched to same genome for this gene family: ' + gf)

            # Loop through again and add sequences to merged alignment.
            for raw_seq_id, raw_seq in raw_seqs.items():
                if raw_seq_id not in raw_to_genome_id.keys():
                    continue
                genome_id = raw_to_genome_id[raw_seq_id]
                raw_seq = raw_seq.upper()
                for seq_i in range(len(raw_seq)):
                    if seq_i not in indices_to_ignore:
                        dataset_merged_alignment[genome_id] += raw_seq[seq_i]

            if first_gene:
                genome_ids = sorted(list(dataset_merged_alignment.keys()))       
                first_gene = False

            else:
                # Sanity checks on full set of genes and genomes for this gene family:
                if genome_ids != sorted(list(represented_genomes)):
                    print(genome_ids)
                    print(sorted(list(represented_genomes)))
                    sys.exit('Mismatch in observed and expected genomes for this gene family: ' + gf)

        # Final sanity check that all sequences are the same length.
        seq_length = len(dataset_merged_alignment[genome_ids[0]])
        for genome_id in genome_ids[1:]:
            if seq_length != len(dataset_merged_alignment[genome_id]):
                sys.exit('Error: not all final core genome alignments are the same length for species ' + species + ' and dataset ' + dataset)
        
        outfile = '/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/merged_aligned_core_genes_strainfacts_and_ref/' + species + '.' + dataset + '.fna'

        write_fasta(dataset_merged_alignment, outfile)
