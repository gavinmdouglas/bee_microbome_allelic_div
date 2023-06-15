#!/usr/bin/python3

import sys
from math import log
import itertools
import gzip


def pairwise_differences(seq1, seq2):

    if len(seq1) != len(seq2):
        sys.exit('Error - compared sequences differ in length.')

    bases = ['A', 'C', 'G', 'T']

    num_comparable_sites = 0
    num_diff = 0
    diff_positions = []

    for i in range(len(seq1)):
        base1 = seq1[i]
        base2 = seq2[i]
        if base1 in bases and base2 in bases:
            num_comparable_sites += 1

            if base1 != base2:
                num_diff += 1
                diff_positions.append(str(i))

    diff_positions_merged = ','.join(diff_positions)

    return((num_comparable_sites, num_diff, diff_positions_merged))


def pairwise_dnds(seq1, seq2):

    possible_sites_seq1 = exp_N_S_sites(seq1)
    possible_sites_seq2 = exp_N_S_sites(seq2)

    possible_sites = ((possible_sites_seq1[0] + possible_sites_seq2[0]) / 2,
                      (possible_sites_seq1[1] + possible_sites_seq2[1]) / 2)

    obs_subs = obs_N_S_subs(seq1, seq2)

    pn = obs_subs[0] / possible_sites[0]

    ps = obs_subs[1] / possible_sites[1]

    dn = sub_prop_to_rate(pn)

    ds = sub_prop_to_rate(ps)

    if ds > 0:
        dnds = dn / ds
    else:
        dnds = float('NaN')

    return((obs_subs[0], possible_sites[0], obs_subs[1], possible_sites[1], dn, ds, dnds))


def obs_N_S_subs(seq1, seq2):
    '''Number of non-synonymous and synonymous substitutions between two sequences (codon-aligned already).'''

    if len(seq1) != len(seq2):
        sys.exit("Stopping - lengths of two input sequences differ.")

    if len(seq1) % 3 != 0:
        sys.exit('Stopping - alignment length not a multiple of three: should be a codon alignment!')

    Nd = 0
    Sd = 0
    bases = ['A', 'C', 'G', 'T']

    for i in range(0, len(seq1) - 2, 3):

        codon1 = seq1[i:i + 3]
        codon2 = seq2[i:i + 3]

        if codon1 != codon2:

            # Skip codons if any non-standard bases.
            non_standard_base = False
            for codon_i in range(3):
                if codon1[codon_i] not in bases or codon2[codon_i] not in bases:
                    non_standard_base = True
                    break
            if non_standard_base:
                continue

            codon1_aa = codon_to_aa[codon1]
            codon2_aa = codon_to_aa[codon2]

            if codon1_aa == 'STOP' or codon2_aa == 'STOP':
                continue

            # Note that there can be multiple subs in same codon.
            # In these cases, the order of subs can change whether they are syn vs non-syn.
            # So, need to try every combination of sub orders and then compute mean numbers
            # of each sub type.
            mismatched_indices = [i for i in range(len(codon1)) if codon1[i] != codon2[i]]

            if len(mismatched_indices) > 1:
                sub_permutations = itertools.permutations(mismatched_indices)
            else:
                sub_permutations = [mismatched_indices]

            num_valid_permutations = 0
            num_potential_permutations = 0
            codon_Sd = 0
            codon_Nd = 0

            for sub_set in sub_permutations:

                num_potential_permutations += 1

                stop_in_path = False
                codon1_mutated = codon1
                tmp_Sd = 0
                tmp_Nd = 0

                for codon_i in sub_set:

                    codon2_base = codon2[codon_i]

                    codon1_mutated = codon1_mutated[:codon_i] + codon2_base + codon1_mutated[codon_i + 1:]

                    if codon_to_aa[codon1_mutated] == 'STOP':
                        stop_in_path = True
                        break

                    if codon_to_aa[codon1_mutated] == codon1_aa:
                        tmp_Sd += 1

                    else:
                        tmp_Nd += 1

                if not stop_in_path:
                    codon_Sd += tmp_Sd
                    codon_Nd += tmp_Nd
                    num_valid_permutations += 1

            # Get average numbers of Sd and Nd (assuming there are multiple possible permutations).
            if num_valid_permutations > 0:
                Sd += codon_Sd / num_valid_permutations
                Nd += codon_Nd / num_valid_permutations
            else:
                # Stop codons in all possible mutational paths.
                if num_potential_permutations == 2:
                    Sd += 0.5
                    Nd += 1.5
                elif num_potential_permutations == 6:
                    Sd += 1.0
                    Nd += 2.0
                else:
                    sys.exit('Stopping - the number of potential permutations is not adding up!')

    return((Nd, Sd))


def exp_N_S_sites(in_seq):
    '''Expected nonsynonymous and synonymous sites per sequence. Ignore codons
    that have any base besides four standard bases. Based on NG86 approach in KaKs_Calculator.
    This includes ignoring any possible stop codon paths and re-scaling the final counts to sum
    to the full sequence length.'''

    for lowercase_base in ['a', 'c', 't', 'g']:
        if lowercase_base in in_seq:
            sys.exit('Stopping - lower-case bases identified in input sequence. Make sure this is converted to all caps!')

    bases = ['A', 'C', 'G', 'T']

    S = 0
    N = 0

    if len(in_seq) % 3 != 0:
        sys.exit('Stopping - sequence length not a multiple of three: should be a codon alignment!')

    for i in range(0, len(in_seq) - 2, 3):

        syn = 0
        stop = 0

        codon = in_seq[i:i + 3]

        # Skip codon if any non-standard bases.
        non_standard_base = False
        for codon_i in range(3):
            if codon[codon_i] not in bases:
                non_standard_base = True
                break
        if non_standard_base:
            continue

        current_aa = codon_to_aa[codon]

        if current_aa == 'STOP':
            continue

        # Note that KaKs Calculator only looks at the first and last codon positions,
        # as syn subs can only occur there. However, I decided to check the second
        # position too, as subs causing stop codons can also occur there.
        for codon_i in range(3):

            codon_base = codon[codon_i]

            for base in bases:

                if base != codon_base:

                    new_codon = codon[:codon_i] + base + codon[codon_i + 1:]

                    if codon_to_aa[new_codon] == current_aa:
                        syn += 1

                    elif codon_to_aa[new_codon] == 'STOP':
                        stop += 1

        S += syn / 3
        N += (3 - (stop / 3) - (syn / 3))

    S /= 2
    N /= 2

    scaling_factor = len(in_seq) / (S + N)

    return((N * scaling_factor, S * scaling_factor))


def sub_prop_to_rate(p):

    if p == 0:
        return(0)

    renorm_val = (4 * p) / 3

    # Quick hack to allow wonky comparisons where the observed proportion of divergent codons is too high to compute a rate.
    if renorm_val >= 1:
        renorm_val = 0.99999

    return(-1 * (3/4) * log(1 - renorm_val))


def read_fasta(filename, cut_header=False, convert_upper=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:

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

    fasta_in.close()

    if convert_upper:
        for seq_name in seq.keys():
            seq[seq_name] = seq[seq_name].upper()

    return seq


codon_to_aa = {"TTT":"F", "TTC":"F",
                "TTA":"L", "TTG":"L",
                "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                "TAT":"Y", "TAC":"Y",
                "TAA":"STOP", "TAG":"STOP",
                "TGT":"C", "TGC":"C",
                "TGA":"STOP",
                "TGG":"W",
                "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                "CAT":"H", "CAC":"H",
                "CAA":"Q", "CAG":"Q",
                "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                "ATT":"I", "ATC":"I", "ATA":"I",
                "ATG":"M",
                "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                "AAT":"N", "AAC":"N",
                "AAA":"K", "AAG":"K",
                "AGT":"S", "AGC":"S",
                "AGA":"R", "AGG":"R",
                "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
                "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                "GAT":"D", "GAC":"D",
                "GAA":"E", "GAG":"E",
                "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
