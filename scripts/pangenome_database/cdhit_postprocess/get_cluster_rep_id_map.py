#!/usr/bin/python3

import sys
from collections import defaultdict


def main():

    print('\t'.join(['cluster', 'gene']))

    with open('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/cdhit_out/pangenome_c0.95_aL0.90_s0.90.clstr', 'r') as cluster_out:
        for cluster_line in cluster_out:
            if cluster_line[0] == '>':
                current_cluster = cluster_line[1:].rstrip()
                continue
            cluster_line_split = cluster_line.split()

            if cluster_line_split[-1] == "*":
                gene = cluster_line_split[2][1:-3]
                print(current_cluster + '\t' + gene)


if __name__ == '__main__':
    main()
