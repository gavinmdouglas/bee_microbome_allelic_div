#!/usr/bin/python3

# Read through CD-HIT output file and print out mapping of
# all clusters to *all* member genes. One link between
# cluster and member per-line.


def main():

    print('\t'.join(['cluster', 'gene']))

    with open('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/pangenome_c0.95_aL0.90_s0.90.clstr', 'r') as cluster_out:
        for cluster_line in cluster_out:
            if cluster_line[0] == '>':
                current_cluster = cluster_line[1:].rstrip()
                continue

            cluster_line_split = cluster_line.split()

            gene = cluster_line_split[2][1:-3]
            print(current_cluster + '\t' + gene)


if __name__ == '__main__':
    main()
