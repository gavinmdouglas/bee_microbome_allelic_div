#!/usr/bin/python3

import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser(

            description="Trim bedfile coordinates by a specified number of "
                        "base-pairs on each end. Any features below the min "
                        "size will not be printed out. Note that this script "
                        "was written expecting a bedfile of gene coordinates "
                        "to be input. Note that this script will only consider "
                        "the first three fields.",

    epilog='''Usage example:

    python trim_bed_coor.py --input BED --trim_size 100 --min_final_size 100 --output BED

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FASTA", type=str,
                        help="Path to input BED (e.g., gene coordinates)",
                        required=True)

    parser.add_argument("-o", "--output", metavar="BED", type=str,
                        help="Path to output BED", required=True)

    parser.add_argument("-t", "--trim_size", metavar="INT", type=int,
                        help="Number of bases to trim from the sides of each feature.",
                        required=True)

    parser.add_argument("-m", "--min_final_size", metavar="INT", type=int,
                        help="Minimum size of features (after trimming) to retain for output.",
                        required=True)

    args = parser.parse_args()

    total_features = 0
    num_features_below_cutoff = 0

    with open(args.output, "w") as OUTPUT_BED:

        with open(args.input, "r") as INPUT_BED:

            for line in INPUT_BED:
                line = line.rstrip()
                line_split = line.split()

                scaffold = line_split[0]
                new_start = int(line_split[1]) + args.trim_size
                new_end = int(line_split[2]) - args.trim_size

                new_size = new_end - new_start

                total_features += 1
                
                if new_size < args.min_final_size:
                    num_features_below_cutoff += 1
                
                else:
                    print(scaffold + "\t" + str(new_start) + "\t" +
                          str(new_end), file = OUTPUT_BED)

    print(str(num_features_below_cutoff) + " of " + str(total_features) +
          " trimmed features were below the min size cut-off and were not output.")

if __name__ == '__main__':
    main()
