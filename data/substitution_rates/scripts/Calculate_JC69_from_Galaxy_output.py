#########################################
# Tim Webster, Arizona State University #
# July 11, 2018                         #
# Edited by Daniel Cotter 7/12/18       #
#########################################

from __future__ import division
from __future__ import print_function
import argparse
import csv
import math


def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(
        description="This script adds JC69 substitution-corrected divergence "
        "as an additional column to the output of Galaxy's Estimate "
        "Substitution Rates on windows.")

    parser.add_argument(
        "--input_file", required=True,
        help="Full path to input file. Should be standard output file from "
        "'estimate substitution rates' tool in Galaxy using windows. Should "
        "have six columns: chrom, start, stop, n_sites, n_diffs, p_distance.")

    parser.add_argument(
        "--output_file", required=True,
        help="Full path to output file.")

    args = parser.parse_args()

    return args


def jc69(p_dist):
    print(p_dist)
    if p_dist == 0:
        return 0
    elif p_dist == 'NA':
        return 'NA'
    elif p_dist == 1:
        return 1
    else:
        term1 = 1.0 - ((4.0 / 3.0) * p_dist)
        jc_dist = -0.75 * math.log(term1)
        return jc_dist


def main():
    args = parse_args()

    with open(args.input_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        with open(args.output_file, "w") as o:
            writer = csv.writer(o, delimiter="\t")
            for i in reader:
                print(i)
                if i[0][0] == "#":
                    continue
                else:
                    if i[5] == 'NA':
                        p = i[5]
                    else:
                        p = float(i[5])
                    jc = jc69(p)
                    i.append(jc)
                    writer.writerow(i)


if __name__ == "__main__":
    main()
