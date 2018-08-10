"""
average_ld_by_window.py
Daniel Cotter
Updated: 3/30/17

The purpose of this script is to take pairwise LD output from plink
and calculate average LD in set windows
-------------------------------------------------------------------------------
"""

# import required packages ----------------------------------------------------
import sys
import csv
import argparse
import ld_analysis

# Set argument variables ------------------------------------------------------
""" input file: pairwise LD file generated in space-delimited table format
    by plink using the command 'with-freqs'
                (essential so the columns are correctly indexed)
    windows: bed file with desired windows for the analysis
    ld_bin_size: length in either direction from the focal position in which
                 we should consider LD (in kilobases)
    output_file: path to an output file
"""
# Parse arguments
parser = argparse.ArgumentParser(description="Filter window anlysis output.")

parser.add_argument("--plink_ld", required=True,
                    help="Path to the input plink LD file. This file should" +
                    "be generated using ")
parser.add_argument("--windows", required=False, default=False,
                    help="Path to the input window file.")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")
parser.add_argument("--binSize", nargs='?', type=int, required=True,
                    help="Size of the bin that LD is calculated in. This " +
                    "extends from the current site being analyzed in each " +
                    "direction. This is in kilobases.")
parser.add_argument("--byRegion", action='store_true', help="The script " +
                    "will analyze the established chrX regions.")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()

# Script ----------------------------------------------------------------------

# Read in windows file and set starting conditions
# Open the window file as bed coordinates so that every
# window_coordinates[i][1] corresponds to the start position
# and [2] to the end position.
if args.windows is not False:
    with open(args.windows, 'rU') as f:
        window_file = list(csv.reader(f, delimiter='\t'))
        for window in window_file:
            window[1] = float(window[1])
            window[2] = float(window[2])

# Open the LD file and perform the analysis line by line
if args.byRegion is True:
    ld_calculations = ld_analysis.LD_loop_byRegion(
        args.plink_ld, args.binSize)
else:
    ld_calculations = ld_analysis.LD_loop(
        args.plink_ld, window_file, args.binSize)

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in ld_calculations:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in ld_calculations:
            writer.writerow(row)
