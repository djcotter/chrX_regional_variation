"""
permute_chrX_windows.py
Daniel Cotter

use the 100kb windows and the byRegion file to test for significance of the
difference between PAR1, PAR2, XTR, and nonPAR.
-------------------------------------------------------------------------------
"""

# import required modules -----------------------------------------------------
import argparse
import sys
import csv

# parse the command line ------------------------------------------------------
parser = argparse.ArgumentParser(description="Calculate windowed diversity" +
                                 " across a chromosome.")

# Parse the command line
parser.add_argument("--byWindow", required=True,
                    help="Windowed diversity calculated across chrX")
parser.add_argument("--byRegion", required=True,
                    help="Region diversity calculated across chrX")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")
parser.add_argument("--repetitions", nargs='?', type=int, default=100,
                    help="number of times the permutation should be run. " +
                    "Default is 100.")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()

# script ----------------------------------------------------------------------

# define the coordinates for the chrX regions
wc = [["PAR1", 60001, 2699520],
      ["nonPAR1", 2699520, 88193855],
      ["XTR", 88193855, 93193855],
      ["nonPAR2", 93193855, 154931044],
      ["PAR2", 154931044, 155260560]]

# read the windowed_diversity file into memory
with open(args.byWindow, 'r') as f:
    windowed_diversity = list(csv.reader(f, delimiter='\t'))

# Check how many windows are in each region and add all diversity to list
diversity = []
PAR1_windows = 0
PAR2_windows = 0
XTR_windows = 0

for window in windowed_diversity:
    # track PAR1 window
    if int(window[1]) < wc[0][2]:
        label = 'PAR1'
    # track PAR2 window
    elif int(window[1]) < wc[4][2] and int(window[1]) > wc[4][1] - 1e5:
        label = 'PAR2'
    # track XTR window
    elif int(window[1]) < wc[2][2] and int(window[1]) > wc[2][1] - 1e5:
        label = 'XTR'
    # track nonPAR
    else:
        label = 'nonPAR'

    # remove windows with no called site information
    if window[4] is not 'NA':
        diversity.append([label, window[3]])

if args.output is not False:
    for line in diversity:
        print(line)
