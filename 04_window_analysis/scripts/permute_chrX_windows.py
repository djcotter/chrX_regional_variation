"""
permute_chrX_windows.py
Daniel Cotter

use the 100kb windows and the byRegion file to test for significance of the
difference between PAR1, PAR2, XTR, and nonPAR.
creates a copy of the diversity by region file and adds a col with P values
-------------------------------------------------------------------------------
"""

# import required modules -----------------------------------------------------
import argparse
import sys
import csv
import numpy as np

# parse the command line ------------------------------------------------------
parser = argparse.ArgumentParser(description="Calculate windowed diversity"
                                 " across a chromosome.")

# Parse the command line
parser.add_argument("--byWindow", required=True,
                    help="Windowed diversity calculated across chrX")
parser.add_argument("--byRegion", required=True,
                    help="Region diversity calculated across chrX")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is"
                    " provided, default is standard out.")
parser.add_argument("--replicates", nargs='?', type=int, default=100,
                    help="number of times the permutation should be run. "
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

# read the byRegion_diversity file into memory
with open(args.byRegion, 'r') as f:
    byRegion_diversity = list(csv.reader(f, delimiter='\t'))

# Check how many windows are in each region and add all diversity to list
diversity = []
PAR1_windows = 0
PAR1_mean = 0
PAR2_windows = 0
PAR2_mean = 0
XTR_windows = 0
XTR_mean = 0
nonPAR_windows = 0
nonPAR_mean = 0

for window in windowed_diversity:
    # skip lines that do not have info or were filtered out
    if window[4] == 'NA' or window[3] == 'NA':
        continue

    # track PAR1 window
    if int(window[1]) < wc[0][2]:
        label = 'PAR1'
        PAR1_windows += 1
        PAR1_mean += float(window[3])
    # track PAR2 window
    elif int(window[1]) < wc[4][2] and int(window[1]) > wc[4][1] - 1e5:
        label = 'PAR2'
        PAR2_windows += 1
        PAR2_mean += float(window[3])
    # track XTR window
    elif int(window[1]) < wc[2][2] and int(window[1]) > wc[2][1] - 1e5:
        label = 'XTR'
        XTR_windows += 1
        XTR_mean += float(window[3])
    # track nonPAR
    else:
        label = 'nonPAR'
        nonPAR_windows += 1
        nonPAR_mean += float(window[3])

    # append diversity information to new array
    diversity.append([label, float(window[3])])

# get the initial values for PAR1, nonPAR, XTR, and PAR2
PAR1_mean = PAR1_mean / PAR1_windows
PAR2_mean = PAR2_mean / PAR2_windows
XTR_mean = XTR_mean / XTR_windows
nonPAR_mean = nonPAR_mean / nonPAR_windows

# add all diversity values to a numpy array
div_vals = np.array([x[1] for x in diversity])

# establish cutoffs for vals in diversity array
PAR1_end = PAR1_windows
nonPAR_end = PAR1_end + nonPAR_windows
XTR_end = nonPAR_end + XTR_windows
PAR2_end = XTR_end + PAR2_windows

# base counts for the permutation test
PAR1_count = 0
PAR2_count = 0
XTR_count = 0

# permute the array and test the difference between the region means
for i in range(args.replicates):
    PAR1_test = np.average(
        np.random.permutation(div_vals)[0:PAR1_end])
    nonPAR_test = np.average(
        np.random.permutation(div_vals)[PAR1_end:nonPAR_end])
    XTR_test = np.average(
        np.random.permutation(div_vals)[nonPAR_end:XTR_end])
    PAR2_test = np.average(
        np.random.permutation(div_vals)[XTR_end:PAR2_end])

    if (PAR1_test - nonPAR_test) >= (PAR1_mean - nonPAR_mean):
        PAR1_count += 1
    if (PAR2_test - nonPAR_test) >= (PAR2_mean - nonPAR_mean):
        PAR2_count += 1
    if (XTR_test - nonPAR_test) >= (XTR_mean - nonPAR_mean):
        XTR_count += 1

# calculate P values
PAR1_p = float(PAR1_count) / args.replicates
if PAR1_p == 0:
    PAR1_p = "< {}".format(1 / args.replicates)

PAR2_p = float(PAR2_count) / args.replicates
if PAR2_p == 0:
    PAR2_p = "< {}".format(1 / args.replicates)

XTR_p = float(XTR_count) / args.replicates
if XTR_p == 0:
    XTR_p = "< {}".format(1 / args.replicates)

# add P values to byRegion file
byRegion_diversity[0].append(PAR1_p)  # PAR1
byRegion_diversity[1].append('------')  # nonPAR
byRegion_diversity[2].append(XTR_p)  # XTR
byRegion_diversity[3].append(PAR2_p)  # PAR2


# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in byRegion_diversity:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in byRegion_diversity:
            writer.writerow(row)
