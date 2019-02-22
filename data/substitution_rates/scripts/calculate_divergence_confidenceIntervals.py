"""
calculate_divergence_confidenceIntervals.py
Daniel Cotter

given a wholeChr or byRegion divergence file and a 100kb divergence file
calculate bootstrapped confidence intervals for each region/wholeChr
-------------------------------------------------------------------------------
"""
# import required modules -----------------------------------------------------
import argparse
import sys
import csv
import numpy as np
import re

# maybe go back and use actual mean ???

# Bootstrap functions ---------------------------------------------------------
def rand_samples(data_array):
    """
    Create a randomly dsitributed set of data based on a given array and
    the length of the given data array
    """
    array_len = len(data_array)
    if array_len > 0:
        indices = np.random.randint(array_len - 1, size=array_len)
        return [data_array[i] for i in indices]
    else:
        return []


def bootstrap_CI_mean(data_array, replicates):
    """
    Bootstraps the data to get a 95% confidence interval for the region/chr
    """
    resamples = []
    for i in range(replicates):
        samples = rand_samples(data_array)
        if samples:
            sums = np.sum(samples, axis=0)
            resamples.append(sums[1] / sums[0])

    return [np.mean(resamples),
            np.nanpercentile(resamples, 2.5),
            np.nanpercentile(resamples, 97.5)]


# parse the command line ------------------------------------------------------
parser = argparse.ArgumentParser(description="Calculate windowed diversity" +
                                 " across a chromosome.")

# Parse the command line
parser.add_argument("--regionFile", required=True,
                    help="byRegion or wholeChr file")
parser.add_argument("--windowFile", required=True,
                    help="100kb window file")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")
parser.add_argument("--replicates", nargs='?', type=int, default=100,
                    help="number of times the permutation should be run. " +
                    "Default is 100.")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()

# script ----------------------------------------------------------------------

# determine if byRegion or wholeChr file
pattern = re.compile(r'byRegion')
if pattern.search(args.regionFile):
    chrX_file = True
    # define the adjusted coordinates for the chrX regions
    # (add telomeres in for the sake of the following window comparisons)
    wc = {"PAR1": [0, 2699520],
          "nonPAR1": [2699520, 88193855],
          "XTR": [88193855, 93193855],
          "nonPAR2": [93193855, 154931044],
          "PAR2": [154931044, 155270560]}
else:
    chrX_file = False

# read the byRegion or wholeChr file into memory
with open(args.regionFile, 'r') as f:
    byRegion = list(csv.reader(f, delimiter='\t'))

if chrX_file:
    # for line in region_file grab the mean for each region
    # this file should be in the order PAR1, nonPAR, XTR, PAR2
    mean_dict = {}
    for line, tag in zip(byRegion, ['PAR1', 'nonPAR', 'XTR', 'PAR2']):
        mean_dict[tag] = float(line[5])

# read the window file into memory
with open(args.windowFile, 'r') as f:
    byWindow = list(csv.reader(f, delimiter='\t'))

# if it is a chrX file
if chrX_file:
    new_wc = []
    nonParVals = []
    for key in wc.keys():
        region = wc[key]
        myVals = []
        for x in byWindow:
            if (int(x[1]) + int(x[2])) / 2 > region[0] and \
                    (int(x[1]) + int(x[2])) / 2 <= region[1]:
                if not x[5] == 'NA':
                    myVals.append([int(x[3]), float(x[4])])
        if key == 'nonPAR1' or key == 'nonPAR2':
            nonParVals = nonParVals + myVals
        else:
            new_wc.append([key] + region +
                          bootstrap_CI_mean(myVals, args.replicates))
    new_wc.append(['nonPAR', wc['nonPAR1'][0], wc['nonPAR2'][1]] +
                  bootstrap_CI_mean(nonParVals, args.replicates))
    results = new_wc
# if not a chrX file
else:
    myVals = []
    for x in byWindow:
        if not x[5] == 'NA':
            myVals.append([int(x[3]), float(x[4])])
    results = byRegion[0][0:3] + bootstrap_CI_mean(myVals, args.replicates)

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in results:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in results:
            writer.writerow(row)
