"""
window_calculations.py
Daniel Cotter

use filtered diversity file, callable sites file, and window file to calculate
diversity in windows across a given chromosome
-------------------------------------------------------------------------------

"""

# import required modules -----------------------------------------------------
import sys
import csv
import argparse

# parse arguments -------------------------------------------------------------
parser = argparse.ArgumentParser(description="Calculate windowed diversity" +
                                 " across a chromosome.")

# Parse the command line
parser.add_argument("--chrX", required=True,
                    help="chrX byRegion values")
parser.add_argument("--chrY", required=True,
                    help="chrY wholeChr file")
parser.add_argument("--chr8", required=True,
                    help="chr8 wholeChr file")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()

# script ----------------------------------------------------------------------
# read input files into script
with open(args.chrX, 'r') as f:
    chrX = list(csv.reader(f, delimiter='\t'))

with open(args.chrY, 'r') as f:
    chrY = list(csv.reader(f, delimiter='\t'))

with open(args.chr8, 'r') as f:
    chr8 = list(csv.reader(f, delimiter='\t'))

# calculate ratios
data = []
for line in chrX:
    data.append(["{}/{}".format(line[0], 'A'),
                 float(line[3]) / float(chr8[0][3])])
data.append(["{}/{}".format(chrY[0][0], 'A'),
             float(chrY[0][3]) / float(chr8[0][3])])

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in data:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in data:
            writer.writerow(row)
