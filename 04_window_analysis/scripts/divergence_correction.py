"""
divergence_correction.py
Daniel Cotter

take a filtered substitution rate file and use it to correct the output of
the window_calculations.py script
-------------------------------------------------------------------------------

"""

# import required modules -----------------------------------------------------
import sys
import csv
import argparse

# parse arguments -------------------------------------------------------------
parser = argparse.ArgumentParser(description="correct for divergence")

# Parse the command line
parser.add_argument("--diversity", required=True,
                    help="Diversity BED file.")
parser.add_argument("--divergence", required=True,
                    help="Divergence BED file of same length.")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()

# Script ----------------------------------------------------------------------
# Open files and read to lists
with open(args.diversity, 'r') as f:
    diversity = list(csv.reader(f, delimiter='\t'))

with open(args.divergence, 'r') as f:
    divergence = list(csv.reader(f, delimiter='\t'))

merged = []
for i in range(len(diversity)):
    a = diversity[i][3]
    b = divergence[i][6]
    if a == 'NA' or b == 'NA':
        newVal = 'NA'
    elif float(a) == 0 or float(b) == 0:
        newVal = 0
    else:
        newVal = float(a) / float(b)
    merged.append([diversity[i][0], diversity[i][1], diversity[i][2],
                   newVal, diversity[i][4], diversity[i][5]])

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in merged:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in merged:
            writer.writerow(row)
