"""
window_calculations.py
Daniel Cotter

filter the window analysis output by removing windows where the number of
callable sites is less than a given percentage of the window size
-------------------------------------------------------------------------------

"""

# import required modules -----------------------------------------------------
import sys
import csv
import argparse

# parse arguments -------------------------------------------------------------
parser = argparse.ArgumentParser(description="Filter window anlysis output.")

parser.add_argument("--input", required=True,
                    help="Path to the input BED file")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")
parser.add_argument("--windowSize", nargs='?', type=int, required=True,
                    help="Size of the windows in the file")
parser.add_argument("--filter", nargs='?', type=float, default=0.1,
                    help="The size of the filter. Default is 0.1 (10%)")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()

# script ----------------------------------------------------------------------
# set the filter to use
filter = args.windowSize * args.filter

# read the input file
with open(args.input, 'r') as f:
    window_file = list(csv.reader(f, delimiter='\t'))

window_file[14][4] == 'NA'
# loop through input and mask windows that don't pass filter
for line in window_file:
    if line[4] == 'NA':
        line[3] = 'NA'
        line[4] = 'NA'
        line[5] = 'NA'
    elif int(line[4]) < filter:
        line[3] = 'NA'
        line[4] = 'NA'
        line[5] = 'NA'

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in window_file:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in window_file:
            writer.writerow(row)
