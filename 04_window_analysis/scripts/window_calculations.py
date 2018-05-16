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
parser.add_argument("--windows", required=True,
                    help="Windowed Interval file defining how the analysis" +
                    " should be split up across the current chromosome.")
parser.add_argument("--callable", required=True,
                    help="BED file consisting of filtered callable sites" +
                    " distributed across the current chromosome.")
parser.add_argument("--diversity", required=True,
                    help="BED file containing all filtered SNPs and pi" +
                    " calculated at each site.")
parser.add_argument("--output", nargs="?", default=True,
                    help="Path to output file. If no output file is" +
                    " provided, default is standard out.")
parser.add_argument("--sliding", action='store_true', help="This indicates" +
                    " that the window file contains sliding windows, and" +
                    " slightly alters the algorithm that is used. This" +
                    " method is more inefficient than one involving" +
                    " nonoverlapping windows")

# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()
print(args)

# Script ----------------------------------------------------------------------
# Open Window Coordinates File and read to a list
with open(args.windows, 'rU') as f:
    window_coordinates = list(csv.reader(f, delimiter='\t'))

# Open diversity file
diversity = csv.reader(open(args.diversity, 'rU'), delimiter='\t')

# Open callable sites file
callable = csv.reader(open(args.callable, 'rU'), delimiter='\t')

# initialize an empty data array
data = []

# loop through the window coordinates file
if args.sliding is not True:

    # the first row of each file is initialized
    c = next(callable)
    d = next(diversity)

    # use a variable to save sites that stradle a window Boundary
    last_window_calls = 0

    for w in window_coordinates:
        # for each window keep track of the total called sites and SNP count
        sum_called = last_window_calls    # if there is overlap over window
        last_window_calls = 0
        count = 0
        pi_sum = 0

        # analyze the callable sites intervals
        while int(c[1]) < int(w[2]):
            if int(c[1]) >= int(w[1]) and int(c[2]) <= int(w[2]):
                # add the length of the callable site to the sum
                called = int(c[2]) - int(c[1])
                sum_called += called
                # set c to next line of callable and return to top of loop
                try:
                    c = next(callable)
                    continue
                except StopIteration:
                    break
            if int(c[1]) >= int(w[1]) and int(c[2]) > int(w[2]):
                called = int(w[2]) - int(c[1])
                sum_called += called
                # remember the length of the site that overlaps the next window
                last_window_calls = int(c[2]) - int(w[2])
                # set c to next line of callable and return to top of loop
                try:
                    c = next(callable)
                    continue
                except StopIteration:
                    break
            else:
                # in case the callable sites file stops before the window file
                # this block will keep advancing the script instead of catching
                # inside of an infinite while loop
                called = 0
                sum_called += 0
                break

        # analyze the diversity by site file
        while int(d[2]) > int(w[1]) and int(d[2]) <= int(w[2]):
            pi_sum += float(d[3])
            count += 1
            try:
                d = next(diversity)
            except StopIteration:
                break

        # use information to get windowed diversity
        if sum_called > 0:
            data.append([w[0], w[1], w[2], float(pi_sum / sum_called),
                        sum_called, count])
        else:
            data.append([w[0], w[1], w[2], pi_sum, "NA"])

# if --sliding is provided in command, use old algorithm
else:
    # loop through the window file
    for w in window_coordinates:
        sum_called = 0
        count = 0
        pi_sum = float(0)

        # loop through the callable sites file and find all in given window
        for c in callable:
            called = 0
            if int(c[1]) >= int(w[1]) and int(c[2]) <= int(w[2]):
                called = int(c[2]) - int(c[1])
                sum_called += called
            if int(c[1]) >= int(w[1]) and int(c[1]) < int(w[2]) and \
                    int(c[2]) > int(w[2]):
                called = int(w[2]) - int(c[1])
                sum_called += called
            if int(c[1]) < int(w[1]) and int(c[2]) <= int(w[2]) and \
                    int(c[2]) > int(w[1]):
                called = int(c[2]) - int(w[1])
                sum_called += called
            else:
                sum_called += 0

        # loop through diversity file and sum all sites in given window
        for d in diversity:
            if int(d[2]) > int(w[1]) and int(d[2]) <= int(w[2]):
                pi_sum += float(d[3])
                count += 1

        # print results to data array
        if sum_called > 0:
            data.append(['chrX', w[1], w[2], float(pi_sum / sum_called),
                         sum_called, count])
        else:
            data.append(['chrX', w[1], w[2], pi_sum, "NA"])


# write the results to output_file or standard out depending in args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in data:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in data:
            writer.writerow(row)
