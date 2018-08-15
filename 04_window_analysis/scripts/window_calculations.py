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
import numpy as np

# parse arguments -------------------------------------------------------------
parser = argparse.ArgumentParser(description="Calculate windowed diversity" +
                                 " across a chromosome.")

# Parse the command line
parser.add_argument("--windows", required=False,
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
parser.add_argument("--chrX_windows", action='store_true', help="This flag" +
                    " indicates that the analysis should be performed on the" +
                    " X chromosome using windows that correspond to PAR1," +
                    " nonPAR, XTR, and PAR2. No windows need to be provided" +
                    " for this analysis.")
parser.add_argument("--replicates", nargs='?', type=int, default=100,
                    help="number of times the bootstrap " +
                    "should be run. Default is 100.")
parser.add_argument("--test", nargs="?", default=False, help="If specified," +
                    " the script will output results of a filter test to" +
                    " the provided path. Test1 calculates pi for the region" +
                    " defined as nonPAR + PARs + XTR. Test2 calculates" +
                    " pi for the region defined as nonPAR + XTR.")


# check that commands are there
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit()

# add arguments to args variable
args = parser.parse_args()


# Bootstrap functions ---------------------------------------------------------
def rand_samples(data_array):
    """
    Create a randomly dsitributed set of data based on a given array and
    the length of the given data array
    """
    array_len = len(data_array)
    if array_len > 0:
        indices = np.random.random_integers(0, array_len - 1, array_len)
        return np.take(data_array, indices)
    else:
        return np.asarray([0])


def bootstrap_CI_mean(data_array, replicates, called_sites):
    """
    Bootstraps the data to get a 95% confidence interval for the
    """
    resamples = []
    for i in range(replicates):
        resamples.append(np.sum(rand_samples(data_array)) / called_sites)

    return [np.nanpercentile(resamples, 2.5),
            np.nanpercentile(resamples, 97.5)]


# Script ----------------------------------------------------------------------
# Open Window Coordinates File and read to a list
if args.chrX_windows is not True:
    with open(args.windows, 'rU') as f:
        window_coordinates = list(csv.reader(f, delimiter='\t'))
else:
    # these windows correspond to GRCh37 (hg19)
    wc = [["PAR1", 60001, 2699520],
          ["nonPAR1", 2699520, 88193855],
          ["XTR", 88193855, 93193855],
          ["nonPAR2", 93193855, 154931044],
          ["PAR2", 154931044, 155260560]]

# Open diversity file
diversity = csv.reader(open(args.diversity, 'rU'), delimiter='\t')

# Open callable sites file
callable = csv.reader(open(args.callable, 'rU'), delimiter='\t')

# initialize an empty data array
data = []
bootstraps = []

# -----------------------------------------------------------------------------
# loop through the predefined regions for the X chromosome
if args.chrX_windows is True:

    # the first row of each file is initialized
    c = next(callable)
    d = next(diversity)
    # test 1 is nonPAR, PAR1, PAR2, XTR
    # test 2 is nonPAR, XTR
    # nonPAR from main results is test3 (with both filtered)
    test1 = {'pi_vals': [], 'called': 0, 'count': 0}
    test2 = {'pi_vals': [], 'called': 0, 'count': 0}

    # use a variable to save sites that straddle a window Boundary
    last_window_calls = 0

    for w in wc:
        # for each window keep track of the total called sites and SNP count
        sum_called = last_window_calls    # if there is overlap over window
        last_window_calls = 0
        count = 0
        pi_vals = []

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
            pi_vals.append(float(d[3]))
            count += 1
            try:
                d = next(diversity)
            except StopIteration:
                break

        # adds the information directly to the wc list
        if w[0] == "nonPAR1":
            nonPAR_vals = pi_vals
            nonPAR_called = sum_called
            nonPAR_count = count

            test1['pi_vals'] += pi_vals
            test1['called'] += sum_called
            test1['count'] += count

            test2['pi_vals'] += pi_vals
            test2['called'] += sum_called
            test2['count'] += count

        elif w[0] == "nonPAR2":
            nonPAR_vals = nonPAR_vals + pi_vals
            nonPAR_called += sum_called
            nonPAR_count += count

            test1['pi_vals'] += pi_vals
            test1['called'] += sum_called
            test1['count'] += count

            test2['pi_vals'] += pi_vals
            test2['called'] += sum_called
            test2['count'] += count

            bootstraps = bootstrap_CI_mean(np.asarray(nonPAR_vals),
                                           args.replicates, nonPAR_called)
            data.append(['nonPAR', wc[1][1], wc[3][2],
                         float(sum(nonPAR_vals) / nonPAR_called),
                         nonPAR_called, nonPAR_count, bootstraps[0],
                         bootstraps[1]])

        else:
            w.append(sum(pi_vals))  # index 3
            w.append(sum_called)  # index 4
            w.append(count)  # index 5
            if sum_called > 0 and len(pi_vals) > 0:
                bootstraps = bootstrap_CI_mean(np.asarray(pi_vals),
                                               args.replicates, sum_called)
                data.append([w[0], w[1], w[2],
                            float(sum(pi_vals) / sum_called),
                            sum_called, count, bootstraps[0],
                            bootstraps[1]])
            else:
                data.append([w[0], w[1], w[2], "NA",
                             w[4], w[5], "NA", "NA"])
            # keep track of tests
            if w[0] == 'PAR1' or w[0] == 'PAR2':
                test1['pi_vals'] += pi_vals
                test1['called'] += sum_called
                test1['count'] += count
            elif w[0] == 'XTR':
                test1['pi_vals'] += pi_vals
                test1['called'] += sum_called
                test1['count'] += count
                test2['pi_vals'] += pi_vals
                test2['called'] += sum_called
                test2['count'] += count

    # format of data is region, start, end, pi, called_sites, num variants
    # then the bounds of the bootstrapped confidence interval
    # the following will reorder the data to PAR nonPAR XTR PAR2
    data_order = [0, 2, 1, 3]
    data = [data[i] for i in data_order]


# -----------------------------------------------------------------------------
# the standard window analysis for windows provided in a non-sliding format
elif args.sliding is not True:

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
        pi_vals = []

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
            pi_vals.append(float(d[3]))
            count += 1
            try:
                d = next(diversity)
            except StopIteration:
                break

        # use information to get windowed diversity
        if sum_called > 0 and len(pi_vals) > 0:
            bootstraps = bootstrap_CI_mean(np.asarray(pi_vals),
                                           args.replicates, sum_called)
            data.append([w[0], w[1], w[2], float(sum(pi_vals) / sum_called),
                        sum_called, count, bootstraps[0],
                        bootstraps[1]])
        else:
            data.append([w[0], w[1], w[2], sum(pi_vals),
                        sum_called, count, "NA", "NA"])

# -----------------------------------------------------------------------------
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
            data.append(['chrX', w[1], w[2], pi_sum, "NA", "NA"])


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

if args.test is not False:
    bootstrap_test1 = bootstrap_CI_mean(np.asarray(test1['pi_vals']),
                                        args.replicates, test1['called'])
    bootstrap_test2 = bootstrap_CI_mean(np.asarray(test2['pi_vals']),
                                        args.replicates, test2['called'])
    data1 = []
    data1.append(['test1', wc[0][1], wc[4][2],
                  float(sum(test1['pi_vals']) / test1['called']),
                  test1['called'], test1['count'],
                  bootstrap_test1[0], bootstrap_test1[1]])
    data1.append(['test2', wc[1][1], wc[3][2],
                  float(sum(test2['pi_vals']) / test2['called']),
                  test2['called'], test2['count'],
                  bootstrap_test2[0], bootstrap_test2[1]])
    with open(args.test, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in data1:
            writer.writerow(row)
