"""
window_calculations.py
Daniel Cotter

use filtered diversity file, callable sites file, and window file to calculate
diversity in windows across a given chromosome
-------------------------------------------------------------------------------

"""

# import required modules -----------------------------------------------------
from sys import argv
import csv

# argument variables ----------------------------------------------------------
script, pop_diversity, called_sites, windows, output_file = argv

# Script ----------------------------------------------------------------------
# Open Window Coordinates File and read to a list
with open(windows, 'rU') as f:
    window_coordinates = list(csv.reader(f, delimiter='\t'))

# Open diversity file
diversity = csv.reader(open(pop_diversity, 'rU'), delimiter='\t')

# Open callable sites file
callable = csv.reader(open(called_sites, 'rU'), delimiter='\t')

# initialize an empty data array
data = []

# the first row of each file is initialized
c = next(callable)
d = next(diversity)

# use a variable to save sites that stradle a window Boundary
last_window_calls = 0

# loop through the window coordinates file
for w in window_coordinates:
    # for each window keep track of the total called sites and count of SNPs
    sum_called = last_window_calls    # if there is overlap from last window
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

    """
    When running through the Y chromosome the script gets caught in an infinite
    loop from what appears to be the for loop advancing to the next iteration
    before it is required. This causes the if statements in the callable block
    to not catch anything and the while loop runs indefinitely until ended. The
    coordinates that this starts at are as follows:

    w:  28800000 28900000   c:  28800008 28800044
    w:  28800000 28900000   c:  28801765 28801841
    w:  28800000 28900000   c:  28802050 28802057
    w:  28900000 29000000   c:  28802050 28802057

    """
    print(c[1] + " < " + w[2])

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

# write the results to output_file
with open(output_file, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    for row in data:
        writer.writerow(row)
