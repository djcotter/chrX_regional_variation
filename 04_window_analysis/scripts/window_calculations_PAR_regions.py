from sys import argv
import csv

script, pop_diversity, called_sites, output_file = argv

# Open the window file as bed coordinates so that every window_coordinates[i][1] corresponds to the start position and [2] is the end position.
# with open(windows, 'rU') as f:
# window_coordinates = list(csv.reader(f, delimiter = '\t'))

# The following line establishes window coordinates that correspond to the various regions across the X chromsome
# These coordinates
window_coordinates = [["PAR1", 60001, 2699520], ["nonPAR1", 2699520, 88193855], [
    "XTR", 88193855, 93193855], ["nonPAR2", 93193855, 154931044], ["PAR2", 154931044, 155260560]]

with open(pop_diversity, 'rU') as f:
    diversity = list(csv.reader(f, delimiter='\t'))

with open(called_sites, 'rU') as f:
    callable = list(csv.reader(f, delimiter='\t'))

data = []

for c in window_coordinates:
    sum_called = 0
    for r in callable:
        called = 0
        if int(r[1]) >= int(c[1]) and int(r[2]) <= int(c[2]):
            called = int(r[2]) - int(r[1])
            sum_called += called
        if int(r[1]) >= int(c[1]) and int(r[1]) < int(c[2]) and int(r[2]) > int(c[2]):
            called = int(c[2]) - int(r[1])
            sum_called += called
        if int(r[1]) < int(c[1]) and int(r[2]) <= int(c[2]) and int(r[2]) > int(c[1]):
            called = int(r[2]) - int(c[1])
            sum_called += called
        else:
            sum_called += 0
    pi_sum = float(0)
    for d in diversity:
        if int(d[2]) > int(c[1]) and int(d[2]) <= int(c[2]):
            pi_sum += float(d[3])
    if sum_called > 0:
        data.append([c[0], c[1], c[2], float(pi_sum), sum_called])
    else:
        data.append(['chrX', c[1], c[2], pi_sum, "NA"])

nonPAR_diversity = (float(data[1][3] + data[3][3])/(data[1][4] + data[3][4]))
nonPAR_called = float(data[1][4] + data[3][4])

formatted_data = [['PAR1', data[0][1], data[0][2], float(data[0][3]/data[0][4]), data[0][4]], ['nonPAR', data[1][1], data[3][2], nonPAR_diversity, nonPAR_called], [
    'XTR', data[2][1], data[2][2], float(data[2][3]/data[2][4]), data[2][4]], ['PAR2', data[4][1], data[4][2], float(data[4][3]/data[4][4]), data[4][4]]]

with open(output_file, 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    for info in formatted_data:
        writer.writerow(info)
    writer.writerow("")
    writer.writerow(
        ["*The nonPAR region spans the XTR but it excludes all data points in the XTR region"])
