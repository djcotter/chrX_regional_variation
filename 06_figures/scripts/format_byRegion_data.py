"""
format_byRegion_data.py
Daniel Cotter

combine results of byRegion and wholeChr analyses to make a table of
populations and their associated diversity in each region
-------------------------------------------------------------------------------

"""

# import required modules -----------------------------------------------------
import sys
import csv
import argparse
from os.path import basename

# parse arguments -------------------------------------------------------------
parser = argparse.ArgumentParser(description="Combine all byRegion and "
                                             "wholeChr data.")

# Parse the command line
parser.add_argument("--chrX_byRegion", nargs='+',
                    help="byRegion data for chrX")
parser.add_argument("--chrY_wholeChr", nargs='+',
                    help="byRegion data for chrX")
parser.add_argument("--chr8_wholeChr", nargs='+',
                    help="byRegion data for chrX")
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
pop_keys = {'ACB': 'AFR',
            'YRI': 'AFR',
            'ASW': 'AFR',
            'ESN': 'AFR',
            'MSL': 'AFR',
            'GWD': 'AFR',
            'LWK': 'AFR',
            'MXL': 'AMR',
            'PUR': 'AMR',
            'CLM': 'AMR',
            'PEL': 'AMR',
            'CHB': 'EAS',
            'JPT': 'EAS',
            'CHS': 'EAS',
            'CDX': 'EAS',
            'KHV': 'EAS',
            'CEU': 'EUR',
            'TSI': 'EUR',
            'FIN': 'EUR',
            'GBR': 'EUR',
            'IBS': 'EUR',
            'GIH': 'SAS',
            'PJL': 'SAS',
            'BEB': 'SAS',
            'STU': 'SAS',
            'ITU': 'SAS'}
temp_data = {}

# append all temp data info to list
for item in args.chrX_byRegion:
    with open(item, 'rU') as f:
        temp = list(csv.reader(f, delimiter='\t'))
    temp_data[basename(item)[0:3]] = []
    for line in temp:
        temp_data[basename(item)[0:3]].append(line[3])
        temp_data[basename(item)[0:3]].append(line[6])
        temp_data[basename(item)[0:3]].append(line[7])

for item in args.chrY_wholeChr:
    with open(item, 'rU') as f:
        temp = list(csv.reader(f, delimiter='\t'))
    temp_data[basename(item)[0:3]].append(float(temp[0][3]))

for item in args.chr8_wholeChr:
    with open(item, 'rU') as f:
        temp = list(csv.reader(f, delimiter='\t'))
    temp_data[basename(item)[0:3]].append(float(temp[0][3]))

data = [['SUPERPOP', 'POP', 'PAR1', 'nonPAR', 'XTR', 'PAR2', 'chrY', 'chr8']]
for key in temp_data:
    newline = [pop_keys[key], key]
    for num in temp_data[key]:
        newline.append(num)
    data.append(newline)

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
