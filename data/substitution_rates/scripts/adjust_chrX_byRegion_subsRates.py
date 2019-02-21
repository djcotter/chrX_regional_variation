"""
adjust_chrX_byRegion_subsRates.py
Daniel Cotter

Takes a chrX byRegion substitution rates file from Galaxt and combines the
two separate nonPAR
-------------------------------------------------------------------------------
"""
# import required modules
import sys
import csv

# parse argvs
script, inputFile, outputFile, part1, part2 = sys.argv
part1 = int(part1)  # the index of the row with nonPAR part1
part2 = int(part2)  # the index of the row with nonPAR part2

# read the file into memory
with open(inputFile, 'r') as f:
    startFile = []
    for line in f:
        startFile.append(line.strip().split('\t'))

L = int(startFile[part1][3]) + int(startFile[part2][3])
N = float(startFile[part1][4]) + float(startFile[part2][4])
newLine = [startFile[0][0],
           startFile[part1][1],
           startFile[part2][2],
           L, N, N / L]

# create a new file with the changed line
newFile = [newLine]
for i, x in enumerate(startFile):
    if not (i == part1 or i == part2):
        newFile.append(x)

# sort the newFile on the chr start position
newFile = sorted(newFile, key=lambda l: int(l[1]))

# output file
with open(outputFile, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
    for line in newFile:
        writer.writerow(line)
