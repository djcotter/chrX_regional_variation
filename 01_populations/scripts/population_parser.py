"""
population_parser.py
Daniel Cotter

Take a reference panel from 1000 genomes and parse by pops and gender
---------------------------------------------------------------------

usage is population_parser.py input_file output_directory
"""

# import statements

from sys import argv
import csv
import pathlib

# -----------------------------------------------------------------------------
# declare argument variables
script, input_file, output_directory = argv

with open(input_file, 'rU') as f:
    samples = list(csv.reader(f, delimiter='\t'))

# pop_dict = {}
# for row in samples:
    # if row[2] == 'super_pop':
    # next
    # elif row[2] in pop_dict:
    # pop_dict[row[2]].append(row)
    # else:
    # pop_dict[row[2]] = [row]

males = []
females = []
for samp in samples:
    if samp[3] == 'male':
        males.append(samp[0])
    elif samp[3] == 'female':
        females.append(samp[0])

# create dictionaries with seperate lists of individuals by pop
pop_dict_males = {}
pop_dict_females = {}
for samp in samples:
    if samp[3] == 'male':
        if samp[2] in pop_dict_males:
            pop_dict_males[samp[2]].append(samp[0])
        else:
            pop_dict_males[samp[2]] = [samp[0]]
    elif samp[3] == 'female':
        if samp[2] in pop_dict_females:
            pop_dict_females[samp[2]].append(samp[0])
        else:
            pop_dict_females[samp[2]] = [samp[0]]

# create dictionaries with seperate lists of individuals by subpop
subpop_dict_males = {}
subpop_dict_females = {}
for samp in samples:
    if samp[3] == 'male':
        if samp[1] in subpop_dict_males:
            subpop_dict_males[samp[1]].append(samp[0])
        else:
            subpop_dict_males[samp[1]] = [samp[0]]
    elif samp[3] == 'female':
        if samp[1] in subpop_dict_females:
            subpop_dict_females[samp[1]].append(samp[0])
        else:
            subpop_dict_females[samp[1]] = [samp[0]]

# print a table with the number of males and females by super_pop
file1 = output_directory + 'pop_table.txt'
with open(file1, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['POP', 'Females', 'Males'])
    for pop in pop_dict_females:
        writer.writerow([pop, len(pop_dict_females[pop]),
                         len(pop_dict_males[pop])])

# print a table with the number of males and females by super_pop
file2 = output_directory + 'subpop_table.txt'
with open(file2, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['POP', 'Females', 'Males'])
    for pop in subpop_dict_females:
        writer.writerow([pop, len(subpop_dict_females[pop]),
                         len(subpop_dict_males[pop])])

# set the output subdirectory and create it if it doesn't exist
outpath1 = output_directory + 'populations/'
pathlib.Path(outpath1).mkdir(parents=True, exist_ok=True)

# loop through files and print parsed lists of individuals
for pop in pop_dict_females:
    temp_path = outpath1 + pop + '_females.txt'
    with open(temp_path, 'w') as file:
        for samp in pop_dict_females[pop]:
            file.write(samp + '\n')
    temp_path2 = outpath1 + pop + '_males.txt'
    with open(temp_path2, 'w') as file:
        for samp in pop_dict_males[pop]:
            file.write(samp + '\n')
    temp_path3 = outpath1 + pop + '_individuals.txt'
    with open(temp_path3, 'w') as file:
        for samp in pop_dict_females[pop]:
            file.write(samp + '\n')
        for samp in pop_dict_males[pop]:
            file.write(samp + '\n')

# set the output subdirectory and create it if it doesn't exist
outpath2 = output_directory + 'subpopulations/'
pathlib.Path(outpath2).mkdir(parents=True, exist_ok=True)

# loop through the files and print parsed lists of individuals
for pop in subpop_dict_females:
    temp_path = outpath2 + pop + '_females.txt'
    with open(temp_path, 'w') as file:
        for samp in subpop_dict_females[pop]:
            file.write(samp + '\n')
    temp_path2 = outpath2 + pop + '_males.txt'
    with open(temp_path2, 'w') as file:
        for samp in subpop_dict_males[pop]:
            file.write(samp + '\n')
    temp_path3 = outpath2 + pop + '_individuals.txt'
    with open(temp_path3, 'w') as file:
        for samp in subpop_dict_females[pop]:
            file.write(samp + '\n')
        for samp in subpop_dict_males[pop]:
            file.write(samp + '\n')

# write a file with a list of all individuals parsed only by genders
temp_path = outpath1 + 'ALL' + '_females.txt'
with open(temp_path, 'w') as file:
    for samp in females:
        file.write(samp + '\n')

temp_path = outpath1 + 'ALL' + '_males.txt'
with open(temp_path, 'w') as file:
    for samp in males:
        file.write(samp + '\n')

temp_path = outpath1 + 'ALL' + '_individuals.txt'
with open(temp_path, 'w') as file:
    for samp in females:
        file.write(samp + '\n')
    for samp in males:
        file.write(samp + '\n')
