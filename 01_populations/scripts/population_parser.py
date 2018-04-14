from sys import argv
import csv

# usage is population_parser.py input_file out_directory
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

file1 = output_directory + 'pop_table.txt'
with open(file1, 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['POP', 'Females', 'Males'])
    for pop in pop_dict_females:
        writer.writerow([pop, len(pop_dict_females[pop]),
                         len(pop_dict_males[pop])])

file2 = output_directory + 'subpop_table.txt'
with open(file2, 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['POP', 'Females', 'Males'])
    for pop in subpop_dict_females:
        writer.writerow([pop, len(subpop_dict_females[pop]),
                         len(subpop_dict_males[pop])])

outpath1 = output_directory + 'pops/'
for pop in pop_dict_females:
    temp_path = outpath1 + pop + '_females'
    with open(temp_path, 'wb') as file:
        for samp in pop_dict_females[pop]:
            file.write(samp + '\n')
    temp_path2 = outpath1 + pop + '_males'
    with open(temp_path2, 'wb') as file:
        for samp in pop_dict_males[pop]:
            file.write(samp + '\n')
    temp_path3 = outpath1 + pop + '_individuals'
    with open(temp_path3, 'wb') as file:
        for samp in pop_dict_females[pop]:
            file.write(samp + '\n')
        for samp in pop_dict_males[pop]:
            file.write(samp + '\n')

outpath2 = output_directory + 'subpops/'
for pop in subpop_dict_females:
    temp_path = outpath2 + pop + '_females'
    with open(temp_path, 'wb') as file:
        for samp in subpop_dict_females[pop]:
            file.write(samp + '\n')
    temp_path2 = outpath2 + pop + '_males'
    with open(temp_path2, 'wb') as file:
        for samp in subpop_dict_males[pop]:
            file.write(samp + '\n')
    temp_path3 = outpath2 + pop + '_individuals'
    with open(temp_path3, 'wb') as file:
        for samp in subpop_dict_females[pop]:
            file.write(samp + '\n')
        for samp in subpop_dict_males[pop]:
            file.write(samp + '\n')
