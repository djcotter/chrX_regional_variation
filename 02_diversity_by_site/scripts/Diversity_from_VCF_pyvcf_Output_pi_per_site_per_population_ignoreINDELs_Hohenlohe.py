# This script takes a list of populations and outputs a file for each
# population listing pi per site
# Written by Tim Webster
# Last updated on April 19, 2016 by Tim Webster
# Edited by Daniel Cotter 4/9/18
# update to Python 3.5 compatibility - 4/12/18


#############################
# Import required modules   #
#############################
"""
- argparse, collections, csv, and sys are part of the standard Python library
- biopython, numpy, sympy, and vcf are easily installed via conda and pip
"""

import argparse
import collections
import csv
import sys
# from Bio import SeqIO
import numpy as np
import sympy
import vcf
from os import path

########################################
# Import files from the command line   #
########################################
parser = argparse.ArgumentParser(description="UPDATE THIS")

# Print help/usage if no arguments are supplied
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)


# Parse the command line
parser.add_argument("--vcf", nargs='?', default=sys.stdin,
                    help="Input VCF file. Can be gzipped. Reads from stdin" +
                    " by default")
parser.add_argument("--population_lists", nargs="*", default=None,
                    help="Default is None. If None, script will treat all" +
                    " individuals as coming from the same population." +
                    " Else, will read files and treat all individuals in" +
                    " each file as coming from the same population.")
parser.add_argument("--chrom_inc", default=None,
                    help="Default is None. Provide the name of the scaffold" +
                    " or chromosome to be analyzed.  If none provided," +
                    " will output all sites in VCF")
parser.add_argument("--out_directory", nargs='?', default="",
                    help="Default is current directory. By specifying output" +
                    " directory, all files will be generated in user-" +
                    "provided location.")
parser.add_argument("--haploid", action='store_true', help="This flag" +
                    " indicates that the chromosome has haploid sites and" +
                    " grabs one allele if a call is haploid (e.g. returns" +
                    " \"A/T\"). Only use this flag if the chromosome being" +
                    " processed is all haploid or partly (e.g. nonPAR of X).")
args = parser.parse_args()
print("Parsing arguments from command line...")
print("")

###########################################################################

# Update this to follow the format of the population list reader below
# and grab names from the command line ###########
print("Reading input files...")
print("")

# Open vcf file and initiate the CyVCF parser
vcf_reader = vcf.Reader(filename=args.vcf)

# Process input population lists
# Includes a number of cleaning and filtering steps to clean up accidental
# whitespace, and read from both vertical and horizonal lists (horizontal have
# to be tab separated) If no input lists are provided, it treats all samples
# in the VCF as coming from the same population

if args.population_lists is not None:
    populations = []
    no_populations = False
    for i in args.population_lists:
        with open(i, "r") as f:
            temp_pop = [item for sublist in list(csv.reader(
                f, delimiter="\t")) for item in sublist]
            temp_pop = [x.strip() for x in temp_pop]
            while "" in temp_pop:
                temp_pop.remove("")
        populations.append([temp_pop, i, []])
else:
    populations = [[vcf_reader.samples, args.vcf, []]]


########################################################################
# Functions to calculate diversity and conduct bootstrap resampling  #
########################################################################

def pi_overall(tot_diff, k, sequence_length):
    """ Calculates mean nucleotide diversity, pi, using equation from Box 1.3
    in Charlesworth and Charleswoth (2010):
    (tot_diff / (k choose 2)) / sequence_length
    where:
        tot_diff = the total number of differences observed in all pairwise
            comparisons
        k = the number of chromosomes sampled
            (e.g. k = 6 for 3 diploid individuals)
        sequence_length = the number of bases analyzed per sequence
            (should be same for all sequences)
    """
    if k == 0:
        return 0
    elif k == 1:
        return 0
    else:
        numerator = float(tot_diff) / ((float(k) * (float(k) - 1)) / 2.0)
        return numerator / float(sequence_length)


def pi_site(allele_count_list):
    """Function calculates pi from Hohenlohe et al. (2010)

    pi = 1 - sum((n_i choose 2) / (n choose 2))
        where:
            n_i is count of allele i in sample
            n is the sum of n_i (allele1 plus allele2 in this case, as we
                assume bi-allelic sites

    inputs:
        allele_count_list is a list of the counts of the different alleles
            present at a site
    assumptions:
        snps
        sympy installed and imported for binomial coefficient calculation
    """
    n_tot = 0
    for i in allele_count_list:
        n_tot += i
    if n_tot == 0:
        return 0
    elif n_tot == 1:
        return 0
    else:
        pi_denom = float(sympy.binomial(n_tot, 2))
        pi_numer = 0.0
        for i in allele_count_list:
            pi_numer += float(sympy.binomial(i, 2))
        return (1.0 - (pi_numer / pi_denom))


def count_diffs(allele_list):
    """ Takes a list or string of alleles to count the number of differences
        among chromosomes.
    Example: For an input site from 3 diploid individuals with genotypes
        (A/T), (T/T), and (C/A), appropriate inputs would be either
        ["A","T","T","T","C","A"] or "ATTTCA".

    Returns:
        The number of pairwise differences as an integer
    """
    diffs = 0
    for index, i in enumerate(allele_list):
        for j in allele_list[index + 1:]:
            if i != j:
                diffs += 1
    return diffs


def rand_sample(data_input, n_sites):
    """Outputs a numpy array of resampled (with replacement) values
        for bootstrap function
    Inputs:
        data_input is the list of values to resample
        n_vals is the number of values to be resampled
    Requirements:
        numpy installed and imported

    Note: IF THE NUMBER OF SITES TO BE RESAMPLED IS LONGER THAN THE DATA_INPUT,
        ALL ADDITIONAL SITES ARE ASSUMED TO HAVE A VALUE OF 0
    """
    dim = len(data_input)
    # create random indices to grab the random sample
    indices = np.random.random_integers(0, n_sites, n_sites)
    # return all random values with indices less than the length of the input
    # table --- all other sites are assumed to be zero
    return np.take(data_input, indices[indices < dim])


def bootstrap_pi_distribution(data_input, n_sites, replicates):
    """ Returns a bootstrapped distribution (with replacement) as a list.
    data_input is the list of values
    n_sites is the total number of callable site (will be longer than
        data_input if the input VCF only contained polymorphic sites)
    replicates is the number of bootstrap replicates
    """
    resamples = []
    n_sites = float(n_sites)
    for i in range(replicates):
        resamples.append(np.sum(rand_sample(data_input, n_sites)) / n_sites)
    return resamples

###########################################################################
# Parse VCF file, calculate pi per site, and count pairwise differences ##
###########################################################################


print("Beginning diversity calculations")
counter = 0
for record in vcf_reader:
    if record.CHROM == args.chrom_inc and not args.haploid:
        for pop in populations:
            allele_list = []
            for indv in pop[0]:
                call = record.genotype(indv)
                if call['GT'] is not None and call['GT'] is not '.' \
                        and record.is_snp is True:
                    # call.gt_bases returns in the format "A/T",
                    # so this grabs the A and
                    #        the T, while skipping the / (or |)
                    if len(call.gt_bases) < 4:
                        allele_list.append(call.gt_bases[0])
                        allele_list.append(call.gt_bases[2])
            # Process allele list and calculate pi and number of differences
            #  pop[2].append([record.CHROM,record.POS,\
            # pi_overall(count_diffs(allele_list), len(allele_list), 1.0)])
            # Code to calculate Hohenlohe et al pi instead
            allele_count = collections.Counter(allele_list)
            pop[2].append([record.CHROM, record.POS, pi_site(
                [allele_count[x] for x in allele_count])])

    if record.CHROM == args.chrom_inc and args.haploid:
        for pop in populations:
            allele_list = []
            for indv in pop[0]:
                call = record.genotype(indv)
                if call['GT'] is not None and call['GT'] is not '.' \
                        and record.is_snp is True:
                    # call.gt_bases returns in the format "A"
                    # when processing haploid samples
                    # this grabs just "A" if only one allele
                    # and grabs both alleles if the format is "A/T"
                    if len(call.gt_bases) == 1:
                        allele_list.append(call.gt_bases)
                    elif len(call.gt_bases) < 4:
                        allele_list.append(call.gt_bases[0])
                        allele_list.append(call.gt_bases[2])
            # Process allele list and calculate pi and number of differences
            #  pop[2].append([record.CHROM,record.POS,\
            # pi_overall(count_diffs(allele_list), len(allele_list), 1.0)])
            # Code to calculate Hohenlohe et al pi instead
            allele_count = collections.Counter(allele_list)
            pop[2].append([record.CHROM, record.POS, pi_site(
                [allele_count[x] for x in allele_count])])

    counter += 1
    if counter % 10000 == 0:
        print("{} records complete...".format(counter))

print("VCF traversal complete")

#############################################
#  Print Output Files - one per population ##
#############################################

for pop in populations:
    out_file = args.out_directory + path.basename(pop[1]) + "_chr" + \
        args.chrom_inc + "_pi_output_by_site.txt"
    with open(out_file, "w") as f:
        w = csv.writer(f, dialect="excel-tab")
        w.writerows(pop[2])
