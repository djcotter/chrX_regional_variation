"""
PAB_variation.Snakefile
Daniel Cotter

analyze diversity and LD across chrX and chrY from 1000 genomes data
---------------------------------------------------------------------

Requires:
    conda
        to activate the included PAB_variation.yml environment
"""

# Import statements -----------------------------------------------------------

import json
from os import path

# Global configurations -------------------------------------------------------

# declare a path to the configuration file
configfile: 'config.yml'

# parse the populations from the provided .json file into arrays of pops codes
POP_CODES = sorted(json.load(open(config['POPS_CODES']))['Populations'])
SUBPOP_CODES = sorted(json.load(open(config['POPS_CODES']))['Subpopulations'])
# select the filter from the configfile that should be used
FILTER = ['filter1']
# link to diversity script
DIVERSITY_SCRIPT = ""

# Rules -----------------------------------------------------------------------

rule all:
    input:
        expand('02_diversity_by_site/results/{pops}_individuals_{chr}' +
               '_pi_output_by_site.txt', pops=POP_CODES + SUBPOP_CODES,
               chr=['chrX', 'chr8']),
        expand('02_diversity_by_site/results/{pops}_males_{chr}' +
               '_pi_output_by_site.txt', pops=POP_CODES + SUBPOP_CODES,
               chr=['chrY']),
        expand('03_filters/results/complete_{chr}_{filter_iter}.bed',
               chr=['chrX', 'chrY', 'chr8'],
               filter_iter=FILTER)

rule parse_populations:
    input:
        panel = 'data/integrated_call_samples_v3.20130502.ALL.panel'
    params:
        out_dir = '01_populations/results/',
        pop_parse = '01_populations/scripts/population_parser.py'
    output:
        out_pops = '01_populations/results/{pops}_{group}'
    shell:
        'python {params.pop_parse} {input.panel} {params.out_dir}'

rule calculate_pi:
    input:
        individuals = expand('01_populations/results/{pops}_individuals',
                             pops=POP_CODES + SUBPOP_CODES),
        males_only = expand('01_populations/results/{pops}_males',
                            pops=POP_CODES + SUBPOP_CODES),
        chrom_file = lambda wildcards: config['chromosomes'][wildcards.chr]
    params:
        calc_pi = path.join(
            "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_",
            "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py"),
        out_dir = '02_diversity_by_site/results/',
        chrom = lambda wildcards: wildcards.chr,
        opts = lambda wildcards: '--population_lists {{input.individuals}}' + \
            ' --chrom_inc X --haploid'if wildcards.chr == 'chrX' else \
            ('--population_lists {{input.males_only}}' +
             ' --chrom_inc Y --haploid' if wildcards.chr == 'chrY'
             else '--population_lists {{input.individuals}} --chrom_inc 8')
    output:
        path.join('02_diversity_by_site/results/{pops}_{group}_{chr}' +
                  '_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.chrom_file} {params.opts} "
        "--out_directory {params.out_dir}"
#     run:
#         commands = ["python {params.calc_pi} --vcf {input.chrX} " +
#                     "--population_lists {input.individuals} --chrom_inc " +
#                     "X --haploid --out_directory {params.out_dir}",
#
#                     "python {params.calc_pi} --vcf {input.chrY} " +
#                     "--population_lists {input.males_only} --chrom_inc " +
#                     "Y --haploid --out_directory {params.out_dir}",
#
#                     "python {params.calc_pi} --vcf {input.chr8} " +
#                     "--population_lists {input.individuals} --chrom_inc " +
#                     "8 --out_directory {params.out_dir}"]
#         if params.chrom == 'chrX':
#             shell(commands[0])
#         elif params.chrom == 'chrY':
#             shell(commands[1])
#         elif params.chrom == 'chr8':
#             shell(commands[2])
#
# rule calculate_pi_chrX:
#     input:
#         individuals = expand('01_populations/results/{pops}_individuals',
#                              pops=POP_CODES + SUBPOP_CODES),
#         chrX = config['chromsome']['chrX']
#     params:
#         calc_pi = path.join(
#             "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_",
#             "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py"),
#         out_dir = '02_diversity_by_site/results/'
#     output:
#         path.join('02_diversity_by_site/results/{pops}_{group}_{chr}',
#                   '_pi_output_by_site.txt')
#     shell:
#         "python {params.calc_pi} --vcf {input.chrX} "
#         "--population_lists {input.individuals} --chrom_inc X "
#         "--haploid --out_directory {params.out_dir}"
#
# rule calculate_pi_chrY:
#     input:
#         males_only = expand('01_populations/results/{pops}_males',
#                             pops=POP_CODES + SUBPOP_CODES),
#         chrY = config['chromsome']['chrY']
#     params:
#         calc_pi = path.join(
#             "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_",
#             "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py"),
#         out_dir = '02_diversity_by_site/results/'
#     output:
#         path.join('02_diversity_by_site/results/{pops}_{group}_{chr}',
#                   '_pi_output_by_site.txt')
#     shell:
#         "python {params.calc_pi} --vcf {input.chrY} "
#         "--population_lists {input.males_only} --chrom_inc Y "
#         "--haploid --out_directory {params.out_dir}"

rule create_filter:
    input:
        lambda wildcards: expand(
            '03_filters/raw_filters/{filter_type}_{chrom}_filter.bed',
            filter_type=config["filters"][wildcards.filter_iter],
            chrom=[wildcards.chr])
    output:
        '03_filters/results/complete_{chr}_{filter_iter}.bed'
    shell:
        "cat {input} | sort -k1,1 -k2,2n | "
        "awk \'BEGIN{{OFS=" "}}{{print $1,$2,$3,$4}}\' | "
        "bedtools merge -i stdin > {output}"

# rule callable_sites_to_BED:
#     input:
#         ""
#     output:
#         ""
#     shell:
#         ""

rule apply_filter:
    input:
        filter = "03_filters/results/complete_{chr}_{filter_iter}.bed",
        callable_sites = "",
        pi_per_site = path.join("02_diversity_by_site/results/",
                                "{pops}_{group}_{chr}_pi_output_by_site.txt")
    output:
        path.join("04_window_analysis/results/",
                  "{pops}_{group}_{chr}_{filter_iter}_{window}_diversity.bed")

rule merge_files:
    input:
        rules.calculate_pi.output
    output:
        out_file = "01_populations/results/results2.html"
    shell:
        "touch 01_populations/results/results2.html"
