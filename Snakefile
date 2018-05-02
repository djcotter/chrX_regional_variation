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
CHR = ['chrX', 'chrY', 'chr8']

# Rules -----------------------------------------------------------------------

rule all:
    input:
        expand('02_diversity_by_site/results/{pops}_{chr_group}' +
               '_pi_output_by_site.txt', pops=POP_CODES + SUBPOP_CODES,
               chr_group=['individuals_chrX', 'individuals_chr8',
                          'males_chrY']),
        expand('03_filters/results/complete_{chr}_{filter_iter}.bed',
               chr=CHR, filter_iter=FILTER)

rule parse_populations:
    input:
        panel = config['panel']
    params:
        out_dir = '01_populations/results/',
        pop_parse = '01_populations/scripts/population_parser.py'
    output:
        out_pops = expand('01_populations/results/{pops}_{group}',
                          pops=POP_CODES + SUBPOP_CODES,
                          group=['males', 'females', 'individuals'])
    shell:
        'python {params.pop_parse} {input.panel} {params.out_dir}'

rule calculate_pi_chrX:
    input:
        individuals = expand('01_populations/results/{pops}_individuals',
                             pops=POP_CODES + SUBPOP_CODES),
        chrX = config['chromosomes']['chrX']
    params:
        calc_pi = "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_" + \
            "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py",
        out_dir = '02_diversity_by_site/results/'
    output:
        expand(path.join('02_diversity_by_site/results',
                         '{pops}_individuals_chrX_pi_output_by_site.txt'),
               pops=POP_CODES + SUBPOP_CODES)
    shell:
        "python {params.calc_pi} --vcf {input.chrX} "
        "--population_lists {input.individuals} --chrom_inc X "
        "--haploid --out_directory {params.out_dir}"

rule calculate_pi_chr8:
    input:
        individuals = expand('01_populations/results/{pops}_individuals',
                             pops=POP_CODES + SUBPOP_CODES),
        chr8 = config['chromosomes']['chr8']
    params:
        calc_pi = "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_" + \
            "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py",
        out_dir = '02_diversity_by_site/results/'
    output:
        expand(path.join('02_diversity_by_site/results',
                         '{pops}_individuals_chr8_pi_output_by_site.txt'),
               pops=POP_CODES + SUBPOP_CODES)
    shell:
        "python {params.calc_pi} --vcf {input.chr8} "
        "--population_lists {input.individuals} --chrom_inc 8 "
        "--out_directory {params.out_dir}"

rule calculate_pi_chrY:
    input:
        males_only = expand('01_populations/results/{pops}_males',
                            pops=POP_CODES + SUBPOP_CODES),
        chrY = config['chromosomes']['chrY']
    params:
        calc_pi = "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_" + \
            "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py",
        out_dir = '02_diversity_by_site/results/'
    output:
        expand(path.join('02_diversity_by_site/results',
                         '{pops}_males_chrY_pi_output_by_site.txt'),
               pops=POP_CODES + SUBPOP_CODES)
    shell:
        "python {params.calc_pi} --vcf {input.chrY} "
        "--population_lists {input.males_only} --chrom_inc Y "
        "--haploid --out_directory {params.out_dir}"

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

rule callable_sites_to_BED:
    input:
        callable = config['callable_sites']
    output:
        temp(callable_sites_{chr}.bed)
    shell:
        ""

rule apply_filter:
    input:
        filter = "03_filters/results/complete_{chr}_{filter_iter}.bed",
        callable_sites = "",
        pi_per_site = path.join("02_diversity_by_site/results/",
                                "{pops}_{group}_{chr}_pi_output_by_site.txt")
    output:
        path.join("04_window_analysis/results/",
                  "{pops}_{chr}_{filter_iter}_{window}_diversity.bed")

rule merge_files:
    input:
        rules.calculate_pi_chrX.output,
        rules.calculate_pi_chr8.output,
        rules.calculate_pi_chrY.output
    output:
        out_file = "01_populations/results/results2.html"
    shell:
        "touch 01_populations/results/results2.html"
