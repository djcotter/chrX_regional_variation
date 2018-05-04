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
POPULATIONS = sorted(json.load(open(config['POP_CODES']))['Populations'])
SUBPOPULATIONS = sorted(json.load(open(config['POP_CODES']))['Subpopulations'])

# select the filter from the configfile that should be used
FILTER = ['filter1']
WINDOW = ['100kb']

# sets the populations to be a list of all pops and subpops
# POPS = POPULATIONS + SUBPOPULATIONS
POPS = 'ALL'

# select a sex to use for analysis of chrX and chr8
# use "males", "females", or "individuals" (for both)
SEX = 'individuals'

# Global variables ------------------------------------------------------------

# defines the chromosomes to be analyzed
CHR = ['chrX', 'chr8', 'chrY']  # script not built for additional chromosomes

# link to diversity script
DIVERSITY_SCRIPT = '02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_' + \
    'Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py'

# takes the provided SEX and combines it with chromosomes to generate group_chr
GROUP = [SEX, SEX, 'males']
GROUP_CHR = [x + '-' + y for x, y in zip(GROUP, CHR)]

# Rules -----------------------------------------------------------------------

rule all:
    input:
        expand('02_diversity_by_site/results/{pops}_{group_chr}' +
               '_pi_output_by_site.txt', pops=POPS,
               group_chr=GROUP_CHR),
        expand('03_filters/results/complete_{chr}_{filter_iter}.bed',
               chr=CHR, filter_iter=FILTER),
        '01_populations/results/results2.html',
        expand('04_window_analysis/results/' +
               '{pops}_{group_chr}_{filter_iter}_{window}_diversity.bed',
               pops=POPS,
               group_chr=GROUP_CHR,
               filter_iter=FILTER, window=WINDOW),
        expand("04_window_analysis/inputs/{chr}_{window}_window.bed",
               chr=CHR, window=WINDOW)

rule parse_populations:
    input:
        panel = config['panel']
    params:
        out_dir = '01_populations/results/',
        pop_parse = '01_populations/scripts/population_parser.py'
    output:
        out_pops = expand('01_populations/results/{pops}_{group}',
                          pops=POPS,
                          group=['males', 'females', 'individuals'])
    shell:
        'python {params.pop_parse} {input.panel} {params.out_dir}'

rule calculate_pi_chrX:
    input:
        groups = expand('01_populations/results/{pops}_{group}',
                        pops=POPS, group=SEX),
        chrX = config['chromosomes']['chrX']
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        expand(path.join('02_diversity_by_site/results',
                         '{pops}_' + SEX + '-chrX_pi_output_by_site.txt'),
               pops=POPS)
    shell:
        "python {params.calc_pi} --vcf {input.chrX} "
        "--population_lists {input.groups} --chrom_inc X "
        "--haploid --out_directory {params.out_dir}"

rule calculate_pi_chr8:
    input:
        groups = expand('01_populations/results/{pops}_{group}',
                        pops=POPS, group=SEX),
        chr8 = config['chromosomes']['chr8']
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        expand(path.join('02_diversity_by_site/results',
                         '{pops}_' + SEX + '-chr8_pi_output_by_site.txt'),
               pops=POPS)
    shell:
        "python {params.calc_pi} --vcf {input.chr8} "
        "--population_lists {input.groups} --chrom_inc 8 "
        "--out_directory {params.out_dir}"

rule calculate_pi_chrY:
    input:
        males_only = expand('01_populations/results/{pops}_{group}',
                            pops=POPS, group='males'),
        chrY = config['chromosomes']['chrY']
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        expand(path.join('02_diversity_by_site/results',
                         '{pops}_males-chrY_pi_output_by_site.txt'),
               pops=POPS)
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

rule split_callable_sites:
    input:
        config['callable_sites']
    params:
        chrom = lambda wildcards: "\"" + wildcards.chr + "\""
    output:
        temp('callable_sites_{chr}.bed')
    shell:
        "cat {input} | awk \'$1 == {params.chrom} {{ print }}\' "
        "> {output}"

rule create_windows:
    input:
        "data/GRCh37_chromosome_coordinates.bed"
    params:
        win = lambda wildcards:
            config["windows"][wildcards.window]["win_size"],
        slide = lambda wildcards: "" if \
            config["windows"][wildcards.window]["overlap"] is False else \
            (" -s " + str(config["windows"][wildcards.window]["slide_size"])),
        chrom = lambda wildcards: "\"" + wildcards.chr + "\""
    output:
        "04_window_analysis/inputs/{chr}_{window}_window.bed"
    shell:
        "cat {input} | awk \'$1 == {params.chrom} {{ print }}\' | "
        "bedtools -b stdin -w {params.win}{params.slide} "
        "> {output}"

rule apply_filter:
    input:
        filter = expand('03_filters/results/complete_{chr}_{filter_iter}.bed',
                        chr=CHR, filter_iter=FILTER),
        callable_sites = expand('callable_sites_{chr}.bed',
                                chr=CHR),
        chrX_pi = rules.calculate_pi_chrX.output,
        chr8_pi = rules.calculate_pi_chr8.output,
        chrY_pi = rules.calculate_pi_chrY.output
    output:
        ""
    shell:
        ""

rule window_analysis:
    input:
        ""
    output:
        path.join('04_window_analysis/results/',
                  '{pop}_{group_chr}_{filter_iter}_{window}_diversity.bed')
    shell:
        ""

rule merge_files:
    input:
        rules.calculate_pi_chrX.output,
        rules.calculate_pi_chr8.output,
        rules.calculate_pi_chrY.output
    output:
        out_file = '01_populations/results/results2.html'
    shell:
        'touch 01_populations/results/results2.html'
