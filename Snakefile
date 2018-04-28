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
from os.path import join

# Global configurations -------------------------------------------------------

configfile: 'config.yml'

POP_CODES = sorted(json.load(open(config['POPS_CODES']))["Populations"])
SUBPOP_CODES = sorted(json.load(open(config['POPS_CODES']))["Subpopulations"])

# Rules -----------------------------------------------------------------------

rule all:
    input:
        "01_populations/results/results.html"

rule parse_populations:
    input:
        panel = "data/integrated_call_samples_v3.20130502.ALL.panel"
    params:
        out_dir = "01_populations/results/",
        pop_parse = "01_populations/scripts/population_parser.py"
    output:
        out_subs = expand('{out_dir}/{pops}_{group}.txt',
                          pops=SUBPOP_CODES,
                          group=["males", "females", "individuals"],
                          out_dir='01_populations/results/subpopulations'),
        out_pops = expand('{out_dir}/{pops}_{group}.txt',
                          pops=POP_CODES,
                          group=["males", "females", "individuals"],
                          out_dir='01_populations/results/populations'),
        pop_table = '01_populations/results/pop_table.txt',
        subpop_table = '01_populations/results/subpop_table.txt'
    shell:
        'python {params.pop_parse} {input.panel} {params.out_dir}'

rule calculate_pi:
    input:
        individuals = expand('{out_dir}/{pops}_{group}.txt',
                             pops=SUBPOP_CODES + POP_CODES,
                             group="individuals",
                             out_dir='01_populations/results/subpopulations'),
        males_only = expand('{out_dir}/{pops}_{group}.txt',
                            pops=SUBPOP_CODES + POP_CODES,
                            group="males",
                            out_dir='01_populations/results/populations')
    params:
        calc_pi = "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_" + \
            "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py",
        chrX = config['chromsomes']['chrX'],
        chrY = config['chromsomes']['chrY'],
        chr8 = config['chromsomes']['chr8'],
        out_dir = '02_diversity_by_site/results',
        individuals = (" ").join([x for x in input.individuals]),
        males_only = (" ").join([x for x in input.males_only])
    output:
        expand('{pop}_chr{chr}_diversity_by_site.txt',
               pops=POP_CODES + SUBPOP_CODES,
               chr=['X', 'Y', '8'])
    shell:
        """
        python {params.calc_pi} --vcf {params.chrX} --population_lists \
        {params.individuals} --chrom_inc chrX --out_directory {params.out_dir}

        python {params.calc_pi} --vcf {params.chrY} --population_lists \
        {params.males_only} --chrom_inc chrY --haploid --out_directory \
        {params.out_dir}

        python {params.calc_pi} --vcf {params.chr8} --population_lists \
        {params.indivduals} --chrom_inc chr8 --out_directory {params.out_dir}
        """

rule merge_files:
    input:
        expand('01_populations/results/subpopulations/{pops}_{group}.txt',
               pops=SUBPOP_CODES,
               group=["males", "females", "individuals"])
    output:
        out_file = "01_populations/results/results.html"
    shell:
        "touch 01_populations/results/results.html"
