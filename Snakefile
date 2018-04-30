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
        "01_populations/results/results2.html"

rule parse_populations:
    input:
        panel = "data/integrated_call_samples_v3.20130502.ALL.panel"
    params:
        out_dir = "01_populations/results/",
        pop_parse = "01_populations/scripts/population_parser.py"
    output:
        out_pops = expand('{out_dir}/{pops}_{group}',
                          pops=SUBPOP_CODES + POP_CODES,
                          group=["males", "females", "individuals"],
                          out_dir='01_populations/results'),
        pop_table = '01_populations/results/pop_table.txt',
        subpop_table = '01_populations/results/subpop_table.txt'
    shell:
        'python {params.pop_parse} {input.panel} {params.out_dir}'

rule calculate_pi:
    input:
        individuals = expand('01_populations/results/{pops}_{group}',
                             pops=SUBPOP_CODES + POP_CODES,
                             group="individuals"),
        males_only = expand('01_populations/results/{pops}_{group}',
                            pops=SUBPOP_CODES + POP_CODES,
                            group="males")
    params:
        calc_pi = "02_diversity_by_site/scripts/Diversity_from_VCF_cyvcf_" + \
            "Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py",
        chrX = config['chromosomes']['chrX'],
        chrY = config['chromosomes']['chrY'],
        chr8 = config['chromosomes']['chr8'],
        out_dir = '02_diversity_by_site/results/',
        individuals = " ".join(expand('{out_dir}/{pops}_{group}',
                                      pops=SUBPOP_CODES + POP_CODES,
                                      group="individuals",
                                      out_dir='01_populations/results')),
        males_only = " ".join(expand('{out_dir}/{pops}_{group}',
                                     pops=SUBPOP_CODES + POP_CODES,
                                     group="males",
                                     out_dir='01_populations/results'))
    output:
        expand('{out_dir}/{pops}_individuals_chr{chr}_pi_output_by_site.txt',
               pops=POP_CODES + SUBPOP_CODES,
               chr=['X', '8'],
               out_dir='02_diversity_by_site/results'),
        expand('{out_dir}/{pops}_males_chrY_pi_output_by_site.txt',
               pops=POP_CODES + SUBPOP_CODES,
               chr='Y',
               out_dir='02_diversity_by_site/results')
    run:
        commands = ["python {params.calc_pi} --vcf {params.chrX} " +
                    "--population_lists {params.individuals} --chrom_inc " +
                    "X --haploid --out_directory {params.out_dir}",
                    "python {params.calc_pi} --vcf {params.chrY} " +
                    "--population_lists {params.males_only} --chrom_inc " +
                    "Y --haploid --out_directory {params.out_dir}",
                    "python {params.calc_pi} --vcf {params.chr8} " +
                    "--population_lists {params.individuals} --chrom_inc " +
                    "8 --out_directory {params.out_dir}"]
        for c in commands:
            shell(c)

rule merge_files:
    input:
        rules.calculate_pi.output
    output:
        out_file = "01_populations/results/results2.html"
    shell:
        "touch 01_populations/results/results2.html"
