"""
PAB_variation.Snakefile
Daniel Cotter

analyze diversity and LD across chrX and chrY from 1000 genomes data
---------------------------------------------------------------------

Requires:
    conda
        to activate the included PAB_variation.yml environment
"""

# import statements

import json
from os.path import join

# Global configurations -------------------------------------------------------

configfile: 'config.yml'

POP_CODES = sorted(json.load(open(config['POPS_CODES']))["Populations"])
SUBPOP_CODES = sorted(json.load(open(config['POPS_CODES']))["Subpopulations"])

rule all:
    input:
        "results.html"

rule parse_populations:
    input:
        panel = "data/integrated_call_samples_v3.20130502.ALL.panel"
    params:
        out_dir = "01_populations/results",
        pop_parse = "01_populations/scripts/population_parser.py"
    output:
        out_subs = expand('{out_dir}/subpopulations/{pops}_{group}.txt',
                          pops=SUBPOP_CODES,
                          group=["_males", "_females", "_individuals"],
                          out_dir={params.out_dir}),
        out_pops = expand('{out_dir}/populations/{pops}_{group}.txt',
                          pops=POP_CODES,
                          group=["_males", "_females", "_individuals"],
                          out_dir={params.out_dir}),
        pop_table = join({params.out_dir}, 'pop_table.txt'),
        subpop_table = join({params.out_dir}, 'subpop_table.txt')
    shell:
        " ".join(params.pop_parse, input.panel, params.out_dir)

rule merge_files:
    input:
        expand('01_populations/results/subpopulations/{pops}_{group}.txt',
               pops=SUBPOP_CODES,
               group=["_males", "_females, _individuals"])
    output:
        out_file = "results.html"
    shell:
        "touch results.html"
