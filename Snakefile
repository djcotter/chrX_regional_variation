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

# Global configurations -------------------------------------------------------

configfile: 'config.yml'

POP_CODES = sorted(json.load(open(confif['POPS_CODES']))["Populations"])
SUBPOP_CODES = sorted(json.load(open(confif['POPS_CODES']))["Subpopulations"])

rule all:
    input:
        "results.html"

rule parse_populations:
    input:
        panel = "data/integrated_call_samples_v3.20130502.ALL.panel"
    params:
        out_dir = "01_populations/results"
    output:
        out_subpops = expand(join(params.out_dir, 'populations/', {pops}, '_',
                                  {group}, '.txt'),
                             pops=SUBPOP_CODES,
                             group=["_males", "_females, _individuals"])
        out_pops = expand(join(params.out_dir, 'populations/', {pops}, '_',
                               {group}, '.txt'),
                          pops=POP_CODES,
                          group=["_males", "_females, _individuals"])
        join(params.out_dir, 'pop_table.txt')
        join(params.out_dir, 'subpop_table.txt')

rule merge_files:
    input:
        expand(join(params.out_dir, 'populations/', {pops}, '_',
                    {group}, '.txt'),
               pops=POP_CODES,
               group=["_males", "_females, _individuals"])
    ouput:
        "results.html"
    shell:
        "touch results.html"
