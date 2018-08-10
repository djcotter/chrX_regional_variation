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
FILTER = ['filter1', 'filter4']
WINDOW = ['100kb']
LD_BIN = ['300kb']

# sets the populations to be a list of all pops and subpops
POPS = POPULATIONS + SUBPOPULATIONS
POPS = 'YRI'

# select a "sex" category to use for analysis of chrX and chr8
# use "males", "females", or "individuals" (for both)
SEX = 'individuals'

# select the pairwise substitution rates to use for divergence correction
CORRECTION = ['uncorrected',
              # 'rheMac2-hg19-corrected',
              # 'calJac3-hg19-corrected',
              'canFam3-hg19-corrected']


# Global variables ------------------------------------------------------------

# defines the chromosomes to be analyzed # chr9 not analyzed
CHR = ['chrX', 'chr8', 'chrY']  # script not built for more chromosomes

# link to diversity script
DIVERSITY_SCRIPT = '02_diversity_by_site/scripts/Diversity_from_VCF_pyvcf_' + \
    'Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py'

# takes the provided SEX and combines it with chromosomes to generate group_chr
# this is only in order to keep group and chr associated in rule all
GROUP = [SEX, SEX, 'males']
GROUP_CHR = [x + '_' + y for x, y in zip(GROUP, CHR)]

# Array containing the number of Mb for chrX ld 06_figures
PLOT_LENGTH = ['15', '156']

# Constrain wildcards globally

wildcard_constraints:
    pop = '[a-zA-Z]{3}'

# Rules -----------------------------------------------------------------------

# Global Rules ----------------------------------------------------------------
rule all:
    input:
        # chrX analyzed by region for all pops
        expand('04_window_analysis/results/' +
               '{pops}_{group_chr}_{filter_iter}_byRegion_{correction}' +
               '_diversity_wPvals.bed',
               pops=POPS, group_chr=SEX + '_chrX', filter_iter=FILTER,
               correction=CORRECTION),
        # windoweded graphs of diversity results
        expand('06_figures/results/' +
               '{pops}_{group_chr}_{filter_iter}_{window}_{correction}' +
               '_diversity.png',
               pops=POPS, group_chr=GROUP_CHR,
               filter_iter=FILTER, window=WINDOW,
               correction=CORRECTION),
        # windowed graphs of males and females across chrX
        expand('06_figures/results/' +
               '{pops}_chrX_malesAndFemales_{filter_iter}_{window}_' +
               '{correction}_diversity.png',
               pops=POPS, filter_iter=FILTER, window=WINDOW,
               correction=CORRECTION),
        # windowed graphs of diversity across the PAB
        # expand('06_figures/results/' +
        #        '{pops}_PAB_{filter_iter}_{window}_{correction}_diversity.png',
        #        pops=POPS, filter_iter=FILTER, window=WINDOW,
        #        correction=CORRECTION),
        # output for ld_window_analysis
        expand('06_figures/results/' +
               '{pops}_{group_chr}_{window}_windows_{ld_bin}_LDbins_' +
               '95bootstrapCI_{plotSize}Mb.png',
               pops=POPS, group_chr="chrX_females",
               window=WINDOW, ld_bin=LD_BIN, plotSize=PLOT_LENGTH),
        # output for diversity split by chr/region
        # expand('06_figures/results/{pop}_{group}_totalDiversity_' +
        #        '{filter_iter}_{correction}_byChrRegion.png',
        #        pop=POPS, filter_iter=FILTER, group=SEX,
        #        correction=CORRECTION),
        # output for ratios tables
        expand('06_figures/results/' +
               '{pop}_{group}_{filter_iter}_{correction}_ratios.txt',
               pop=POPS, group=SEX, filter_iter=FILTER, correction=CORRECTION),
        # ld vs pi correlation plots
        expand('06_figures/results/' +
               '{pops}_{group_chr}_{window}_windows_{correction}' +
               '_{filter_iter}_{ld_bin}_LDbin_correlation.png',
               pops=POPS, group_chr="chrX_females", correction=CORRECTION,
               window=WINDOW, ld_bin=LD_BIN, filter_iter=FILTER),
        # diversity ratios figure
        expand('06_figures/results/' +
               'subpops_X-A_PAR-A_XTR-A_ratios_{correction}_{filter_iter}.png',
               correction=CORRECTION, filter_iter=FILTER),
        # demography corrected PARvX figure
        expand('06_figures/results/' +
               'subpops_X-PAR_demography_corrected_ratios_{correction}' +
               '_{filter_iter}.png', correction=CORRECTION,
               filter_iter=FILTER)

rule parse_populations:
    input:
        panel = config['panel']
    params:
        out_dir = '01_populations/results/',
        pop_parse = '01_populations/scripts/population_parser.py'
    output:
        out_pops = expand('01_populations/results/{pops}_{group}',
                          pops=POPULATIONS + SUBPOPULATIONS,
                          group=['males', 'females', 'individuals'])
    shell:
        'python {params.pop_parse} {input.panel} {params.out_dir}'

rule subset_VCF:
    input:
        pop_file = '01_populations/results/{pop}_{group}',
        vcf_file = lambda wildcards: config['chromosomes'][wildcards.chr]
    output:
        temp(path.join('data', 'subset_{chr}_{pop}_{group}.vcf'))
    shell:
        "bcftools view -S {input.pop_file} {input.vcf_file} > {output}"

rule calculate_pi_autosomes:
    input:
        group = '01_populations/results/{pop}_{group}',
        vcf = path.join('data', 'subset_{chr}_{pop}_{group}.vcf')
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/',
        chrom = lambda wildcards: wildcards.chr[3:]
    wildcard_constraints:
        chr = 'chr[0-9]+'
    output:
        path.join('02_diversity_by_site/results',
                  '{pop}_{group}_{chr}_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc {params.chrom} "
        "--out_directory {params.out_dir}"

# rule calculate_pi_chrX:
#     input:
#         group = '01_populations/results/{pop}_{group}',
#         vcf = path.join('data', 'subset_chrX_{pop}_{group}.vcf')
#     params:
#         calc_pi = DIVERSITY_SCRIPT,
#         out_dir = '02_diversity_by_site/results/'
#     output:
#         path.join('02_diversity_by_site/results',
#                   '{pop}_{group}_chrX_pi_output_by_site.txt')
#     shell:
#         "python {params.calc_pi} --vcf {input.vcf} "
#         "--population_lists {input.group} --chrom_inc X "
#         "--haploid --out_directory {params.out_dir}"

rule subset_AND_filter_TEST:
    input:
        path.join('data', 'subset_{chr}_{pop}_{group}.vcf')
    output:
        path.join('data', 'filtered_vcf_{chr}_{pop}_{group}.vcf')
    shell:
        "bcftools view -Ou -m2 -M2 -v snps {input} | bcftools view -Ov "
        "--min-ac 1:minor > {output}"

rule calculate_pi_chrX_TEST1:
    input:
        group = '01_populations/results/{pop}_{group}',
        vcf = lambda wildcards: config['chromosomes']['chrX']
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        path.join('02_diversity_by_site/results',
                  '{pop}_{group}_chrX_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc X "
        "--haploid --out_directory {params.out_dir}"

rule calculate_pi_chrY:
    input:
        group = '01_populations/results/{pop}_{group}',
        vcf = path.join('data', 'subset_chrY_{pop}_{group}.vcf')
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        path.join('02_diversity_by_site/results',
                  '{pop}_{group}_chrY_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc Y "
        "--haploid --out_directory {params.out_dir}"

rule create_filter:
    input:
        lambda wildcards: expand(
            '03_filters/raw_filters/{filter_type}',
            filter_type=config["filters"][wildcards.filter_iter])
    params:
        chrom = lambda wildcards: wildcards.chr
    output:
        '03_filters/results/complete_{chr}_{filter_iter}.bed'
    shell:
        "cat {input} | grep {params.chrom} | sort -k1,1 -k2,2n | "
        "cut -f1,2,3 | bedtools merge -i stdin > {output}"

rule split_callable_sites:
    input:
        config['callable_sites']
    params:
        chrom = lambda wildcards: "\"" + wildcards.chr + "\""
    output:
        temp('data/callable_sites_{chr}.bed')
    shell:
        "grep {params.chrom} {input} > {output}"

rule create_windows:
    input:
        'data/hg19_chromosome_coordinates.bed'
    params:
        win = lambda wildcards:
            config["windows"][wildcards.window]["win_size"],
        slide = lambda wildcards: "" if \
            config["windows"][wildcards.window]["overlap"] is False else \
            (" -s " + str(config["windows"][wildcards.window]["slide_size"])),
        chrom = lambda wildcards: "\"" + wildcards.chr + "\""
    output:
        '04_window_analysis/inputs/{chr}_{window}_window.bed'
    shell:
        "grep {params.chrom} {input} | "
        "bedtools makewindows -b stdin -w {params.win}{params.slide} "
        "> {output}"

rule convert_diverstiy_to_bed:
    input:
        path.join('02_diversity_by_site', 'results',
                  '{pop}_{group}_{chr}' +
                  '_pi_output_by_site.txt')
    params:
        bedConvert = path.join('02_diversity_by_site', 'scripts',
                               'bedConvert.py')
    output:
        temp(path.join('02_diversity_by_site', 'results',
                       '{pop}_{group}_{chr}' +
                       '_pi_output_by_site.bed'))
    shell:
        "python {params.bedConvert} {input} {output}"

rule filter_callable_sites:
    input:
        filter = '03_filters/results/complete_{chr}_{filter_iter}.bed',
        callable_sites = 'data/callable_sites_{chr}.bed'
    output:
        path.join('04_window_analysis', 'inputs',
                  'callable_sites_{chr}_{filter_iter}.bed')
    shell:
        "bedtools subtract -a {input.callable_sites} -b {input.filter} "
        "> {output}"

rule filter_diversity_by_site:
    input:
        diversity_by_site = path.join('02_diversity_by_site', 'results',
                                      '{pop}_{group}_{chr}' +
                                      '_pi_output_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_{filter_iter}.bed')
    output:
        path.join('04_window_analysis', 'inputs',
                  '{pop}_{group}_{chr}_{filter_iter}' +
                  '_pi_by_site.bed')
    shell:
        "bedtools intersect -a {input.diversity_by_site} "
        "-b {input.filtered_callable} > {output}"

# Window Analysis -------------------------------------------------------------
rule window_analysis:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{pop}_{group}_{chr}_{filter_iter}' +
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_' +
                                      '{filter_iter}.bed'),
        windows = path.join('04_window_analysis', 'inputs',
                            '{chr}_{window}_window.bed')
    wildcard_constraints:
        # this regular expression matches things like '100kb' or '1Mb'
        # it is used as a way to allow all other window wildcards other than
        # 'byRegion' to default to this rule
        window = '[0-9]+[A-Za-z]+'
    params:
        window_calcs = '04_window_analysis/scripts/window_calculations.py',
        slide = lambda wildcards: '' if \
            config["windows"][wildcards.window]["overlap"] is False else \
            '--sliding ',
        replicates = 1000
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                  '_diversity.unfiltered.bed')
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --windows {input.windows} "
        "{params.slide}--replicates {params.replicates} --output {output}"

rule window_analysis_byRegion:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{pop}_{group}_{chr}_{filter_iter}' +
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_' +
                                      '{filter_iter}.bed')
    wildcard_constraints:
        window = 'byRegion'
    params:
        window_calcs = '04_window_analysis/scripts/window_calculations.py',
        replicates = 1000
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                  '_diversity.bed')
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --chrX_windows "
        "--replicates {params.replicates} --output {output}"

rule filter_windows_by_callable_sites:
    input:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                  '_diversity.unfiltered.bed')
    params:
        script = path.join('04_window_analysis', 'scripts',
                           'filter_windows_byCallableSites.py'),
        winSize = lambda wildcards:
            config["windows"][wildcards.window]["win_size"]
    wildcard_constraints:
        # this regular expression matches things like '100kb' or '1Mb'
        # it is used as a way to allow all other window wildcards other than
        # 'byRegion' to default to this rule
        window = '[0-9]+[A-Za-z]+'
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                  '_diversity.bed')
    shell:
        "python {params.script} --input {input} --windowSize {params.winSize} "
        "--filter 0.1 --output {output}"

rule create_pseudo_uncorrected_divergence:
    output:
        temp(path.join('data', 'substitution_rates',
                       '{chr}_uncorrected_{filter_iter}_{window}' +
                       '_substitution_rates.txt'))
    shell:
        "touch {output}"

rule divergence_JukesCantor69_correction:
    input:
        path.join('data', 'substitution_rates',
                  '{chr}_{correction}_{filter_iter}_{window}' +
                  '_substitution_rates.txt')
    params:
        script = path.join('data', 'substitution_rates', 'scripts',
                           'Calculate_JC69_from_Galaxy_output.py'),
        corrected = lambda wildcards: False if wildcards.correction == \
            "uncorrected" else True
    output:
        temp(path.join('data', 'substitution_rates',
                       '{chr}_{correction}_{filter_iter}_{window}' +
                       '_substitution_rates_JC69.txt'))
    run:
        if params.corrected is True:
            shell("python {params.script} --input_file {input} " +
                  "--output_file {output}")
        else:
            shell("mv {input} {output}")

rule windowed_divergence_correction:
    input:
        diversity = path.join('04_window_analysis', 'results',
                              '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                              '_diversity.bed'),
        divergence = path.join('data', 'substitution_rates',
                               '{chr}_{correction}_{filter_iter}_{window}' +
                               '_substitution_rates_JC69.txt')
    params:
        corrected = lambda wildcards: False if wildcards.correction == \
            "uncorrected" else True,
        script = '04_window_analysis/scripts/divergence_correction.py'
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}_' +
                  '{correction}_diversity.bed')
    run:
        if params.corrected:
            shell("python {params.script} --diversity {input.diversity} " +
                  "--divergence {input.divergence} --output {output}")
        else:
            shell("cp {input.diversity} {output}")

rule permute_chrX_regions:
    input:
        byWindow_100kb = path.join('04_window_analysis', 'results',
                                   '{pop}_{group}_{chr}_{filter_iter}' +
                                   '_100kb_{correction}_diversity.bed'),
        byRegion = path.join('04_window_analysis', 'results',
                             '{pop}_{group}_{chr}_{filter_iter}' +
                             '_byRegion_{correction}_diversity.bed')
    params:
        permutation_script = path.join('04_window_analysis', 'scripts',
                                       'permute_chrX_windows.py'),
        replicates = 10000
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_byRegion_{correction}' +
                  '_diversity_wPvals.bed')
    shell:
        "python {params.permutation_script} --byRegion {input.byRegion} "
        "--byWindow {input.byWindow_100kb} --replicates {params.replicates} "
        "--output {output}"

rule plot_windowed_diversity:
    input:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}_{correction}' +
                  '_diversity.bed')
    params:
        R_script = path.join('06_figures', 'scripts',
                             'plot_windowed_diversity.R'),
        chrom = lambda wildcards: wildcards.chr,
        height = lambda wildcards: config[wildcards.correction][wildcards.chr]
    output:
        path.join('06_figures', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}_{correction}' +
                  '_diversity.png')
    shell:
        "Rscript {params.R_script} -i {input} -o {output} -c {params.chrom} "
        "--maxHeight {params.height}"

rule plot_sex_specific_chrX_windows:
    input:
        chrX_males = path.join('04_window_analysis', 'results',
                               '{pop}_males_chrX_{filter_iter}_{window}' +
                               '_{correction}_diversity.bed'),
        chrX_females = path.join('04_window_analysis', 'results',
                                 '{pop}_females_chrX_{filter_iter}_{window}' +
                                 '_{correction}_diversity.bed')
    params:
        R_script = path.join('06_figures', 'scripts',
                             'plot_chrX_maleFemale_diversity.R'),
        height = lambda wildcards: config[wildcards.correction]['chrX']
    output:
        path.join('06_figures', 'results',
                  '{pop}_chrX_malesAndFemales_{filter_iter}_{window}' +
                  '_{correction}_diversity.png')
    shell:
        "Rscript {params.R_script} --males {input.chrX_males} --females "
        "{input.chrX_females} --maxHeight {params.height} -o {output}"

rule plot_PAB_diversity:
    input:
        chrX_males = path.join('04_window_analysis', 'results',
                               '{pop}_males_chrX_{filter_iter}_{window}' +
                               '_{correction}_diversity.bed'),
        chrX_females = path.join('04_window_analysis', 'results',
                                 '{pop}_females_chrX_{filter_iter}_{window}' +
                                 '_{correction}_diversity.bed'),
        chrY = path.join('04_window_analysis', 'results',
                         '{pop}_males_chrY_{filter_iter}_{window}' +
                         '_{correction}_diversity.bed')
    params:
        R_script = path.join('06_figures', 'scripts',
                             'plot_PAB_diversity.R'),
        height = lambda wildcards: config[wildcards.correction]['chrX']
    output:
        path.join('06_figures', 'results',
                  '{pop}_PAB_{filter_iter}_{window}_{correction}' +
                  '_diversity.png')
    shell:
        "Rscript {params.R_script} --chrX_females {input.chrX_females} "
        "--chrX_males {input.chrX_males} --chrY {input.chrY} "
        "--maxHeight {params.height} --output {output}"

# Analysis by Region ----------------------------------------------------------

rule get_wholeChr_bed:
    input:
        path.join('data', 'hg19_chromosome_coordinates.bed')
    params:
        chrom = lambda wildcards: wildcards.chr
    output:
        temp(path.join('data', '{chr}_wholeChr.bed'))
    shell:
        'grep {params.chrom} {input} > {output}'

rule window_analysis_wholeChr:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{pop}_{group}_{chr}_{filter_iter}' +
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_' +
                                      '{filter_iter}.bed'),
        windows = path.join('data', '{chr}_{window}.bed')
    wildcard_constraints:
        window = 'wholeChr'
    params:
        window_calcs = '04_window_analysis/scripts/window_calculations.py'
    output:
        temp(path.join('04_window_analysis', 'results',
                       '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                       '_diversity.bed'))
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --windows {input.windows} "
        "--output {output}"

# rule plot_diversity_byRegion_byWholeChr:
#     input:
#         chrX = path.join('04_window_analysis', 'results',
#                          '{pop}_{group}_chrX_{filter_iter}_byRegion' +
#                          '_{correction}_diversity.bed'),
#         chrY = path.join('04_window_analysis', 'results',
#                          '{pop}_males_chrY_{filter_iter}_wholeChr' +
#                          '_{correction}_diversity.bed'),
#         chr8 = path.join('04_window_analysis', 'results',
#                          '{pop}_{group}_chr8_{filter_iter}_wholeChr' +
#                          '_{correction}_diversity.bed')
#     params:
#         R_script = ''
#     output:
#         path.join('06_figures', 'results',
#                   '{pop}_{group}_totalDiversity_{filter_iter}_{correction}' +
#                   '_byChrRegion.png')
#     shell:
#         "touch {output}"

rule calculate_A_ratios:
    input:
        chrX = path.join('04_window_analysis', 'results',
                         '{pop}_{group}_chrX_{filter_iter}_byRegion' +
                         '_{correction}_diversity.bed'),
        chrY = path.join('04_window_analysis', 'results',
                         '{pop}_males_chrY_{filter_iter}_wholeChr' +
                         '_{correction}_diversity.bed'),
        chr8 = path.join('04_window_analysis', 'results',
                         '{pop}_{group}_chr8_{filter_iter}_wholeChr' +
                         '_{correction}_diversity.bed')
    params:
        script = '06_figures/scripts/ratios_table.py'
    output:
        path.join('06_figures', 'results',
                  '{pop}_{group}_{filter_iter}_{correction}_ratios.txt')
    shell:
        "python {params.script} --chrX {input.chrX} --chrY {input.chrY} "
        "--chr8 {input.chr8} --output {output}"

rule prepare_subpop_ratio_plotting_data:
    input:
        chrX = lambda wildcards: expand(path.join('04_window_analysis',
                                                  'results',
                                                  '{pops}_{Group}_chrX_' +
                                                  '{filter}_byRegion' +
                                                  '_{correction_div}_' +
                                                  'diversity.bed'),
                                        pops=SUBPOPULATIONS, Group=SEX,
                                        filter=wildcards.filter_iter,
                                        correction_div=wildcards.correction),
        chrY = lambda wildcards: expand(path.join('04_window_analysis',
                                                  'results',
                                                  '{pops}_males_chrY_' +
                                                  '{filter}_wholeChr' +
                                                  '_{correction_div}_' +
                                                  'diversity.bed'),
                                        pops=SUBPOPULATIONS,
                                        filter=wildcards.filter_iter,
                                        correction_div=wildcards.correction),
        chr8 = lambda wildcards: expand(path.join('04_window_analysis',
                                                  'results',
                                                  '{pops}_{Group}_chr8_' +
                                                  '{filter}_wholeChr' +
                                                  '_{correction_div}_' +
                                                  'diversity.bed'),
                                        pops=SUBPOPULATIONS, Group=SEX,
                                        filter=wildcards.filter_iter,
                                        correction_div=wildcards.correction)
    output:
        path.join('06_figures', 'results',
                  'subpops_{filter_iter}_{correction}_ratios_table.txt')
    params:
        script = path.join('06_figures', 'scripts',
                           'format_byRegion_data.py')
    shell:
        "python {params.script} --chrX_byRegion {input.chrX} "
        "--chr8_wholeChr {input.chr8} --chrY_wholeChr {input.chrY} "
        "--output {output}"

rule plot_A_Ratios_across_subpops:
    input:
        path.join('06_figures', 'results',
                  'subpops_{filter_iter}_{correction}_ratios_table.txt')
    output:
        path.join('06_figures', 'results',
                  'subpops_X-A_PAR-A_XTR-A_ratios_{correction}' +
                  '_{filter_iter}.png')
    params:
        Rscript = path.join('06_figures', 'scripts',
                            'plot_subpop_diversity_ratios.R')
    shell:
        "Rscript {params.Rscript} --subpops_data {input} -o {output}"

rule plot_relative_PARvX_ratios:
    input:
        path.join('06_figures', 'results',
                  'subpops_{filter_iter}_{correction}_ratios_table.txt')
    output:
        path.join('06_figures', 'results',
                  'subpops_X-PAR_demography_corrected_ratios_{correction}' +
                  '_{filter_iter}.png')
    params:
        Rscript = path.join('06_figures', 'scripts',
                            'relative_XA_ratios_bySubpop.R')
    shell:
        "Rscript {params.Rscript} --subpops_data {input} -o {output}"

# LD analysis -----------------------------------------------------------------
rule cythonize_ld_script:
    input:
        ld_analysis = '05_ld_windows/scripts/ld_analysis.pyx',
        setup = '05_ld_windows/scripts/setup.py'
    output:
        '05_ld_windows/scripts/ld_analysis.c'
    shell:
        "python {input.setup} build_ext --inplace"

rule subset_VCF_for_LD:
    input:
        pop_file = '01_populations/results/{pop}_{group}',
        vcf_file = lambda wildcards: config['chromosomes'][wildcards.chr]
    output:
        temp(path.join('data', 'subset_LD_{chr}_{pop}_{group}.vcf'))
    shell:
        "bcftools view -S {input.pop_file} {input.vcf_file} > {output}"


rule filter_vcf:
    input:
        vcf = path.join('data', 'subset_LD_{chr}_{pop}_{group}.vcf'),
        targets = path.join('03_filters', 'results',
                            'complete_{chr}_{filter_iter}.bed')
    output:
        temp(path.join('data', 'subset_LD_{chr}_{pop}_{group}_{filter_iter}' +
                       '_snpsONLY-mac-filtered.recode.vcf'))
    params:
        prefix = path.join('data', 'subset_LD_{chr}_{pop}_{group}_' +
                           '{filter_iter}_snpsONLY-mac-filtered'),
        filter = lambda wildcards: wildcards.filter_iter
    shadow: "full"
    shell:
        "cat {input.targets} | sed 's/^chr//' > temp_{params.filter}.bed && "
        "bcftools view -Ou -m2 -M2 -v snps {input.vcf} | bcftools view -Ov "
        "--min-ac 1:minor | vcftools --vcf - "
        "--exclude-bed temp_{params.filter}.bed --recode --out {params.prefix}"

rule calculate_ld:
    input:
        path.join('data', 'subset_LD_{chr}_{pop}_{group}_{filter_iter}' +
                  '_snpsONLY-mac-filtered.recode.vcf')
    params:
        out_path = path.join('data', '{pop}_{chr}_{group}_{filter_iter}' +
                             '_filtered_ld_R2')
    output:
        temp(path.join('data', '{pop}_{chr}_{group}_{filter_iter}' +
                       '_filtered_ld_R2.ld'))
    shadow: "full"
    threads: 4
    shell:
        "plink2 --vcf {input} --r2 with-freqs --memory 8000 "
        "--threads {threads} "
        "--ld-window 700000 --ld-window-kb 1100 --out {params.out_path}"

rule ld_window_analysis:
    input:
        LD = path.join('data', '{pop}_{chr}_{group}_{filter_iter}' +
                       '_filtered_ld_R2.ld'),
        windows = path.join('04_window_analysis', 'inputs',
                            '{chr}_{window}_window.bed'),
        script = path.join('05_ld_windows', 'scripts', 'ld_analysis.c')
    params:
        LD_bin = lambda wildcards: config['ld_bins'][wildcards.ld_bin]['size'],
        script = '05_ld_windows/scripts/average_ld_by_window.py'
    output:
        path.join('05_ld_windows', 'results', '{pop}_{chr}_{group}_' +
                  '{window}_windows_{filter_iter}' +
                  '_{ld_bin}_LDbins_95bootstrapCI.txt')
    shell:
        "python {params.script} --plink_ld {input.LD} "
        "--windows {input.windows} --binSize {params.LD_bin} "
        "--output {output}"

rule plot_ld_windows:
    input:
        path.join('05_ld_windows', 'results', '{pop}_{chr}_{group}_' +
                  '{window}_windows_{filter_iter}' +
                  '_{ld_bin}_LDbins_95bootstrapCI.txt')
    output:
        path.join('06_figures', 'results',
                  '{pop}_{chr}_{group}_{window}_windows_{filter_iter}' +
                  '_{ld_bin}_LDbins_95bootstrapCI_{LDplotSize}Mb.png')
    params:
        R_script = path.join('06_figures', 'scripts',
                             'plot_LD_bins.R'),
        winSize = lambda wildcards:
            config["windows"][wildcards.window]["win_size"],
        plot_size = lambda wildcards: wildcards.LDplotSize
    shell:
        "Rscript {params.R_script} -i {input} --winSize {params.winSize} "
        "--zoom {params.plot_size} -o {output} --maxHeight 1"

rule plot_ld_pi_correlation:
    input:
        ld = path.join('05_ld_windows', 'results', '{pop}_{chr}_{group}_' +
                       '{window}_windows_{filter_iter}' +
                       '_{ld_bin}_LDbins_95bootstrapCI.txt'),
        pi = path.join('04_window_analysis', 'results',
                       '{pop}_{group}_{chr}_{filter_iter}_{window}_' +
                       '{correction}_diversity.bed')
    params:
        R_script = path.join('06_figures', 'scripts',
                             'ld_pi_correlation.R')
    output:
        path.join('06_figures', 'results',
                  '{pop}_{chr}_{group}_{window}_windows_{correction}' +
                  '_{filter_iter}_{ld_bin}_LDbin_correlation.png')
    shell:
        "Rscript {params.R_script} --LD {input.ld} --diversity {input.pi} "
        "--output {output}"
