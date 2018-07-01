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
LD_BIN = ['300kb']

# sets the populations to be a list of all pops and subpops
POPS = POPULATIONS + SUBPOPULATIONS
POPS = 'ALL'

# select a "sex" category to use for analysis of chrX and chr8
# use "males", "females", or "individuals" (for both)
SEX = 'individuals'

# select the pairwise substitution rates to use for divergence correction
CORRECTION = ['uncorrected', 'rheMac2-hg19-corrected']

# Global variables ------------------------------------------------------------

# defines the chromosomes to be analyzed
CHR = ['chrX', 'chr8', 'chrY']  # script not built for additional chromosomes

# link to diversity script
DIVERSITY_SCRIPT = '02_diversity_by_site/scripts/Diversity_from_VCF_pyvcf_' + \
    'Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe.py'

# takes the provided SEX and combines it with chromosomes to generate group_chr
# this is only in order to keep group and chr associated in rule all
GROUP = [SEX, SEX, 'males']
GROUP_CHR = [x + '_' + y for x, y in zip(GROUP, CHR)]

# Rules -----------------------------------------------------------------------

# Global Rules ----------------------------------------------------------------
rule all:
    input:
        # chrX analyzed by region for all pops
        expand('04_window_analysis/results/' +
               '{pops}_{group_chr}_{filter_iter}_byRegion_{correction}' +
               '_diversity.bed',
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
        expand('06_figures/results/' +
               '{pops}_PAB_{filter_iter}_{window}_{correction}_diversity.png',
               pops=POPS, filter_iter=FILTER, window=WINDOW,
               correction=CORRECTION),
        # output for ld_window_analysis
        expand('06_figures/results/' +
               '{pops}_{group_chr}_{window}_windows_{ld_bin}_LDbins_' +
               '95bootstrapCI.png',
               pops=POPS, group_chr="chrX_females",
               window=WINDOW, ld_bin=LD_BIN)

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

rule subset_VCF:
    input:
        pop_file = '01_populations/results/{pop}_{group}',
        vcf_file = lambda wildcards: config['chromosomes'][wildcards.chr]
    output:
        temp(path.join('data', 'subset_{chr}_{pop}_{group}.vcf'))
    shell:
        "bcftools view -S {input.pop_file} {input.vcf_file} > {output}"

rule calculate_pi_chr8:
    input:
        group = '01_populations/results/{pop}_{group}',
        vcf = path.join('data', 'subset_chr8_{pop}_{group}.vcf')
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        path.join('02_diversity_by_site/results',
                  '{pop}_{group}_chr8_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc 8 "
        "--out_directory {params.out_dir}"

rule calculate_pi_chrX:
    input:
        group = '01_populations/results/{pop}_{group}',
        vcf = path.join('data', 'subset_chrX_{pop}_{group}.vcf')
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
            '03_filters/raw_filters/{filter_type}_{chrom}_filter.bed',
            filter_type=config["filters"][wildcards.filter_iter],
            chrom=[wildcards.chr])
    output:
        '03_filters/results/complete_{chr}_{filter_iter}.bed'
    shell:
        "cat {input} | sort -k1,1 -k2,2n | "
        "awk \'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3}}\' | "
        "bedtools merge -i stdin > {output}"

rule split_callable_sites:
    input:
        config['callable_sites']
    params:
        chrom = lambda wildcards: "\"" + wildcards.chr + "\""
    output:
        'data/callable_sites_{chr}.bed'
    shell:
        "cat {input} | awk \'$1 == {params.chrom} {{ print }}\' "
        "> {output}"

rule create_windows:
    input:
        'data/GRCh37_chromosome_coordinates.bed'
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
        "cat {input} | awk \'$1 == {params.chrom} {{ print }}\' | "
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
            '--sliding '
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                  '_diversity.unfiltered.bed')
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --windows {input.windows} "
        "{params.slide}--output {output}"

rule window_analysis_byRegion:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{pop}_{group}_{chr}_{filter_iter}' +
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_' +
                                      '{filter_iter}.bed')
    wildcard_constraints:
        window = "byRegion"
    params:
        window_calcs = '04_window_analysis/scripts/window_calculations.py'
    output:
        temp(path.join('04_window_analysis',
                       '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                       '_diversity.bed'))
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --chrX_windows --output {output}"

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

rule windowed_divergence_correction:
    input:
        diversity = path.join('04_window_analysis', 'results',
                              '{pop}_{group}_{chr}_{filter_iter}_{window}' +
                              '_diversity.bed'),
        divergence = path.join('data', 'substitution_rates',
                               '{chr}_{correction}_{filter_iter}_{window}' +
                               '_substitution_rates.txt')
    params:
        corrected = lambda wildcards: False if wildcards.correction == \
            "uncorrected" else True,
        script = '04_window_analysis/scripts/divergence_correction.py'
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}_' +
                  '{correction}_diversity.bed')
    run:
        if params.corrected is True:
            shell("python {params.script} --diversity {input.diversity} " +
                  "--divergence {input.divergence} --output {output}")
        else:
            shell("cp {input.diversity} {output}")

rule permute_chrX_regions:
    input:
        byWindow_100kb = path.join('04_window_analysis', 'results',
                                   '{pop}_{group}_{chr}_{filter_iter}' +
                                   '_100kb_{correction}_diversity.bed'),
        byRegion = path.join('04_window_analysis',
                             '{pop}_{group}_{chr}_{filter_iter}' +
                             '_byRegion_diversity.bed')
    params:
        permutation_script = path.join('04_window_analysis', 'scripts',
                                       'permute_chrX_windows.py'),
        replicates = 10000
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_byRegion_{correction}' +
                  '_diversity.bed')
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
        chrom = lambda wildcards: wildcards.chr
    output:
        path.join('06_figures', 'results',
                  '{pop}_{group}_{chr}_{filter_iter}_{window}_{correction}' +
                  '_diversity.png')
    shell:
        "Rscript {params.R_script} -i {input} -o {output} -c {params.chrom}"

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
                             'plot_chrX_maleFemale_diversity.R')
    output:
        path.join('06_figures', 'results',
                  '{pop}_chrX_malesAndFemales_{filter_iter}_{window}' +
                  '_{correction}_diversity.png')
    shell:
        "Rscript {params.R_script} --males {input.chrX_males} --females "
        "{input.chrX_females} -o {output}"

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
                             'plot_PAB_diversity.R')
    output:
        path.join('06_figures', 'results',
                  '{pop}_PAB_{filter_iter}_{window}_{correction}' +
                  '_diversity.png')
    shell:
        "Rscript {params.R_script} --chrX_females {input.chrX_females} "
        "--chrX_males {input.chrX_males} --chrY {input.chrY} "
        "--output {output}"

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
        path.join('data', 'subset_LD_{chr}_{pop}_{group}.vcf')
    output:
        temp(path.join('data', 'subset_{chr}_{pop}_{group}' +
                       '_snpsONLY_mac_filtered.vcf'))
    shadow: "shallow"
    shell:
        "bcftools view -m2 -M2 -v snps {input} | vcftools --vcf - "
        "--mac 1 --recode --out {output}"

rule calculate_ld:
    input:
        path.join('data', 'subset_{chr}_{pop}_{group}' +
                  '_snpsONLY_mac_filtered.vcf')
    output:
        temp(path.join('data', '{pop}_{chr}_{group}_filtered_ld_R2.ld'))
    shadow: "shallow"
    shell:
        "plink2 --vcf {input} --memory 4000 --r2 with-freqs "
        "--ld-window 700000 --ld-window-kb 800 --out {output}"

rule ld_window_analysis:
    input:
        LD = path.join('data', '{pop}_{chr}_{group}_filtered_ld_R2.ld'),
        windows = path.join('04_window_analysis', 'inputs',
                            '{chr}_{window}_window.bed'),
        script = path.join('05_ld_windows', 'scripts', 'ld_analysis.c')
    params:
        LD_bin = lambda wildcards: config['ld_bins'][wildcards.ld_bin]['size']
    output:
        path.join('05_ld_windows', 'results', '{pop}_{chr}_{group}_' +
                  '{window}_windows_{ld_bin}_LDbins_95bootstrapCI.txt')
    shell:
        "python average_ld_by_window.py --plink_ld {input.LD} "
        "--windows {input.windows} --binSize {params.LD_bin} "
        "--output {output}"

rule plot_ld_windows:
    input:
        path.join('05_ld_windows', 'results', '{pop}_{chr}_{group}_' +
                  '{window}_windows_{ld_bin}_LDbins_95bootstrapCI.txt')
    output:
        path.join('06_figures', 'results',
                  '{pop}_{chr}_{group}_{window}_windows_{ld_bin}' +
                  '_LDbins_95bootstrapCI.png')
    params:
        R_script = path.join('06_figures', 'scripts',
                             'plot_LD_bins.R'),
        winSize = lambda wildcards:
            config["windows"][wildcards.window]["win_size"],
    shell:
        "Rscript {params.R_script} -i {input} --winSize {params.winSize} "
        "--zoom 15 -o {output}"
