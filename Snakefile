"""
chrX_regional_variation
Snakefile
Daniel Cotter

analyze diversity and LD across chrX and chrY from 1000 genomes data
---------------------------------------------------------------------

Requires:
    conda
        to activate the included PAB_variation.yml environment
    snakemake
        included in the Conda environment
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
FILTER = ['filter4']
WINDOW = ['100kb']
LD_BIN = ['300kb']

# sets the populations to be a list of all pops and subpops
POPS = POPULATIONS + SUBPOPULATIONS

POPS = 'YRI'

# select a "sex" category to use for analysis of chrX and chr8
# use "males", "females", or "individuals" (for both)
SEX = 'individuals'

# select the pairwise substitution rates to use for divergence correction
CORRECTION = ['canFam3-hg19-corrected']


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

# Array containing the number of Mb for chrX ld 07_figures
PLOT_LENGTH = ['15', '156']

# Constrain wildcards globally

wildcard_constraints:
    pop = '[a-zA-Z]{3}'

# Rules -----------------------------------------------------------------------

# All -------------------------------------------------------------------------
rule all:
    input:
        fig1_a = expand(path.join('results', 'figures',
                                  '{pops}_{group}_chrX_{filter_iter}_{window}_'
                                  '{correction}_diversity.pdf'),
                        pops='YRI', group=SEX,
                        filter_iter=FILTER, window=WINDOW,
                        correction=CORRECTION),
        fig1_b = expand(path.join('results', 'figures',
                                  'subpops_X-A_PAR-A_XTR-A_ratios_'
                                  '{correction}_{filter_iter}.pdf'),
                        correction=CORRECTION,
                        filter_iter=FILTER),
        fig2 = path.join('results', 'figures',
                         'divergence_ratios_byFilter.pdf'),
        fig3 = expand(path.join('results', 'figures',
                                'LD_allSuperpops_byRegion_'
                                '{ld_bin}_LDbins_{filter_iter}.pdf'),
                      ld_bin=LD_BIN, filter_iter=FILTER),
        fig4a = expand(path.join('results', 'figures',
                                 'allPops_{correction}_diversityRatios_'
                                 'wDistanceFromGenes_unnormalized.pdf'),
                       correction=CORRECTION),
        fig4b = expand(path.join('results', 'figures',
                                 'allPops_{correction}_diversityRatios_'
                                 'wDistanceFromGenes_{denomPop}_'
                                 'normalized.pdf'),
                       correction=CORRECTION,
                       denomPop='MSL'),
        figS1a = expand(path.join('results', 'figures',
                                  '{pops}_chrX_{group}_{window}_windows_'
                                  '{correction}_{filter_iter}_{ld_bin}_'
                                  'LDbin_correlation.pdf'),
                        pops=POPS, group='females',
                        correction=CORRECTION, window=WINDOW,
                        ld_bin=LD_BIN, filter_iter='filter1'),
        figS1b = expand(path.join('results', 'figures',
                                  '{pops}_chrX_{group}_{window}_windows_'
                                  '{correction}_{filter_iter}_{ld_bin}_'
                                  'LDbin_correlation.pdf'),
                        pops=POPS, group='females',
                        correction=CORRECTION, window=WINDOW,
                        ld_bin=LD_BIN, filter_iter='filter4'),
        figS2 = expand(path.join('results', 'figures',
                                 'subpops_X-PAR_X-A_demography-corrected'
                                 '-ratios_{correction}_{filter_iter}.pdf'),
                       correction=CORRECTION, filter_iter=FILTER),
        figS3 = expand(path.join('results', 'figures',
                                 'LD_allSubpops_byRegion_'
                                 '{ld_bin}_LDbins_{filter_iter}.pdf'),
                       ld_bin=LD_BIN, filter_iter=FILTER),
        figS4A_D = expand(path.join('results', 'figures',
                                    'allPops_{correction}_diversityRatios_'
                                    'wDistanceFromGenes_{denomPop}_'
                                    'normalized.pdf'),
                          correction=CORRECTION,
                          denomPop=['TSI', 'PJL', 'KHV', 'PUR']),
        figures_S5 = path.join('results', 'figures', 'recombination_sim.pdf'),
        tableS2 = expand(path.join('results', 'tables',
                                   'allPops_{group}_{filter_iter}_'
                                   '{correction}_diversity_byRegion.csv'),
                         group=SEX, filter_iter=FILTER, correction=CORRECTION),
        tableS3 = expand(path.join('results', 'tables',
                                   'YRI_{correction}_diversity_'
                                   'withDistanceFromGenes.csv'),
                         correction=CORRECTION),
        tableS4 = path.join('results', 'tables',
                            'gene_density_by_region.csv')

# Global Rules ----------------------------------------------------------------
# create files for all populations stratified by sex
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

# subset the main VCF files for each of the populations we study
rule subset_VCF:
    input:
        pop_file = '01_populations/results/{Pop}_{group}',
        vcf_file = lambda wildcards: config['chromosomes'][wildcards.chr]
    output:
        temp(path.join('data', 'subset_{chr}_{Pop}_{group}.vcf'))
    shell:
        "bcftools view -S {input.pop_file} {input.vcf_file} > {output}"

# calculate diversity for the autosome (chr8 here)
rule calculate_pi_autosomes:
    input:
        group = '01_populations/results/{Pop}_{group}',
        vcf = path.join('data', 'subset_{chr}_{Pop}_{group}.vcf')
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/',
        chrom = lambda wildcards: wildcards.chr[3:]
    wildcard_constraints:
        chr = 'chr[0-9]+'
    output:
        path.join('02_diversity_by_site/results',
                  '{Pop}_{group}_{chr}_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc {params.chrom} "
        "--out_directory {params.out_dir}"

# calculate diversity on the X chromosome
rule calculate_pi_chrX:
    input:
        group = '01_populations/results/{Pop}_{group}',
        vcf = path.join('data', 'subset_chrX_{Pop}_{group}.vcf')
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        path.join('02_diversity_by_site/results',
                  '{Pop}_{group}_chrX_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc X "
        "--haploid --out_directory {params.out_dir}"

# calculate diversity on the Y chromosome
rule calculate_pi_chrY:
    input:
        group = '01_populations/results/{Pop}_{group}',
        vcf = path.join('data', 'subset_chrY_{Pop}_{group}.vcf')
    params:
        calc_pi = DIVERSITY_SCRIPT,
        out_dir = '02_diversity_by_site/results/'
    output:
        path.join('02_diversity_by_site/results',
                  '{Pop}_{group}_chrY_pi_output_by_site.txt')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.group} --chrom_inc Y "
        "--haploid --out_directory {params.out_dir}"

# intersect bed files to create an overlapping filer for downstream analysis
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

# get the callable sites file for a given chromosome
rule split_callable_sites:
    input:
        config['callable_sites']
    params:
        chrom = lambda wildcards: "\"" + wildcards.chr + "\""
    output:
        temp('data/callable_sites_{chr}.bed')
    shell:
        "grep {params.chrom} {input} > {output}"

# create a bed file with the requested windows tiling a chromosome
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

# convert diversity output to bed format
rule convert_diverstiy_to_bed:
    input:
        path.join('02_diversity_by_site', 'results',
                  '{Pop}_{group}_{chr}_pi_output_by_site.txt')
    params:
        bedConvert = path.join('02_diversity_by_site', 'scripts',
                               'bedConvert.py')
    output:
        temp(path.join('02_diversity_by_site', 'results',
                       '{Pop}_{group}_{chr}_pi_output_by_site.bed'))
    shell:
        "python {params.bedConvert} {input} {output}"

# remove any sites from the filter from the callable sites file
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

# grab only callable sites from the diversity file
rule filter_diversity_by_site:
    input:
        diversity_by_site = path.join('02_diversity_by_site', 'results',
                                      '{Pop}_{group}_{chr}'
                                      '_pi_output_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_{filter_iter}.bed')
    output:
        path.join('04_window_analysis', 'inputs',
                  '{Pop}_{group}_{chr}_{filter_iter}_pi_by_site.bed')
    shell:
        "bedtools intersect -a {input.diversity_by_site} "
        "-b {input.filtered_callable} > {output}"

# Window Analysis -------------------------------------------------------------
# calculate diversity in windows across a chromosome
rule window_analysis:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{Pop}_{group}_{chr}_{filter_iter}'
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_'
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
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}'
                  '_diversity.unfiltered.bed')
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --windows {input.windows} "
        "{params.slide}--replicates {params.replicates} --output {output}"

# calculate diversity in regions across the X chromosome
# see results/tables/region_coordinates.txt for coordinates
rule window_analysis_byRegion:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{Pop}_{group}_{chr}_{filter_iter}'
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_'
                                      '{filter_iter}.bed')
    wildcard_constraints:
        window = 'byRegion'
    params:
        window_calcs = '04_window_analysis/scripts/window_calculations.py',
        replicates = 1000
    output:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}'
                  '_diversity.bed')
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --chrX_windows "
        "--replicates {params.replicates} --output {output}"

# filter the output by removing windows with too few callable sites
rule filter_windows_by_callable_sites:
    input:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}'
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
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}'
                  '_diversity.bed')
    shell:
        "python {params.script} --input {input} --windowSize {params.winSize} "
        "--filter 0.1 --output {output}"

# create a fake divergence file for uncorrected output
rule create_pseudo_uncorrected_divergence:
    output:
        temp(path.join('data', 'substitution_rates',
                       '{chr}_uncorrected_{filter_iter}_{window}'
                       '_substitution_rates.txt'))
    shell:
        "touch {output}"

# correct substitution rate files from data/substitution_rates folder
rule divergence_JukesCantor69_correction:
    input:
        path.join('data', 'substitution_rates',
                  '{chr}_{correction}_{filter_iter}_{window}'
                  '_substitution_rates.txt')
    params:
        script = path.join('data', 'substitution_rates', 'scripts',
                           'Calculate_JC69_from_Galaxy_output.py'),
        corrected = lambda wildcards: False if wildcards.correction == \
            "uncorrected" else True
    output:
        temp(path.join('data', 'substitution_rates',
                       '{chr}_{correction}_{filter_iter}_{window}'
                       '_substitution_rates_JC69.txt'))
    run:
        if params.corrected is True:
            shell("python {params.script} --input_file {input} "
                  "--output_file {output}")
        else:
            shell("mv {input} {output}")

# divide diversity by divergence by window or return uncorrected output
rule windowed_divergence_correction:
    input:
        diversity = path.join('04_window_analysis', 'results',
                              '{Pop}_{group}_{chr}_{filter_iter}_{window}'
                              '_diversity.bed'),
        divergence = path.join('data', 'substitution_rates',
                               '{chr}_{correction}_{filter_iter}_{window}'
                               '_substitution_rates_JC69.txt')
    params:
        corrected = lambda wildcards: False if wildcards.correction == \
            "uncorrected" else True,
        script = '04_window_analysis/scripts/divergence_correction.py'
    output:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}_'
                  '{correction}_diversity.bed')
    run:
        if params.corrected:
            shell("python {params.script} --diversity {input.diversity} "
                  "--divergence {input.divergence} --output {output}")
        else:
            shell("cp {input.diversity} {output}")

# permute chrX regions to calculate CIs
rule permute_chrX_regions:
    input:
        byWindow_100kb = path.join('04_window_analysis', 'results',
                                   '{Pop}_{group}_{chr}_{filter_iter}'
                                   '_100kb_{correction}_diversity.bed'),
        byRegion = path.join('04_window_analysis', 'results',
                             '{Pop}_{group}_{chr}_{filter_iter}'
                             '_byRegion_{correction}_diversity.bed')
    params:
        permutation_script = path.join('04_window_analysis', 'scripts',
                                       'permute_chrX_windows.py'),
        replicates = 10000
    output:
        path.join('04_window_analysis', 'results', 'subpops_wPvals',
                  '{Pop}_{group}_{chr}_{filter_iter}_byRegion_{correction}'
                  '_diversity_wPvals.bed')
    shell:
        "python {params.permutation_script} --byRegion {input.byRegion} "
        "--byWindow {input.byWindow_100kb} --replicates {params.replicates} "
        "--output {output}"

# plot diversity by windows
rule plot_windowed_diversity:
    input:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}_{correction}'
                  '_diversity.bed')
    params:
        R_script = path.join('07_figures', 'scripts',
                             'plot_windowed_diversity.R'),
        chrom = lambda wildcards: wildcards.chr,
        height = lambda wildcards: config[wildcards.correction][wildcards.chr]
    output:
        path.join('07_figures', 'results',
                  '{Pop}_{group}_{chr}_{filter_iter}_{window}_{correction}'
                  '_diversity.{ext}')
    shell:
        "Rscript {params.R_script} -i {input} -o {output} -c {params.chrom} "
        "--maxHeight {params.height}"

# plot diversity on the X chromosome stratified by sex -- UNUSED
rule plot_sex_specific_chrX_windows:
    input:
        chrX_males = path.join('04_window_analysis', 'results',
                               '{Pop}_males_chrX_{filter_iter}_{window}'
                               '_{correction}_diversity.bed'),
        chrX_females = path.join('04_window_analysis', 'results',
                                 '{Pop}_females_chrX_{filter_iter}_{window}'
                                 '_{correction}_diversity.bed')
    params:
        R_script = path.join('07_figures', 'scripts',
                             'plot_chrX_maleFemale_diversity.R'),
        height = lambda wildcards: config[wildcards.correction]['chrX']
    output:
        path.join('07_figures', 'results',
                  '{Pop}_chrX_malesAndFemales_{filter_iter}_{window}'
                  '_{correction}_diversity.{ext}')
    shell:
        "Rscript {params.R_script} --males {input.chrX_males} --females "
        "{input.chrX_females} --maxHeight {params.height} -o {output}"

# plot diversity across the pseudoautosomal bounday -- UNUSED
rule plot_PAB_diversity:
    input:
        chrX_males = path.join('04_window_analysis', 'results',
                               '{Pop}_males_chrX_{filter_iter}_{window}'
                               '_{correction}_diversity.bed'),
        chrX_females = path.join('04_window_analysis', 'results',
                                 '{Pop}_females_chrX_{filter_iter}_{window}'
                                 '_{correction}_diversity.bed'),
        chrY = path.join('04_window_analysis', 'results',
                         '{Pop}_males_chrY_{filter_iter}_{window}'
                         '_{correction}_diversity.bed')
    params:
        R_script = path.join('07_figures', 'scripts',
                             'plot_PAB_diversity.R'),
        height = lambda wildcards: config[wildcards.correction]['chrX']
    output:
        path.join('07_figures', 'results',
                  '{Pop}_PAB_{filter_iter}_{window}_{correction}'
                  '_diversity.{ext}')
    shell:
        "Rscript {params.R_script} --chrX_females {input.chrX_females} "
        "--chrX_males {input.chrX_males} --chrY {input.chrY} "
        "--maxHeight {params.height} --output {output}"

# Analysis by Region ----------------------------------------------------------
# get a bed file defining the whole chromosome
rule get_wholeChr_bed:
    input:
        path.join('data', 'hg19_chromosome_coordinates.bed')
    params:
        chrom = lambda wildcards: wildcards.chr
    output:
        temp(path.join('data', '{chr}_wholeChr.bed'))
    shell:
        'grep {params.chrom} {input} > {output}'

# use the whole chromosome as the window for the window analysis
rule window_analysis_wholeChr:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{Pop}_{group}_{chr}_{filter_iter}'
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_'
                                      '{filter_iter}.bed'),
        windows = path.join('data', '{chr}_{window}.bed')
    wildcard_constraints:
        window = 'wholeChr'
    params:
        window_calcs = '04_window_analysis/scripts/window_calculations.py',
        replicates = 1000
    output:
        temp(path.join('04_window_analysis', 'results',
                       '{Pop}_{group}_{chr}_{filter_iter}_{window}'
                       '_diversity.bed'))
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --windows {input.windows} "
        "--replicates {params.replicates} --output {output}"

# calculate ratios for downstream analysis
rule calculate_A_ratios:
    input:
        chrX = path.join('04_window_analysis', 'results',
                         '{Pop}_{group}_chrX_{filter_iter}_byRegion'
                         '_{correction}_diversity.bed'),
        chrY = path.join('04_window_analysis', 'results',
                         '{Pop}_males_chrY_{filter_iter}_wholeChr'
                         '_{correction}_diversity.bed'),
        chr8 = path.join('04_window_analysis', 'results',
                         '{Pop}_{group}_chr8_{filter_iter}_wholeChr'
                         '_{correction}_diversity.bed')
    params:
        script = '07_figures/scripts/ratios_table.py'
    output:
        path.join('07_figures', 'results',
                  '{Pop}_{group}_{filter_iter}_{correction}_ratios.txt')
    shell:
        "python {params.script} --chrX {input.chrX} --chrY {input.chrY} "
        "--chr8 {input.chr8} --output {output}"

# prepare a summary file for all ratio data for plotting
rule prepare_subpop_ratio_plotting_data:
    input:
        chrX = lambda wildcards: \
            expand(path.join('04_window_analysis', 'results',
                             '{pops}_{Group}_chrX_{filter}_byRegion'
                             '_{correction_div}_diversity.bed'),
                   pops=SUBPOPULATIONS, Group=SEX,
                   filter=wildcards.filter_iter,
                   correction_div=wildcards.correction),
        chrY = lambda wildcards: \
            expand(path.join('04_window_analysis', 'results',
                             '{pops}_males_chrY_{filter}_wholeChr'
                             '_{correction_div}_diversity.bed'),
                   pops=SUBPOPULATIONS, filter=wildcards.filter_iter,
                   correction_div=wildcards.correction),
        chr8 = lambda wildcards: \
            expand(path.join('04_window_analysis', 'results',
                             '{pops}_{Group}_chr8_{filter}_wholeChr'
                             '_{correction_div}_diversity.bed'),
                   pops=SUBPOPULATIONS, Group=SEX,
                   filter=wildcards.filter_iter,
                   correction_div=wildcards.correction)
    output:
        path.join('07_figures', 'results',
                  'subpops_{filter_iter}_{correction}_ratios_table.txt')
    params:
        script = path.join('07_figures', 'scripts',
                           'format_byRegion_data.py')
    shell:
        "python {params.script} --chrX_byRegion {input.chrX} "
        "--chr8_wholeChr {input.chr8} --chrY_wholeChr {input.chrY} "
        "--output {output}"

# plot the ratios from the summary file above
rule plot_A_Ratios_across_subpops:
    input:
        path.join('07_figures', 'results',
                  'subpops_{filter_iter}_{correction}_ratios_table.txt')
    output:
        path.join('07_figures', 'results',
                  'subpops_X-A_PAR-A_XTR-A_ratios_{correction}'
                  '_{filter_iter}.{ext}')
    params:
        Rscript = path.join('07_figures', 'scripts',
                            'plot_subpop_diversity_ratios.R')
    shell:
        "Rscript {params.Rscript} --subpops_data {input} -o {output}"

# plor the ratios from the summary file above
rule plot_relative_XvPAR_XvA_ratios:
    input:
        path.join('07_figures', 'results',
                  'subpops_{filter_iter}_{correction}_ratios_table.txt')
    output:
        path.join('07_figures', 'results',
                  'subpops_X-PAR_X-A_demography-corrected-ratios'
                  '_{correction}_{filter_iter}.{ext}')
    params:
        Rscript = path.join('07_figures', 'scripts',
                            'relative_XA_ratios_bySubpop.R')
    shell:
        "Rscript {params.Rscript} --subpops_data {input} -o {output}"

# move diversity data into a single folder
rule move_chr8_diversity_data:
    input:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_chr8_{filter_iter}_wholeChr'
                  '_{correction}_diversity.bed')
    output:
        path.join('04_window_analysis', 'results', 'chr8_data',
                  '{Pop}_{group}_chr8_{filter_iter}_wholeChr'
                  '_{correction}_diversity.bed')
    shell:
        "cp {input} {output}"

# create a supplemental table for diversity across all populations
# grab chrX data from output of permutation script
rule create_supp_table_allPops:
    input:
        chrX = lambda wildcards: expand(
            path.join(
                '04_window_analysis', 'results',
                'subpops_wPvals',
                '{Pop}_{group}_chrX_{filter_iter}'
                '_byRegion_{correction}_diversity_wPvals.bed'),
            pop=SUBPOPULATIONS, group=SEX,
            filter_iter=wildcards.filter_iter,
            correction=wildcards.correction),
        chr8 = lambda wildcards: expand(
            path.join('04_window_analysis', 'results', 'chr8_data',
                      '{pops}_{Group}_chr8_{filter}_wholeChr'
                      '_{correction_div}_diversity.bed'),
            pops=SUBPOPULATIONS, Group=SEX,
            filter=wildcards.filter_iter,
            correction_div=wildcards.correction)
    params:
        Rscript = path.join('04_window_analysis', 'scripts',
                            'make_diversity_tables_allPops.R'),
        folder1 = path.join('04_window_analysis', 'results', 'subpops_wPvals'),
        folder2 = path.join('04_window_analysis', 'results', 'chr8_data')
    output:
        path.join('04_window_analysis', 'results',
                  'allPops_{group}_{filter_iter}_{correction}_diversity_'
                  'byRegion.csv')
    shell:
        "Rscript {params.Rscript} --folder1 {params.folder1} "
        "--folder2 {params.folder2} -o {output}"

# move files for creating the distance from genes supplmental table
rule move_supp_table_files:
    input:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_{chr}_{filter_iter}_wholeChr'
                  '_{correction}_diversity.bed')
    output:
        path.join('04_window_analysis', 'results', 'diversity_byFilter',
                  '{Pop}_{group}_{chr}_{filter_iter}_wholeChr'
                  '_{correction}_diversity.bed')
    shell:
        "cp {input} {output}"

# move files for creating the distance from genes supplmental table
rule move_chrX_supp_table_files:
    input:
        path.join('04_window_analysis', 'results',
                  '{Pop}_{group}_chrX_{filter_iter}_byRegion'
                  '_{correction}_diversity.bed')
    output:
        path.join('04_window_analysis', 'results', 'diversity_byFilter',
                  '{Pop}_{group}_chrX_{filter_iter}_byRegion'
                  '_{correction}_diversity.bed')
    shell:
        "cp {input} {output}"

# use files to create distance from genes supplemental table
rule create_supp_table_distanceFromGenes:
    input:
        chrX = lambda wildcards: expand(
            path.join(
                '04_window_analysis', 'results',
                'diversity_byFilter',
                '{Pop}_{group}_chrX_{filter_iter}'
                '_byRegion_{correction}_diversity.bed'),
            pop='YRI', group=SEX,
            filter_iter=['filter' + str(i) for i in range(1, 8)],
            correction=['uncorrected', wildcards.correction]),
        chrY = lambda wildcards: expand(
            path.join(
                '04_window_analysis', 'results',
                'diversity_byFilter',
                '{Pop}_{group}_chrY_{filter_iter}'
                '_wholeChr_{correction}_diversity.bed'),
            pop='YRI', group='males',
            filter_iter=['filter' + str(i) for i in range(1, 8)],
            correction=['uncorrected', wildcards.correction]),
        chr8 = lambda wildcards: expand(
            path.join(
                '04_window_analysis', 'results',
                'diversity_byFilter',
                '{Pop}_{group}_chrY_{filter_iter}'
                '_wholeChr_{correction}_diversity.bed'),
            pop='YRI', group='males',
            filter_iter=['filter' + str(i) for i in range(1, 8)],
            correction=['uncorrected', wildcards.correction])
    params:
        Rscript = path.join('04_window_analysis', 'scripts',
                            'make_diversity_table_distanceFromGenes.R'),
        folder = path.join('04_window_analysis', 'results',
                           'diversity_byFilter')
    output:
        path.join('04_window_analysis', 'results',
                  'YRI_{correction}_diversity_withDistanceFromGenes.csv')
    shell:
        "Rscript {params.Rscript} --folder {params.folder} -o {output}"

# LD analysis -----------------------------------------------------------------
# compile the LD window calculation modules using cython
rule cythonize_ld_script:
    input:
        ld_analysis = '05_ld_windows/scripts/ld_analysis.pyx',
        setup = '05_ld_windows/scripts/setup.py'
    output:
        '05_ld_windows/scripts/ld_analysis.c'
    shell:
        "python {input.setup} build_ext --inplace"

# subset the VCF for each population of interest
rule subset_VCF_for_LD:
    input:
        pop_file = '01_populations/results/{Pop}_{group}',
        vcf_file = lambda wildcards: config['chromosomes'][wildcards.chr]
    output:
        temp(path.join('data', 'subset_LD_{chr}_{Pop}_{group}.vcf'))
    shell:
        "bcftools view -S {input.pop_file} {input.vcf_file} > {output}"

# filter the VCF for snpsonly and minor allele count > 1
rule filter_vcf_for_LD:
    input:
        vcf = path.join('data', 'subset_LD_{chr}_{Pop}_{group}.vcf'),
        targets = path.join('03_filters', 'results',
                            'complete_{chr}_{filter_iter}.bed')
    output:
        path.join('data', 'subset_LD_{chr}_{Pop}_{group}_{filter_iter}'
                  '_snpsONLY-mac-filtered.recode.vcf')
    params:
        prefix = path.join('data', 'subset_LD_{chr}_{Pop}_{group}_'
                           '{filter_iter}_snpsONLY-mac-filtered'),
        filter = lambda wildcards: wildcards.filter_iter
    shadow: "full"
    shell:
        "cat {input.targets} | sed 's/^chr//' > temp_{params.filter}.bed && "
        "bcftools view -Ou -m2 -M2 -v snps {input.vcf} | bcftools view -Ov "
        "--min-ac 1:minor | vcftools --vcf - "
        "--exclude-bed temp_{params.filter}.bed --recode --out {params.prefix}"

# calculate LD using plink at all variants within 1 MB of each other
rule calculate_ld:
    input:
        path.join('data', 'subset_LD_{chr}_{Pop}_{group}_{filter_iter}'
                  '_snpsONLY-mac-filtered.recode.vcf')
    params:
        out_path = path.join('data', '{Pop}_{chr}_{group}_{filter_iter}'
                             '_filtered_ld_R2')
    output:
        temp(path.join('data', '{Pop}_{chr}_{group}_{filter_iter}'
                       '_filtered_ld_R2.ld'))
    shadow: "full"
    threads: 4
    shell:
        "plink --vcf {input} --r2 with-freqs --memory 8000 "
        "--threads {threads} "
        "--ld-window 700000 --ld-window-kb 1100 --out {params.out_path}"

# calculate average LD by window
rule ld_window_analysis:
    input:
        LD = path.join('data', '{Pop}_{chr}_{group}_{filter_iter}'
                       '_filtered_ld_R2.ld'),
        windows = path.join('04_window_analysis', 'inputs',
                            '{chr}_{window}_window.bed'),
        script = path.join('05_ld_windows', 'scripts', 'ld_analysis.c')
    wildcard_constraints:
        window = '[0-9]+[A-Za-z]+'
    params:
        LD_bin = lambda wildcards: config['ld_bins'][wildcards.ld_bin]['size'],
        script = '05_ld_windows/scripts/average_ld_by_window.py'
    output:
        path.join('05_ld_windows', 'results', '{Pop}_{chr}_{group}_'
                  '{window}_windows_{filter_iter}'
                  '_{ld_bin}_LDbins_95bootstrapCI.txt')
    shell:
        "python {params.script} --plink_ld {input.LD} "
        "--windows {input.windows} --binSize {params.LD_bin} "
        "--output {output}"

# calculate average LD by window for each chrX region
rule ld_window_analysis_byRegion:
    input:
        LD = path.join('data', '{Pop}_{chr}_{group}_{filter_iter}'
                       '_filtered_ld_R2.ld'),
        script = path.join('05_ld_windows', 'scripts', 'ld_analysis.c')
    wildcard_constraints:
        window = 'byRegion'
    params:
        LD_bin = lambda wildcards: config['ld_bins'][wildcards.ld_bin]['size'],
        script = '05_ld_windows/scripts/average_ld_by_window.py'
    output:
        path.join('05_ld_windows', 'results', '{Pop}_{chr}_{group}_'
                  '{window}_windows_{filter_iter}'
                  '_{ld_bin}_LDbins_95bootstrapCI.txt')
    shell:
        "python {params.script} --plink_ld {input.LD} "
        "--byRegion --binSize {params.LD_bin} "
        "--output {output}"

# calculate average LD across whole chromosomes
rule ld_window_analysis_wholeChr_Y_autosomes:
    input:
        LD = path.join('data', '{Pop}_{chr}_{group}_{filter_iter}'
                       '_filtered_ld_R2.ld'),
        windows = path.join('data', '{chr}_{window}.bed'),
        script = path.join('05_ld_windows', 'scripts', 'ld_analysis.c')
    wildcard_constraints:
        window = 'wholeChr'
    params:
        LD_bin = lambda wildcards: config['ld_bins'][wildcards.ld_bin]['size'],
        script = '05_ld_windows/scripts/average_ld_by_window.py'
    output:
        path.join('05_ld_windows', 'results', '{Pop}_{chr}_{group}_'
                  '{window}_windows_{filter_iter}'
                  '_{ld_bin}_LDbins_95bootstrapCI.txt')
    shell:
        "python {params.script} --plink_ld {input.LD} --windows "
        "{input.windows} --binSize {params.LD_bin} --output {output}"

# plot LD by window -- UNUSED
rule plot_ld_windows:
    input:
        path.join('05_ld_windows', 'results', '{Pop}_{chr}_{group}_'
                  '{window}_windows_{filter_iter}'
                  '_{ld_bin}_LDbins_95bootstrapCI.txt')
    output:
        path.join('07_figures', 'results',
                  '{Pop}_{chr}_{group}_{window}_windows_{filter_iter}'
                  '_{ld_bin}_LDbins_95bootstrapCI_{LDplotSize}Mb.{ext}')
    params:
        R_script = path.join('07_figures', 'scripts',
                             'plot_LD_bins.R'),
        winSize = lambda wildcards:
            config["windows"][wildcards.window]["win_size"],
        plot_size = lambda wildcards: wildcards.LDplotSize
    shell:
        "Rscript {params.R_script} -i {input} --winSize {params.winSize} "
        "--zoom {params.plot_size} -o {output} --maxHeight 1"

# plot the correlation between LD and diversity
rule plot_ld_pi_correlation:
    input:
        ld = path.join('05_ld_windows', 'results', '{Pop}_{chr}_{group}_'
                       '{window}_windows_{filter_iter}'
                       '_{ld_bin}_LDbins_95bootstrapCI.txt'),
        pi = path.join('04_window_analysis', 'results',
                       '{Pop}_{group}_{chr}_{filter_iter}_{window}_'
                       '{correction}_diversity.bed')
    params:
        R_script = path.join('07_figures', 'scripts',
                             'ld_pi_correlation.R'),
        distance_filtered = lambda wildcards: \
            config["filter_descriptions"][wildcards.filter_iter]
    output:
        path.join('07_figures', 'results',
                  '{Pop}_{chr}_{group}_{window}_windows_{correction}'
                  '_{filter_iter}_{ld_bin}_LDbin_correlation.{ext}')
    shell:
        "Rscript {params.R_script} --LD {input.ld} --diversity {input.pi} "
        "--filter {params.distance_filtered} --output {output}"

# plot the correlation between LD and diversity without filtering -- UNUSED
rule plot_ld_pi_correlation_noFilter:
    input:
        ld = path.join('05_ld_windows', 'results', '{Pop}_{chr}_{group}_'
                       '{window}_windows_filter0'
                       '_{ld_bin}_LDbins_95bootstrapCI.txt'),
        pi = path.join('04_window_analysis', 'results',
                       '{Pop}_{group}_{chr}_{filter_iter}_{window}_'
                       '{correction}_diversity.bed')
    params:
        R_script = path.join('07_figures', 'scripts',
                             'ld_pi_correlation.R'),
        distance_filtered = lambda wildcards: \
            config["filter_descriptions"][wildcards.filter_iter]
    output:
        path.join('07_figures', 'results',
                  '{Pop}_{chr}_{group}_{window}_windows_{correction}'
                  '_{filter_iter}_{ld_bin}_LDbin_correlation_'
                  'noFilterLD.{ext}')
    shell:
        "Rscript {params.R_script} --LD {input.ld} --diversity {input.pi} "
        "--filter {params.distance_filtered} --output {output}"

# Plot Results ----------------------------------------------------------------
# prepare LD data across all pops for chrX and chr8 to be plotted
rule prepare_LD_data_allPops:
    input:
        chrX = lambda wildcards: \
            expand(path.join('05_ld_windows', 'results',
                             '{pops}_chrX_females_byRegion_windows_{filter}'
                             '_{LDbin}_LDbins_95bootstrapCI.txt'),
                   pops=SUBPOPULATIONS, filter=wildcards.filter_iter,
                   LDbin=wildcards.ld_bin),
        chr8 = lambda wildcards: \
            expand(path.join('05_ld_windows', 'results',
                             '{pops}_chr8_individuals_wholeChr_windows_'
                             '{filter}_{LDbin}_LDbins_95bootstrapCI.txt'),
                   pops=SUBPOPULATIONS, filter=wildcards.filter_iter,
                   LDbin=wildcards.ld_bin)
    params:
        script = path.join('07_figures', 'scripts',
                           'format_byRegion_LD_data.py')
    output:
        path.join('07_figures', 'results',
                  'LD_all_subpops_{filter_iter}_{ld_bin}_LDbins_'
                  '95bootstrapCI.txt')
    shell:
        "python {params.script} --chrX_byRegion {input.chrX} --chr8_wholeChr "
        "{input.chr8} --output {output}"

# prepare LD data across all superpops for chrX and chr8 to be plotted
rule prepare_LD_data_allSuperpops:
    input:
        chrX = lambda wildcards: \
            expand(path.join('05_ld_windows', 'results',
                             '{pops}_chrX_females_byRegion_windows_{filter}'
                             '_{LDbin}_LDbins_95bootstrapCI.txt'),
                   pops=POPULATIONS, filter=wildcards.filter_iter,
                   LDbin=wildcards.ld_bin),
        chr8 = lambda wildcards: \
            expand(path.join('05_ld_windows', 'results',
                             '{pops}_chr8_individuals_wholeChr_windows_'
                             '{filter}_{LDbin}_LDbins_95bootstrapCI.txt'),
                   pops=POPULATIONS, filter=wildcards.filter_iter,
                   LDbin=wildcards.ld_bin)
    params:
        script = path.join('07_figures', 'scripts',
                           'format_byRegion_LD_data.py')
    output:
        path.join('07_figures', 'results',
                  'LD_all_superpops_{filter_iter}_{ld_bin}_LDbins_'
                  '95bootstrapCI.txt')
    shell:
        "python {params.script} --chrX_byRegion {input.chrX} --chr8_wholeChr "
        "{input.chr8} --output {output}"

# plot ratios of divergence to between regions with different filers
# This file is provided in the data directory and is composed of ratios of
# similarly provided files
rule plot_divergence_ratios_byFilter:
    input:
        path.join('data', 'substitution_rates',
                  'divergence_ratios_JC69-corrected_filter1-5.txt')
    params:
        script = path.join('07_figures', 'scripts',
                           'divergence_diversity_ratios_byFilter.R')
    output:
        path.join('07_figures', 'results',
                  'divergence_ratios_byFilter.{ext}')
    shell:
        "Rscript {params.script} --divergenceRatios {input} -o {output}"

# plot LD by region for a given population
rule plot_LD_byRegion:
    input:
        chrX = path.join('05_ld_windows', 'results',
                         '{Pop}_chrX_females_byRegion_windows_{filter_iter}_'
                         '{ld_bin}_LDbins_95bootstrapCI.txt'),
        chr8 = path.join('05_ld_windows', 'results',
                         '{Pop}_chr8_individuals_wholeChr_windows_'
                         '{filter_iter}_{ld_bin}_LDbins_95bootstrapCI.txt')
    params:
        script = path.join('07_figures', 'scripts', 'plot_LD_byRegion.R')
    output:
        path.join('07_figures', 'results', '{Pop}_LD_byRegion_{ld_bin}_'
                  'LDbins_{filter_iter}.{ext}')
    shell:
        "Rscript {params.script} --chrX {input.chrX} "
        "--chr8 {input.chr8} -o {output}"

# plot LD by region for all populations
rule plot_LD_byRegion_allPops:
    input:
        path.join('07_figures', 'results',
                  'LD_all_subpops_{filter_iter}_{ld_bin}_LDbins_'
                  '95bootstrapCI.txt')
    params:
        script = path.join('07_figures', 'scripts', 'plot_LD_byPop.R')
    output:
        path.join('07_figures', 'results', 'LD_allSubpops_byRegion_{ld_bin}_'
                  'LDbins_{filter_iter}.{ext}')
    shell:
        "Rscript {params.script} --subpops_data {input} -o {output}"

# plot LD by region grouped for all superpopulations
rule plot_LD_byRegion_allSuperpops:
    input:
        path.join('07_figures', 'results',
                  'LD_all_superpops_{filter_iter}_{ld_bin}_LDbins_'
                  '95bootstrapCI.txt')
    params:
        script = path.join('07_figures', 'scripts', 'plot_LD_bySuperpop.R')
    output:
        path.join('07_figures', 'results', 'LD_allSuperpops_byRegion_'
                  '{ld_bin}_LDbins_{filter_iter}.{ext}')
    shell:
        "Rscript {params.script} --pops_data {input} -o {output}"

# plot diversity ratios as a function of distance from genes
rule plot_distance_fromGenes_unnormalized:
    input:
        filter1 = path.join('07_figures', 'results',
                            'subpops_filter1_{correction}_ratios_table.txt'),
        filter2 = path.join('07_figures', 'results',
                            'subpops_filter2_{correction}_ratios_table.txt'),
        filter3 = path.join('07_figures', 'results',
                            'subpops_filter3_{correction}_ratios_table.txt'),
        filter4 = path.join('07_figures', 'results',
                            'subpops_filter4_{correction}_ratios_table.txt'),
        filter5 = path.join('07_figures', 'results',
                            'subpops_filter5_{correction}_ratios_table.txt')
    params:
        script = path.join('07_figures', 'scripts',
                           'distanceFromGenes_allPops.R')
    output:
        o1 = path.join('07_figures', 'results',
                       'allPops_{correction}_diversityRatios'
                       '_wDistanceFromGenes_unnormalized.{ext}')
    shell:
        "Rscript {params.script} --filter1 {input.filter1} --filter2 "
        "{input.filter2} --filter3 {input.filter3} --filter4 {input.filter4} "
        "--filter5 {input.filter5} --output1 {output.o1}"

# plot diversity ratios as a function of distance from genes
# after specifying a population to normalize to (default = "MSL")
rule plot_distance_fromGenes_normalized:
    input:
        filter1 = path.join('07_figures', 'results',
                            'subpops_filter1_{correction}_ratios_table.txt'),
        filter2 = path.join('07_figures', 'results',
                            'subpops_filter2_{correction}_ratios_table.txt'),
        filter3 = path.join('07_figures', 'results',
                            'subpops_filter3_{correction}_ratios_table.txt'),
        filter4 = path.join('07_figures', 'results',
                            'subpops_filter4_{correction}_ratios_table.txt'),
        filter5 = path.join('07_figures', 'results',
                            'subpops_filter5_{correction}_ratios_table.txt')
    params:
        script = path.join('07_figures', 'scripts',
                           'distanceFromGenes_allPops.R'),
        denom_pop = lambda wildcards: wildcards.denomPop
    output:
        o2 = path.join('07_figures', 'results',
                       'allPops_{correction}_diversityRatios_'
                       'wDistanceFromGenes_{denomPop}_normalized.{ext}')
    shell:
        "Rscript {params.script} --filter1 {input.filter1} --filter2 "
        "{input.filter2} --filter3 {input.filter3} --filter4 {input.filter4} "
        "--filter5 {input.filter5} --output2 {output.o2} "
        "--denom_pop {params.denom_pop}"

# FST Analyses ----------------------------------------------------------------
# calculate pairwise FST in regions and output a table -- UNUSED
rule calculate_FST:
    input:
        chrX = path.join('data', 'subset_LD_chrX_ALL_females_{filter_iter}'
                         '_snpsONLY-mac-filtered.recode.vcf'),
        chr8 = path.join('data', 'subset_LD_chr8_ALL_females_{filter_iter}'
                         '_snpsONLY-mac-filtered.recode.vcf'),
        pop_data = config['panel']
    params:
        script = path.join('06_FST', 'scripts', 'calculateFST_byRegion.R'),
        comparisons = lambda wildcards: wildcards.comparisons
    output:
        tempDir = temp(
            directory(path.join('06_FST',
                                '{comparisons}_{filter_iter}_temp'))),
        FST_table = path.join('06_FST', 'results',
                              '{comparisons}_{filter_iter}_FST_byRegion.csv')
    shell:
        "mkdir {output.tempDir} && "
        "Rscript {params.script} --chrX_data {input.chrX} "
        "--chr8_data {input.chr8} --pop_data {input.pop_data} "
        "--temp_directory {output.tempDir} --comparisons {params.comparisons} "
        "--output {output.FST_table}"

# misc ------------------------------------------------------------------------
# calculate the fraction of sequence tiled by gene annotations
rule calculate_gene_density:
    input:
        path.join('03_filters', 'raw_filters',
                  'UCSC_refseq_gencode_gene_transcripts_hg19.bed')
    params:
        script = path.join('04_window_analysis', 'scripts', 'gene_density.R')
    output:
        path.join('04_window_analysis', 'results',
                  'gene_density_by_region.csv')
    shell:
        "Rscript {params.script} --input {input} --output {output}"

# plot simulation results
rule plot_sim_results:
    input:
        path.join('08_simulations', 'results', 'SLiM_output.txt')
    output:
        path.join('07_figures', 'results', 'recombination_sim.pdf')
    params:
        script = path.join('08_simulations', 'scripts', 'plot_simulations.R')
    shell:
        "Rscript {params.script} --input {input} --output {output}"

# move results to final folder ------------------------------------------------
# move all figures to results/figures folder (uploaded to github)
rule move_figure_output:
    input:
        path.join('07_figures', 'results', '{file}')
    output:
        path.join('results', 'figures', '{file}')
    shell:
        "cp {input} {output}"

# move all tables to results/tables folder (uploaded to github)
rule move_table_output:
    input:
        path.join('04_window_analysis', 'results', '{file}')
    output:
        path.join('results', 'tables', '{file}')
    shell:
        "cp {input} {output}"
