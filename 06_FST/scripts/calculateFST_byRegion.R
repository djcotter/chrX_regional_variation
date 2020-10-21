# install any required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("SNPRelate", quietly = TRUE)) {
  BiocManager::install("SNPRelate", update = FALSE)
}
if (!requireNamespace("SeqArray", quietly = TRUE)) {
  BiocManager::install("SeqArray", update = FALSE)
}

# source required packages
library(SNPRelate)
library(SeqArray)
library(optparse)
library(tidyverse)

# parse arguments from the command line
option_list = list(
  make_option(c('--chrX_data'), type='character', default=NULL,
              help="path to chrX VCF"),
  make_option(c('--chr8_data'), type='character', default=NULL,
              help="path to chr8 VCF"),
  make_option(c('--pop_data'), type='character', default=NULL,
              help="path to panel file"),
  make_option(c('--temp_directory'), type='chracter', default=NULL,
              help="path to temporary output directory"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output files"),
  make_option(c('--comparisons'), type='character', default='SUPERPOPS',
              help='Compare "POPS" or "SUPERPOPS"'),
  make_option(c('--num_samples'), type='integer', default=NA,
              help='number to downsample populations. set to 0 to use all pops.')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$chrX_data) || is.null(opt$chr8_data) || is.null(opt$pop_data) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and Output files must be specified.", call.=FALSE)
}

if(is.null(opt$temp_directory)) {
  print_help(opt_parser)
  stop("Temp directory must be specified.", call.=FALSE)
}

if( !opt$comparisons %in% c('POPS', 'SUPERPOPS')) {
  print_help(opt_parser)
  stop("Comparison must either be 'POPS' or 'SUPERPOPS'.", call.=FALSE)
}

chr8.vcf = opt$chr8_data
chrX.vcf = opt$chrX_data
tempdir = opt$temp_directory
pop_panel = opt$pop_data

# Read pop data from panel
all_samples <- read.table(pop_panel, header=T)
pop_comparisons <- if (opt$comparisons == 'SUPERPOPS') {
  pop_combos <- as.data.frame(t(combn(as.character(unique(all_samples$super_pop)), 2)), stringsAsFactors = F)
  if ( is.na(opt$num_samples) ) {opt$num_samples = 150}
} else if (opt$comparisons == 'POPS') {
  pop_combos <- as.data.frame(t(combn(as.character(unique(all_samples$pop)), 2)), stringsAsFactors = F)
  if ( is.na(opt$num_samples) ) {opt$num_samples = 30}
}

# Read chrX into GDS file for filtering
chrX.seqfile.path = seqVCF2GDS(chrX.vcf, paste(tempdir, 'chrX_seqGDS.gds', sep='/'))
chrX.seqfile = seqOpen(chrX.seqfile.path)

# PAR1 -----
seqSetFilterChrom(chrX.seqfile, 'X', from.bp=0, to.bp = 2699520)
seqGDS2SNP(chrX.seqfile, paste(tempdir, 'PAR1_snpGDS.gds', sep='/'))
seqResetFilter(chrX.seqfile)

# nonPAR ----
seqSetFilterChrom(chrX.seqfile, c('X', 'X'), from.bp=c(2699520, 93193855), to.bp = c(88193855, 154931044))
seqGDS2SNP(chrX.seqfile, paste(tempdir, 'nonPAR_snpGDS.gds', sep='/'))
seqResetFilter(chrX.seqfile)

# XTR ----
seqSetFilterChrom(chrX.seqfile, 'X', from.bp=88193855, to.bp = 93193855)
seqGDS2SNP(chrX.seqfile, paste(tempdir, 'XTR_snpGDS.gds', sep='/'))
seqClose(chrX.seqfile)

snpgdsVCF2GDS(chr8.vcf, paste(tempdir, 'chr8_snpGDS.gds', sep='/'))

main_df = NULL
for (i in seq(1:length(pop_combos$V1))) {
  pop1 = pop_combos$V1[i]; pop2 = pop_combos$V2[i]
  regions = c('PAR1', 'nonPAR', 'XTR', 'chr8')
  for (region in regions) {
    fileName = paste(region, '_snpGDS.gds', sep='')
    filePath = paste(tempdir, fileName, sep='/')
    genofile = snpgdsOpen(filePath)
    # use differnt columns if analyzing super pops versus pops
    # the following blocks run in either case
    if (opt$comparisons == 'SUPERPOPS') {
      if (opt$num_samples == 0) {
        pop_subset <- all_samples %>% subset(super_pop %in% c(pop1, pop2)) %>% 
          subset(gender=='female') %>% pull(sample) %>% as.character()
      } else {
        pop1_subset <- all_samples %>% subset(super_pop %in% c(pop1)) %>% subset(gender=='female') %>% 
          pull(sample) %>% sample(opt$num_samples, replace = F) %>% as.character()
        pop2_subset <- all_samples %>% subset(super_pop %in% c(pop2)) %>% subset(gender=='female') %>% 
          pull(sample) %>% sample(opt$num_samples, replace = F) %>% as.character()
        pop_subset <- c(pop1_subset, pop2_subset)
      }
      sample_df <- all_samples %>% subset(sample %in% pop_subset)
      sample_subset <- sample_df %>% pull(sample) %>% as.character()
      pop.info <- sample_df %>% pull(super_pop) %>% as.character() %>% as.factor()
    } else if (opt$comparisons == 'POPS') {
      if (opt$num_samples == 0) {
        pop_subset <- all_samples %>% subset(pop %in% c(pop1, pop2)) %>% 
          subset(gender=='female') %>% pull(sample) %>% as.character()
      } else {
        pop1_subset <- all_samples %>% subset(pop %in% c(pop1)) %>% subset(gender=='female') %>% 
          pull(sample) %>% sample(opt$num_samples, replace = F) %>% as.character()
        pop2_subset <- all_samples %>% subset(pop %in% c(pop2)) %>% subset(gender=='female') %>% 
          pull(sample) %>%sample(opt$num_samples, replace = F) %>% as.character()
        pop_subset <- c(pop1_subset, pop2_subset)
      }
      sample_df <- all_samples %>% subset(sample %in% pop_subset)
      sample_subset <- sample_df %>% pull(sample) %>% as.character()
      pop.info <- sample_df %>% pull(pop) %>% as.character() %>% as.factor()
    }
    fst = snpgdsFst(genofile, sample.id = sample_subset, population = pop.info, method="W&C84", autosome.only = FALSE)
    df = data.frame(Comparison=paste(pop1, pop2, sep='-'), region=region, fst=fst$Fst)
    main_df = rbind(main_df, df)
    snpgdsClose(genofile)
  }
}

write.csv(main_df, file=opt$output, row.names=F, quote=T)
