suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))

option_list = list(
  make_option(c('--folder'), type='character', default=NULL,
              help="path to folder with all files"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output files")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$folder) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and Output files must be specified.", call.=FALSE)
}

# Define Pop Codes --------------------------------------------
pop_codes <- c(rep('AFR', 7), rep('AMR', 4), 
               rep('EAS', 5), rep('EUR', 5), 
               rep('SAS', 5)) %>%
  set_names(c(
    "ACB", "YRI", "ASW", "ESN", "MSL", "GWD", "LWK",
    "MXL", "PUR", "CLM", "PEL",
    "CHB", "JPT", "CHS", "CDX", "KHV",
    "CEU", "TSI", "FIN", "GBR", "IBS",
    "GIH", "PJL", "BEB", "STU", "ITU"
  ))

# Colect input files ------------------------------------------

main_df <- NULL
my_files <- dir(opt$folder, "*.bed", full.names=TRUE)
for (file in my_files) {
  filename <- tools::file_path_sans_ext(basename(file)) %>% str_split('_')
  filename <- filename[[1]]
  pop <- filename[1]
  superpop <- pop_codes[pop]
  df <- read.table(file, sep='\t', col.names = c('region', 'start', 'stop', 
                                                 'pi', 'callable_sites', 'variants',
                                                 'ci_l', 'ci_h', 'pVal'),
                   stringsAsFactors = F)
  
  df$pop <- pop; df$superpop <- superpop
  new_df <- df %>% select(superpop, pop, region, pi, pVal) %>%
    rename(`Super Population` = superpop, `Population`=pop,
          `Region`=region, `Diversity (pi)`=pi, `P Value (vs nonPAR)`=pVal)
  main_df = rbind(main_df, new_df)
}

main_df <- main_df %>% arrange(`Super Population`, Population, Region)
write.csv(main_df, file=opt$output, row.names=F, quote=T)
