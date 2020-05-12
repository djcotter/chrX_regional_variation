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

# Define filter Codes --------------------------------------------
filter_codes <- c('0kb', '1kb', '5kb', '10kb', '20kb', '50kb', '100kb') %>%
  set_names(c('filter1', 'filter2', 'filter3', 'filter4',
              'filter5', 'filter6', 'filter7'))

# Colect input files ------------------------------------------

main_df <- NULL
my_files <- dir(opt$folder, "*.bed", full.names=TRUE)
for (file in my_files) {
  filename <- tools::file_path_sans_ext(basename(file)) %>% str_split('_')
  filename <- filename[[1]]
  chr <- filename[3]
  correction <- filename[6]
  filter <- filter_codes[filename[4]]
  df <- read.table(file, sep='\t', col.names = c('region', 'start', 'stop', 
                                                 'pi', 'callable_sites', 'variants',
                                                 'ci_l', 'ci_h'),
                   stringsAsFactors = F)
  
  df$correction <- correction; df$filter <- filter
  new_df <- df %>% select(filter, region, correction, callable_sites, variants, pi) %>%
    rename(`Distance from Genes` = filter, Region=region, `Divergence Correction`=correction,
            `Callable Sites`=callable_sites, `Variant Sites`=variants, `Diversity (pi)`=pi)
  main_df = rbind(main_df, new_df)
}

main_df <- main_df %>% spread(`Divergence Correction`, `Diversity (pi)`) %>% 
  rename('Uncorrected Diversity'=uncorrected, 'Corrected Diversity (canFam3)'=`canFam3-hg19-corrected`) %>%
  select(Region, `Distance from Genes`, `Callable Sites`, `Variant Sites`, `Uncorrected Diversity`, `Corrected Diversity (canFam3)`) %>%
  mutate(Region=fct_relevel(Region, 'PAR1', 'nonPAR', 'XTR', 'PAR2', 'chrY', 'chr8'),
         `Distance from Genes`=fct_relevel(`Distance from Genes`, '0kb', '1kb', '5kb', '10kb', '20kb', '50kb', '100kb')) %>%
  arrange(Region, `Distance from Genes`)
write.csv(main_df, file=opt$output, row.names=F, quote=T)
