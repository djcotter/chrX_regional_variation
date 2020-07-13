library(tidyverse)
library(optparse)

option_list = list(
  make_option(c('--input'), type='character', default=NULL,
              help="path to input file"),
  make_option(c('--output'), type='character', default=NULL,
              help="path to windowed diversity file")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and Output files must be specified.", call.=FALSE)
}

region_lengths <- data.frame(region=c("PAR1", "nonPAR", "XTR", "PAR2", "chr8"), 
                            region_lengths=c(2639519, 147231524, 5000000, 329516, 146364022))

# "PAR1": [60001, 2699520]
# "nonPAR1": [2699520, 88193855]
# "XTR": [88193855, 93193855]
# "nonPAR2": [93193855, 154931044]
# "PAR2": [154931044, 155260560]

df <- read.table(opt$input, col.names = c('chr', 'start', 'stop')) %>% 
  subset(chr == 'chrX' | chr=='chr8')

df_chr8 <- df %>% subset(chr=='chr8') %>% 
  mutate(label='chr8')

df_chrX <- df %>% subset(chr=='chrX') %>% 
  mutate(label=ifelse(start >= 60001 & stop <= 2699520, 'PAR1',
                      ifelse(start >= 2699520 & stop <= 88193855, 'nonPAR',
                             ifelse(start >= 88193855 & stop <= 93193855, 'XTR',
                                    ifelse(start >= 93193855 & stop <= 154931044, 'nonPAR',
                                           ifelse(start >= 154931044 & stop <= 155260560, 'PAR2', 'NA'))))))

df_chrX_straddle <- df_chrX %>% subset(label=='NA') %>% 
  mutate(start=ifelse(start >60001 & start < 2699520, paste(start, 2699520, sep=":"),
                      ifelse(start > 2699520 & start < 88193855, paste(start, 88193855, sep=":"),
                             ifelse(start > 88193855 & start < 93193855, paste(start, 93193855, sep=":"),
                                    ifelse(start > 93193855 & start < 154931044, paste(start, 154931044, sep=":"),
                                           ifelse(start > 154931044 & start < 155260560, paste(start, 155260560, sep=":"), 'NA')))))) %>%
  mutate(stop=ifelse(stop > 2699520 & stop < 88193855, paste(2699520, stop, sep=":"),
                     ifelse(stop > 88193855 & stop < 93193855, paste(88193855, stop, sep=":"),
                            ifelse(stop > 93193855 & stop < 154931044, paste(93193855, stop, sep=":"),
                                   ifelse(stop > 154931044 & stop < 155260560, paste(154931044, stop, sep=":"), 'NA')))))

df_chrX_straddle <- rbind(data.frame(chr='chrX', pos=df_chrX_straddle$start), data.frame(chr='chrX', pos=df_chrX_straddle$stop)) %>%
  separate(pos, sep=':', into=c('start', 'stop'), convert = TRUE) %>% 
  mutate(label=ifelse(start >= 60001 & stop <= 2699520, 'PAR1',
                      ifelse(start >= 2699520 & stop <= 88193855, 'nonPAR',
                             ifelse(start >= 88193855 & stop <= 93193855, 'XTR',
                                    ifelse(start >= 93193855 & stop <= 154931044, 'nonPAR',
                                           ifelse(start >= 154931044 & stop <= 155260560, 'PAR2', 'NA'))))))

df_chrX <- df_chrX %>% subset(label != 'NA') %>% bind_rows(df_chrX_straddle)

main_df <- rbind(df_chr8, df_chrX) %>% 
  select(-chr) %>% 
  rename(region=label) %>% 
  mutate(length=stop-start) %>% 
  group_by(region) %>% 
  summarise(totalLength=sum(length)) %>%
  inner_join(region_lengths) %>%
  mutate(percent = totalLength/region_lengths) %>%
  rename(Region=region, `Region Length`=region_lengths,
         `Total Length of Genes`=totalLength, `Fraction of Region that is Coding`=percent) %>%
  select(1,3,2,4)

write.csv(main_df, file=opt$output, row.names = F)
