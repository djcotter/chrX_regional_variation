suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL,
              help="path to input file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('-c', '--chromosome'), type='character', default=NULL,
              help='identifier of the chromosome to be plotted'),
  make_option(c('--width'), type='double', default=12.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=7.0,
              help='height for figure'),
  make_option(c('-u', '--units'), type='character', default='in',
              help='units for the figure', metavar="['in', 'cm', 'mm']")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$units = tolower(opt$units)

if(is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and Output files must be specified.", call.=FALSE)
}

if( !(opt$units == 'in' || opt$units == 'mm' || opt$units == 'cm') ) {
  print_help(opt_parser)
  stop("Units must be provided as 'in', 'mm', or 'cm'", call.=FALSE)
}

# declare window size for use plotting values
win_size <- 1000000 # 1 Mb
win_type <- 'Mb' # 1 Mb

# data frame should consist of cols: chr, start, end, value, called_sites, snp_count
df <- read.delim(opt$input, header=FALSE)
df$position <- sapply(df$V2, function(x){((x+(x+100000))/2)/win_size})
df$pi <- ifelse(df$V5 == 'NA', NA, df$V4)

if(opt$chromosome == 'chrX') {
  df$colors <- sapply(df$position, function(x){
    if(x<=2.699520){'red'}
    else if(x<= 93.193855 && x>= 88.193855){'blue'}
    else if(x >= 154.931044){'red'}
    else{'black'}
  })
}else if(opt$chromosome == 'chrY') {

} else {
  df$colors <- sapply(df$position, function(x){'black'})
}

# create the ggplot
p1 <- ggplot(df, aes(x=position, y=pi))

p1 <- p1 + geom_point(col=df$colors) +
  labs(list(x='Position (' + win_type +')',
            y=expression(paste('Diversity (', pi, ')')))) +
  theme(axis.title=element_text(size=14))
ggsave(file=opt$output, height=opt$height, width=opt$width, units=opt$units)
