suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL,
              help="path to input file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('-c', '--chrom'), type='character', default=NULL,
              help='identifier of the chromosome to be plotted'),
  make_option(c('--width'), type='double', default=12.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=7.0,
              help='height for figure'),
  make_option(c('--units'), type='character', default='in',
              help='units for the figure', metavar="['in', 'cm', 'mm']"),
  make_option(c('--maxHeight'), type='double', default=1,
              help='max height for the plot (ylim).')
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
  stop("Units must be provided as 'in', 'mm', or 'cm'.", call.=FALSE)
}

if( !(opt$width > 0 && opt$height > 0) ){
  print_help(opt_parser)
  stop("Width and Height must be a number greater than 0.")
}

# declare window size for use plotting values
win_size <- 1000000 # 1 Mb
win_type <- 'Mb' # 1 Mb

# data frame should consist of cols: chr, start, end, value, called_sites, snp_count
df <- read.delim(opt$input, header=FALSE)
df$position <- sapply(df$V2, function(x){((x+(x+100000))/2)/win_size})
df$pi <- ifelse(df$V5 == 'NA', NA, df$V4)

if(opt$chrom == 'chrX') {
  df$colors <- sapply(df$position, function(x){
    if(x<=2.699520){'red'}
    else if(x<= 93.193855 && x>= 88.193855){'blue'}
    else if(x >= 154.931044){'red'}
    else{'black'}
    })
  max_height <- 0.0045
}else if(opt$chrom == 'chrY') {
  df$colors <- sapply(df$position, function(x){'black'})
  max_height <- 0.0003
} else {
  df$colors <- sapply(df$position, function(x){'black'})
  max_height <- 0.0071
}

#max_height = 0.06

# create the ggplot
p1 <- ggplot(df, aes(x=position, y=pi)) + ylim(0,opt$maxHeight)
#p1 <- p1 + geom_errorbar(aes(ymin=V7, ymax=V8))
p1 <- p1 + geom_point(col=df$colors) +
  labs(list(x=paste('Position (', win_type, ')'),
            y=expression(paste('Diversity (', pi, ')')))) +
  theme(axis.title=element_text(size=14))
ggsave(file=opt$output, height=opt$height, width=opt$width, units=opt$units)
