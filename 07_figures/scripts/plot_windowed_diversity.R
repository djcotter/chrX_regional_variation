suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))

option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL,
              help="path to input file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('-c', '--chrom'), type='character', default=NULL,
              help='identifier of the chromosome to be plotted'),
  make_option(c('--width'), type='double', default=7,
              help='width for figure'),
  make_option(c('--height'), type='double', default=3,
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
p1 <- ggplot(df, aes(x=position, y=pi)) + theme_pubr()

if(opt$chrom == 'chrX') {
  p1 <- p1 +  geom_rect(xmin=0, xmax=2.899520, ymin=0, ymax=opt$maxHeight, fill='grey', alpha=0.3) +
    geom_rect(xmin=88.193855, xmax=93.993855, ymin=0, ymax=opt$maxHeight, fill='grey', alpha=0.3) +
    geom_rect(xmin=154.131044, xmax=156.341, ymin=0, ymax=opt$maxHeight, fill='grey', alpha=0.3) +
    labs(x=paste('X Chromosome Position (',win_type,')', sep=""))
} else {
  p1 <- p1 + labs(x=paste(opt$chrom,' Position (',win_type,')', sep=""))
}

p1 <- p1 + ylim(0,opt$maxHeight) +
  geom_point(col=df$colors) +
  labs(y=expression(bold(paste('Diversity (', pi, ')')))) +
  theme(axis.title=element_text(size=16, face="bold"))

ggsave(file=opt$output, height=opt$height, width=opt$width, units=opt$units)
