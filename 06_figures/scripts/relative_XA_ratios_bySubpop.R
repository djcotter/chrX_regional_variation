suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('--subpops_data'), type='character', default=NULL,
              help="path to preprepared subpops data file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
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

if(is.null(opt$subpops_data) || is.null(opt$output)) {
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

# read in the data
sample_data <- read.delim(opt$subpops_data, header=TRUE)

# divide by most diverse pop
POP_max <- levels(factor(sample_data$POP, levels=sample_data$POP[order(sample_data$chr8, sample_data$PAR1, sample_data$nonPAR, sample_data$XTR, sample_data$chrY, decreasing = TRUE)]))[1]

denom <- subset(sample_data,POP==POP_max)
data <- subset(sample_data,POP!=POP_max)

data$PAR1 <- data$PAR1 / denom$PAR1
data$nonPAR <- data$nonPAR / denom$nonPAR
data$XTR <- data$XTR / denom$XTR
data$PAR2 <- data$PAR2 / denom$PAR2
data$chrY <- data$chrY / denom$chrY
data$chr8 <- data$chr8 / denom$chr8

# prepare new diversity ratio columns
data$X_PAR1 <- data$nonPAR / data$PAR1

# plot data on a line
data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$X_PAR1)])
p1 = ggplot(data, aes(x=POP, y=X_PAR1, color=SUPERPOP)) + geom_point() + geom_hline(yintercept = 1)
p1 = p1 + coord_cartesian(ylim=c(0.6,1.2)) + theme(axis.text.x = element_text(angle=45,hjust=1))
p1 = p1 + labs(color='Super\nPopulation', x='Population', y=expression("Relative X"[pi] / "PAR"[pi]))

ggsave(plot = p1, file=opt$output, height=opt$height, width=opt$width, units=opt$units)