suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))

option_list = list(
  make_option(c('--subpops_data'), type='character', default=NULL,
              help="path to preprepared subpops data file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=12.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=12.0,
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
data <- read.delim(opt$subpops_data, header=TRUE)

#prepare new diversity ratio columns
data$PAR1_A <- data$PAR1 / data$chr8
data$X_A <- data$nonPAR / data$chr8
data$XTR_A <- data$XTR / data$chr8


data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$PAR1_A)])
p_PAR1 = ggplot(data, aes(x=POP, y=PAR1_A, fill=SUPERPOP)) + geom_col()
p_PAR1 = p_PAR1 + theme(axis.text.x = element_text(angle=45,hjust=1))
p_PAR1 = p_PAR1 + labs(fill='Super\nPopulation', x='Population', y=expression("PAR1"[pi] / "A"[pi]))

data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$X_A)])
p_X = ggplot(data, aes(x=POP, y=X_A, fill=SUPERPOP)) + geom_col()
p_X = p_X + theme(axis.text.x = element_text(angle=45,hjust=1))
p_X = p_X + labs(fill='Super\nPopulation', x='Population', y=expression("X"[pi] / "A"[pi]))

data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$XTR_A)])
p_XTR = ggplot(data, aes(x=POP, y=XTR_A, fill=SUPERPOP)) + geom_col()
p_XTR = p_XTR + theme(axis.text.x = element_text(angle=45,hjust=1))
p_XTR = p_XTR + labs(fill='Super\nPopulation', x='Population', y=expression("XTR"[pi] / "A"[pi]))

p1 = ggarrange(p_X, p_PAR1, p_XTR, ncol=1, nrow=3, align='v', common.legend = TRUE, legend='right')

ggsave(plot = p1, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
