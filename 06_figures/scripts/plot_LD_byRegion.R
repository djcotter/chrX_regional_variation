suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))

option_list = list(
  make_option(c('--chrX'), type='character', default=NULL,
              help="path to chrX regions data"),
  make_option(c('--chr8'), type='character', default=NULL,
              help="path to chr8 data"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=5.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=5.0,
              help='height for figure'),
  make_option(c('--units'), type='character', default='in',
              help='units for the figure', metavar="['in', 'cm', 'mm']"),
  make_option(c('--maxHeight'), type='double', default=1,
              help='max height for the plot (ylim).')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$units = tolower(opt$units)

if(is.null(opt$chrX) || is.null(opt$chr8) || is.null(opt$output)) {
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
chrX_data <- read.delim(opt$chrX, header=FALSE)
chrX_data$regions <- c("PAR1", "nonPAR", "XTR", "PAR2")
chr8_data <- read.delim(opt$chr8, header=FALSE)
chr8_data$regions <- "chr8"

data <- rbind(chrX_data, chr8_data)
data$regions <- factor(data$regions, levels=c('PAR1', 'XTR', 'nonPAR', 'PAR2', 'chr8'))

p1 <- ggplot(data, aes(x=regions)) + theme_pubr() +
  geom_col(aes(y=V3), color='black') + 
  geom_errorbar(aes(ymin=V4, ymax=V5), width=0.75, size=0.8) + 
  labs(x = "Regions", y=expression(bold("Average R"^"2"))) +
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  coord_cartesian(ylim=c(0.3, 0.65)) + 
  theme(axis.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=14))

ggsave(plot = p1, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
