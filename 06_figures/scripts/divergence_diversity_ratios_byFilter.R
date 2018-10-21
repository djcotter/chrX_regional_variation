suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))

option_list = list(
  make_option(c('--divergenceRatios'), type='character', default=NULL,
              help="path to preprepared divergence ratios summary"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=6.0,
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

if(is.null(opt$divergenceRatios) || is.null(opt$output)) {
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
data <- read.delim(opt$divergenceRatios, header=TRUE)

data$filter <- factor(data$filter, levels=c("0 kb", "1 kb", "5 kb", "10 kb", "20 kb"))

data$region <- factor(data$region, levels=c("PAR1", "XTR", "nonPAR", "PAR2", "chr8"))
data$region2 <- factor(data$region2, levels=c("PAR1/A", "XTR/A", "X/A", "PAR2/A", "blank"))

p1 <- ggplot(data, aes(x=region, y=divergence)) + theme_pubr() +
  geom_col(aes(fill=filter), position="dodge", color='black') + 
  labs(x="", y="Divergence", fill="Filter") + scale_fill_grey() + 
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  theme(axis.title.y = element_text(size=14, face = "bold")) +
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=10, face= "bold"))

p2 <- ggplot(data, aes(x=region2, y=div_ratio)) + theme_pubr() +
  geom_col(aes(fill=filter), position="dodge", color='black') + 
  labs(x="", y="Divergence Ratio", fill="Filter") + scale_fill_grey() + 
  scale_x_discrete("Genomic Position", labels=c("PAR1/A", "XTR/A", "chrX/A", "PAR2/A", ""),
                   limits=c("PAR1/A", "XTR/A", "X/A", "PAR2/A", "PAR1/A")) +
  geom_hline(yintercept=1) +
  theme(axis.title.x = element_text(size=14, face = "bold"),
        axis.title.y = element_text(size=14, face = "bold")) +
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=10, face= "bold"))

p <- ggarrange(p1, p2, ncol=1, nrow=2, align="v", labels=c("a)", "b)"), common.legend = TRUE, legend="right" )

ggsave(plot = p, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
