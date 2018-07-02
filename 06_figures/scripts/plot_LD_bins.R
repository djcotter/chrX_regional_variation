suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL,
              help="path to LD_data input file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output files"),
  make_option(c('--winSize'), type='integer', default=100000,
              help='size of the window included in the LD_data file.'),
  make_option(c('--zoom'), type='integer', default=0,
              help='number of megabases that should be inluded in the plot'),
  make_option(c('--width'), type='double', default=12.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=7.0,
              help='height for figure'),
  make_option(c('--units'), type='character', default='in',
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
  stop("Units must be provided as 'in', 'mm', or 'cm'.", call.=FALSE)
}

if( !(opt$width > 0 && opt$height > 0) ){
  print_help(opt_parser)
  stop("Width and Height must be a number greater than 0.")
}

LD_data <- read.csv(opt$input, header = FALSE, sep = "\t")

LD_data$adj_position <- sapply(LD_data$V1, function(x){(((x+(x+opt$winSize))/2)/1000000)})

LD_data$region <- sapply(LD_data$Adj_position, function(x){if(x<=2.699){"PAR1"}else if(x >= 88.1 & x <= 93.1){"XTR"}else if(x>=154.9){"PAR2"}else{"nonPAR"}})
LD_data$region <- factor(LD_data$region, levels=c("PAR1", "nonPAR", "XTR", "PAR2"))

p1 = ggplot(LD_data, aes(x = adj_position))
p1 = p1 + geom_point(aes(y=V3, colour=region)) + geom_errorbar(aes(ymin=V4, ymax=V5)) 
p1 = p1 + scale_color_manual(values=c("PAR1" = "red", "nonPAR" = "black", "XTR" = "blue","PAR2" =  "red"))
p1 = p1 + labs(list(x = 'Postion (Mb)', y = expression(paste('Average ', R^2)), colour = "Region"))
p1 = p1 + theme(axis.title=element_text(size=12), legend.title=element_text(face='bold', size=10), legend.text=element_text(size=8))

if(opt$zoom > 0) {
  p1 = p1 + coord_cartesian(xlim=c(0,opt$zoom))
}

ggsave(plot=p2, file=opt$output, width=opt$width, height=opt$height, units=opt$units)
