suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('--LD'), type='character', default=NULL,
              help="path to average LD by window file"),
  make_option(c('--diversity'), type='character', default=NULL,
              help="path to windowed diversity file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output files"),
  make_option(c('--winSize'), type='integer', default=100000,
              help='size of the window included in the LD_data file.'),
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

if(is.null(opt$LD) || is.null(opt$diversity) || is.null(opt$output)) {
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

LD <- read.delim(opt$LD, header=FALSE)
pi <- read.delim(opt$diversity, header=FALSE)

pi$adj_position <- sapply(pi$V2, function(x){(((x+(x+opt$winSize))/2)/1000000)})
pi$region <- sapply(pi$adj_position, function(x){if(x<=2.699){"PAR1"}else if(x >= 88.1 & x <= 93.1){"XTR"}else if(x>=154.9){"PAR2"}else{"nonPAR"}})

data <- data.frame (chr="chrX", adj_position=pi$adj_position, region=pi$region, R_squared=LD$V3, pi=pi$V4)
data <- na.omit(data)

model <- lm(R_squared ~ pi, data=data)
label1 <- paste("R^2 ==", summary(model)$r.squared, sep=" ")
group.colors = c(PAR1="red", nonPAR="black", XTR="blue", PAR2="red")

p1 = ggplot(data, aes(x=R_squared, y=pi))  
p1 = p1 + geom_point(aes(color=data$region)) + scale_colour_manual(name="Region", breaks=c("PAR1", "nonPAR", "XTR", "PAR2"), values=group.colors) 
p1 = p1 + xlab(bquote('Linkage Disequilibrium (Average '* R^2*')')) + ylab(bquote('Diversity ('*pi*')')) + annotate("text", x=0.65, y=max(data$pi)-mean(data$pi)/2, label=label1, parse=TRUE, fontface=2, size=5)
p1 = p1 + geom_smooth(aes(x=data$R_squared, y=data$pi), method=lm, color="green")

ggsave(plot = p1, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
