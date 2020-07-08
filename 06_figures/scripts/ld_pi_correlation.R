suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))

option_list = list(
  make_option(c('--LD'), type='character', default=NULL,
              help="path to average LD by window file"),
  make_option(c('--diversity'), type='character', default=NULL,
              help="path to windowed diversity file"),
  make_option(c('--LD2'), type='character', default=NULL,
              help="path to average LD by window file"),
  make_option(c('--diversity2'), type='character', default=NULL,
              help="path to windowed diversity file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output files"),
  make_option(c('--filter'), type='character', default=NULL,
              help="how many kb were filtered with the provided files"),
  make_option(c('--winSize'), type='integer', default=100000,
              help='size of the window included in the LD_data file.'),
  make_option(c('--width'), type='double', default=8.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=8.0,
              help='height for figure'),
  make_option(c('--units'), type='character', default='in',
              help='units for the figure', metavar="['in', 'cm', 'mm']"),
  make_option(c('--maxHeight'), type='double', default=1,
              help='max height for the plot (ylim).')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$units = tolower(opt$units)

if(is.null(opt$LD) || is.null(opt$diversity) || is.null(opt$output) || is.null(opt$filter)) {
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
label2 <- paste("P ==", summary(model)$coefficients[2,4], sep=" ")

# LD2 <- read.delim(opt$LD2, header=FALSE)
# pi2 <- read.delim(opt$diversity2, header=FALSE)
# 
# pi2$adj_position <- sapply(pi2$V2, function(x){(((x+(x+opt$winSize))/2)/1000000)})
# pi2$region <- sapply(pi2$adj_position, function(x){if(x<=2.699){"PAR1"}else if(x >= 88.1 & x <= 93.1){"XTR"}else if(x>=154.9){"PAR2"}else{"nonPAR"}})
# 
# data2 <- data.frame (chr="chrX", adj_position=pi2$adj_position, region=pi2$region, R_squared=LD2$V3, pi=pi2$V4)
# data2 <- na.omit(data2)
# 
# model2 <- lm(R_squared ~ pi, data=data2)
# label2 <- paste("R^2 ==", summary(model2)$r.squared, sep=" ")

group.colors = c(PAR1="red", nonPAR="black", XTR="blue", PAR2="red")
group.shapes = c(PAR1=17, nonPAR=20, XTR=15, PAR2=19)
group.sizes = c(PAR1=3, nonPAR=2, XTR=3, PAR2=3)

p1 <- ggplot(data, aes(x=R_squared, y=pi)) + theme_pubr() + 
  geom_point(aes(shape=region, color=region, size=region)) + 
  scale_colour_manual(name="Region", 
                      breaks=c("PAR1", "nonPAR", "XTR", "PAR2"), 
                      values=group.colors) + 
  scale_shape_manual(name="Region",
                     breaks=c("PAR1", "nonPAR", "XTR", "PAR2"),
                     values=group.shapes) + 
  scale_size_manual(name="Region",
                    breaks=c("PAR1", "nonPAR", "XTR", "PAR2"),
                    values=group.sizes) +
  xlab(bquote(bold('Linkage Disequilibrium (Average '* R^2*')'))) + 
  ylab(bquote(bold('Diversity ('*pi*')'))) + 
  coord_cartesian(xlim=c(0.35,0.85), ylim=c(0,0.008)) +
  annotate("text", x=0.65, y=0.006, 
           label=label1, parse=TRUE, fontface=2, size=5.5) +
  annotate("text", x=0.65, y=0.0055,
           label=label2, parse=TRUE, fontface=2, size=5.5) +
  geom_smooth(aes(x=data$R_squared, y=data$pi), 
              method=lm, color="green", size=1) +
  theme(axis.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=14),
        legend.title = element_text(size=16, face='bold'),
        legend.text = element_text(size=14, face='bold')) +
  labs(title=paste(opt$filter, "filtered from genes", sep=" ")) + 
  theme(title=element_text(size=18, face='bold')) 

p1
# p2 <- ggplot(data2, aes(x=R_squared, y=pi)) + theme_pubr() +
#   geom_point(aes(shape=region, color=region, size=region)) +
#   scale_colour_manual(name="Region",
#                       breaks=c("PAR1", "nonPAR", "XTR", "PAR2"),
#                       values=group.colors) +
#   scale_shape_manual(name="Region",
#                      breaks=c("PAR1", "nonPAR", "XTR", "PAR2"),
#                      values=group.shapes) +
#   scale_size_manual(name="Region",
#                     breaks=c("PAR1", "nonPAR", "XTR", "PAR2"),
#                     values=group.sizes) +
#   xlab(bquote(bold('Linkage Disequilibrium (Average '* R^2*')'))) +
#   ylab(bquote(bold('Diversity ('*pi*')'))) +
#   annotate("text", x=0.65, y=max(data2$pi)-mean(data2$pi)/2,
#            label=label2, parse=TRUE, fontface=2, size=7) +
#   geom_smooth(aes(x=data2$R_squared, y=data2$pi),
#               method=lm, color="green", size=1) +
#   theme(axis.title = element_text(size=16, face="bold"),
#         axis.text = element_text(size=14),
#         legend.title = element_text(size=16, face='bold'),
#         legend.text = element_text(size=14, face='bold')) +
#   labs(title="10 kb from Genes") +
#   theme(title=element_text(size=18, face='bold'))

# p <- ggarrange(p1, p2, ncol=1, nrow=2, align='v', common.legend = TRUE, legend='right')

p <- p1

ggsave(plot = p, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
