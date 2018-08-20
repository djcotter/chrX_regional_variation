suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(gtable))

# discrete colorblind pallete to use
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pop_levels <- c("LWK", "GWD", "ACB", "ASW", "YRI", 
                "ESN", "CEU", "IBS", "TSI", "FIN",
                "GBR", "GIH", "BEB", "ITU", "STU",
                "PJL", "JPT", "CDX", "CHB", "CHS", 
                "KHV", "PUR", "CLM", "MXL", "PEL")

option_list = list(
  make_option(c('--subpops_data'), type='character', default=NULL,
              help="path to preprepared subpops data file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=7.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=4.0,
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
data$X_A <- data$nonPAR / data$chr8
data$XTR_A <- data$XTR / data$chr8
data$PAR1_A <- data$PAR1 / data$chr8

# plot data on a line
data$SUPERPOP <- factor(data$SUPERPOP, levels=c('AFR', 'EUR', 'SAS', 'EAS', 'AMR'))
data$POP <- factor(data$POP, levels=pop_levels)

p1 = ggplot(data, aes(x=POP, color=SUPERPOP)) + geom_point(aes(y=X_PAR1, shape="X_PAR1"), size=2, stroke=1.5) + theme_pubr() +
  coord_cartesian(ylim=c(0.6,1.2)) + labs(color='Super\nPopulation', x='Population', y='Relative Ratios') + 
  geom_point(aes(y=X_A, shape="X_A"), size=2, stroke=1.5) + geom_hline(yintercept = 1) + 
  scale_color_manual(values=cbbPalette) +
  scale_shape_manual("Ratio", breaks=c("X_A", "X_PAR1"), 
                     labels=c(expression(bold("X"[pi] / "A"[pi])), expression(bold("X"[pi] / "PAR"[pi]))), 
                     values=c(4,1)) + theme(legend.text.align = 0) + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  guides(color=FALSE, shape=guide_legend(direction='vertical',title=NULL)) + 
  theme(legend.position=c(0.8,0.85), legend.text=element_text(size=14,face='bold'),
        legend.background=element_rect(fill = 'lightgrey')) +
  theme(axis.title.x = element_text(size=16, face = "bold"),
        axis.title.y = element_text(size=16, face = "bold"))

a <- text_grob("Africa", size=16, face='bold')
b <- text_grob("Europe", size=16, face='bold')
c <- text_grob("S. Asia", size=16, face='bold')
d <- text_grob("E. Asia", size=16, face='bold')
e <- text_grob("Amer.", size=16, face='bold')

gt <- gtable_row('name-row', widths = unit(c(6/25,5/25,5/25,5/25,4/25), 'npc'), grobs=list(a,b,c,d,e))

p <- ggarrange(as_ggplot(gt), p1, ncol=1, nrow=2, align='v', heights = c(0.2,1))


ggsave(plot = p, file=opt$output, height=opt$height, width=opt$width, units=opt$units)

# -----------------------------------------------------------------
# p2 = ggplot(data, aes(x=POP, y=X_A, color=SUPERPOP)) + geom_point() + geom_hline(yintercept = 1)
# p2 = p2 + coord_cartesian(ylim=c(0.6,1.2)) + theme(axis.text.x = element_text(angle=45,hjust=1))
# p2 = p2 + labs(color='Super\nPopulation', x='Population', y=expression("Relative X"[pi] / "A"[pi]))
# 
# data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$XTR_A)])
# p3 = ggplot(data, aes(x=POP, y=XTR_A, color=SUPERPOP)) + geom_point() + geom_hline(yintercept = 1)
# p3 = p3 + coord_cartesian(ylim=c(0.6,1.2)) + theme(axis.text.x = element_text(angle=45,hjust=1))
# p3 = p3 + labs(color='Super\nPopulation', x='Population', y=expression("Relative XTR"[pi] / "A"[pi]))
# 
# 
# data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$PAR1_A)])
# p4 = ggplot(data, aes(x=POP, y=PAR1_A, color=SUPERPOP)) + geom_point() + geom_hline(yintercept = 1)
# p4 = p4 + coord_cartesian(ylim=c(0.6,1.2)) + theme(axis.text.x = element_text(angle=45,hjust=1))
# p4 = p4 + labs(color='Super\nPopulation', x='Population', y=expression("Relative PAR1"[pi] / "A"[pi]))
