suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(gtable))

# discrete colorblind pallete to use
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pop_levels <- c("LWK", "GWD", "ACB", "ASW", "MSL", 
                "YRI", "ESN", "CEU", "IBS", "TSI", 
                "FIN", "GBR", "GIH", "BEB", "ITU", 
                "STU", "PJL", "JPT", "CDX", "CHB", 
                "CHS", "KHV", "PUR", "CLM", "MXL", 
                "PEL")

option_list = list(
  make_option(c('--subpops_data'), type='character', default=NULL,
              help="path to preprepared subpops data file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=7.0,
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
data <- read.delim(opt$subpops_data, header=TRUE)

#prepare new diversity ratio columns
data$PAR1_A <- data$PAR1 / data$chr8
data$PAR1_A_l <- data$PAR1_l / data$chr8_l
data$PAR1_A_h <- data$PAR1_h / data$chr8_h

data$X_A <- data$nonPAR / data$chr8
data$X_A_l <- data$nonPAR_l / data$chr8_l
data$X_A_h <- data$nonPAR_h / data$chr8_h

data$XTR_A <- data$XTR / data$chr8
data$XTR_A_l <- data$XTR_l / data$chr8_l
data$XTR_A_h <- data$XTR_h / data$chr8_h

data$SUPERPOP <- factor(data$SUPERPOP, levels=c('AFR', 'EUR', 'SAS', 'EAS', 'AMR'))
data$POP <- factor(data$POP, levels=pop_levels)

p_PAR1 = ggplot(data, aes(x=POP, y=PAR1_A, fill=SUPERPOP)) + geom_col(color='black') + theme_pubr() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  labs(x='Population', y=expression(bold("PAR1"[pi] / "A"[pi]))) + 
  geom_errorbar(aes(ymin=PAR1_A_l, ymax=PAR1_A_h), width=0.5, size=0.8) + 
  coord_cartesian(ylim=c(0,1.35)) + geom_hline(yintercept=1.0) + 
  geom_hline(yintercept=0.75, linetype='dashed') + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_text(size=16, face="bold")) +
  scale_fill_manual(values=cbbPalette) + theme(plot.margin = unit(c(1,1,1,1), 'mm'))


#data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$XTR_A)])
p_XTR = ggplot(data, aes(x=POP, y=XTR_A, fill=SUPERPOP)) + geom_col(color='black') + theme_pubr() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  labs(fill='Super\nPopulation', x='Population', y=expression(bold("XTR"[pi] / "A"[pi]))) + 
  geom_errorbar(aes(ymin=XTR_A_l, ymax=XTR_A_h), width=0.5, size=0.8) + 
  coord_cartesian(ylim=c(0,1.35)) + geom_hline(yintercept=1.0) + 
  geom_hline(yintercept=0.75, linetype='dashed') + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_text(size=16, face="bold")) +
  scale_fill_manual(values=cbbPalette) + theme(plot.margin = unit(c(1,1,1,1), 'mm'))

#data$POP <- factor(data$POP, levels=data$POP[order(data$SUPERPOP, data$X_A)])
p_X = ggplot(data, aes(x=POP, y=X_A, fill=SUPERPOP)) + geom_col(color='black') + theme_pubr() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  labs(fill='Super\nPopulation', x='Population', y=expression(bold("X"[pi] / "A"[pi]))) + 
  geom_errorbar(aes(ymin=X_A_l, ymax=X_A_h), width=0.5, size=0.8) + 
  coord_cartesian(ylim=c(0,1.35)) + geom_hline(yintercept=1.0) + 
  geom_hline(yintercept=0.75, linetype='dashed') + 
  theme(axis.title.x = element_text(size=16, face = "bold"),
        axis.title.y = element_text(size=16, face = "bold")) +
  scale_fill_manual(values=cbbPalette) + theme(plot.margin = unit(c(1,1,1,1), 'mm'))

a <- text_grob("Africa", size=16, face='bold')
b <- text_grob("Europe", size=16, face='bold')
c <- text_grob("S. Asia", size=16, face='bold')
d <- text_grob("E. Asia", size=16, face='bold')
e <- text_grob("Amer.", size=16, face='bold')

gt <- gtable_row('name-row', widths = unit(c(7/26,5/26,5/26,5/26,4/26), 'npc'), grobs=list(a,b,c,d,e))

p1 = ggarrange(as_ggplot(gt), p_PAR1, p_XTR, p_X, ncol=1, nrow=4, align='v',
               common.legend = TRUE, legend='none',
               heights = c(0.2,1,1,1.3))

ggsave(plot = p1, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
