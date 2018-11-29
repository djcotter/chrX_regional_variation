suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(gtable))
suppressPackageStartupMessages(require(grid))
suppressPackageStartupMessages(require(tidyr))

# discrete colorblind pallete to use
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pop_levels <- c("LWK", "GWD", "ACB", "ASW", "MSL", "YRI",
                "ESN", "CEU", "IBS", "TSI", "FIN",
                "GBR", "GIH", "BEB", "ITU", "STU",
                "PJL", "JPT", "CDX", "CHB", "CHS",
                "KHV", "PUR", "CLM", "MXL", "PEL")

option_list = list(
  make_option(c('--subpops_data'), type='character', default=NULL,
              help="path to preprepared subpops data file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=10.0,
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

AFR_pops <- c("LWK", "GWD", "ACB", "ASW", "MSL", "YRI", "ESN")
AFR_color <- rep(as.character(cbbPalette[1]), length(AFR_pops))
EUR_pops <- c("CEU", "IBS", "TSI", "FIN", "GBR")
EUR_color <- rep(as.character(cbbPalette[2]), length(EUR_pops))
SAS_pops <- c("GIH", "BEB", "ITU", "STU", "PJL")
SAS_color <- rep(as.character(cbbPalette[3]), length(SAS_pops))
EAS_pops <- c("JPT", "CDX", "CHB", "CHS", "KHV")
EAS_color <- rep(as.character(cbbPalette[4]), length(EAS_pops))
AMR_pops <- c("PUR", "CLM", "MXL", "PEL")
AMR_color <- rep(as.character(cbbPalette[5]), length(AMR_pops))

# begin script --------------------------------------------------------------------------------------------------------

data <- read.delim(opt$subpops_data, header=TRUE)

data <- data %>% unite(col="PAR1",PAR1,PAR1_l,PAR1_h,sep="/")
data <- data %>% unite(col="PAR2",PAR2,PAR2_l,PAR2_h,sep="/")
data <- data %>% unite(col="XTR",XTR,XTR_l,XTR_h,sep="/")
data <- data %>% unite(col="nonPAR",nonPAR,nonPAR_l,nonPAR_h,sep="/")
data <- data %>% unite(col="chr8",chr8,chr8_l,chr8_h,sep="/")

data <- data %>% gather(Region, LD, PAR1:chr8)
data <- data %>% separate(LD, c("LD", "LD_l", "LD_h"), sep="/", convert=TRUE)
data$Region <- factor(data$Region, c("PAR1", "XTR", "nonPAR", "PAR2", "chr8"))

p_AFR <- ggplot(subset(data, SUPERPOP=="AFR"), aes(x=Region, group=POP)) + theme_pubr() +
  geom_col(aes(y=LD, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(0,0.7)) + scale_fill_manual("", values=AFR_color) +
  geom_errorbar(aes(ymin=LD_l, ymax=LD_h), position=position_dodge(width=0.9), width=0.75, size=0.8) +
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  coord_cartesian(ylim=c(0.3, 0.7)) +
  theme(axis.title=element_blank(), axis.text=element_text(size=14, angle=45, hjust=1))

p_EUR <- ggplot(subset(data, SUPERPOP=="EUR"), aes(x=Region, group=POP)) + theme_pubr() +
  geom_col(aes(y=LD, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(0,0.7)) + scale_fill_manual("", values=EUR_color) +
  geom_errorbar(aes(ymin=LD_l, ymax=LD_h), position=position_dodge(width=0.9), width=0.75, size=0.8) +
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  coord_cartesian(ylim=c(0.3, 0.7)) +
  theme(axis.title=element_blank(), axis.text=element_text(size=14, angle=45, hjust=1))

p_EAS <- ggplot(subset(data, SUPERPOP=="EAS"), aes(x=Region, group=POP)) + theme_pubr() +
  geom_col(aes(y=LD, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(0,0.7)) + scale_fill_manual("", values=EAS_color) +
  geom_errorbar(aes(ymin=LD_l, ymax=LD_h), position=position_dodge(width=0.9), width=0.75, size=0.8) +
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  coord_cartesian(ylim=c(0.3, 0.7)) +
  theme(axis.title=element_blank(), axis.text=element_text(size=14, angle=45, hjust=1))

p_SAS <- ggplot(subset(data, SUPERPOP=="SAS"), aes(x=Region, group=POP)) + theme_pubr() +
  geom_col(aes(y=LD, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(0,0.7)) + scale_fill_manual("", values=SAS_color) +
  geom_errorbar(aes(ymin=LD_l, ymax=LD_h), position=position_dodge(width=0.9), width=0.75, size=0.8) +
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  coord_cartesian(ylim=c(0.3, 0.7)) +
  theme(axis.title=element_blank(), axis.text=element_text(size=14, angle=45, hjust=1))

temp_AMR_data <- subset(data, SUPERPOP=="AMR")
temp_AMR_data$POP <- factor(temp_AMR_data$POP, levels=c("PUR", "CLM", "MXL", "PEL"))
p_AMR <- ggplot(temp_AMR_data, aes(x=Region, group=POP)) + theme_pubr() +
  geom_col(aes(y=LD, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(0,0.7)) + scale_fill_manual("", values=AMR_color) +
  geom_errorbar(aes(ymin=LD_l, ymax=LD_h), position=position_dodge(width=0.9), width=0.75, size=0.8) +
  scale_x_discrete(labels=c("PAR1", "XTR", "chrX", "PAR2", "chr8")) +
  coord_cartesian(ylim=c(0.3, 0.7)) +
  theme(axis.title=element_blank(), axis.text=element_text(size=14, angle=45, hjust=1))

# create grobs for SUPERPOP labels
a <- text_grob("Africa", size=20, face='bold')
b <- text_grob("Europe", size=20, face='bold')
c <- text_grob("S. Asia", size=20, face='bold')
d <- text_grob("E. Asia", size=20, face='bold')
e <- text_grob("Amer.", size=20, face='bold')
gt <- gtable_row('name-row', widths = unit(c(.2,.2,.2,.2,.2), 'npc'), grobs=list(a,b,c,d,e))

p <- ggarrange(as_ggplot(a), as_ggplot(b), as_ggplot(c), as_ggplot(d), as_ggplot(e), 
               p_AFR, p_EUR, p_SAS, p_EAS, p_AMR,
               ncol=5, nrow=2, legend='none', align='v',
               heights=c(0.1,1))

blank <- rectGrob(gp=gpar(col="white"))
xAxis <- text_grob("Region", size=16, face='bold')
yAxis <- text_grob(expression(bold("Average R"^"2")), size=16, face='bold', rot=90)
gt <- ggarrange(blank, yAxis, ncol=1, nrow=2, heights= c(0.1,1))

p2 <- ggarrange(gt, p, blank, xAxis, ncol=2, nrow=2,
                heights=c(1,0.1), widths=c(0.05,1))

ggsave(plot = p2, file=opt$output, height=opt$height, width=opt$width, units=opt$units)
