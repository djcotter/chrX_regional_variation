suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(gtable))
suppressPackageStartupMessages(require(grid))

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

option_list = list(
  make_option(c('--filter1'), type='character', default=NULL,
              help="path to filter1 ratios data"),
  make_option(c('--filter2'), type='character', default=NULL,
              help="path to filter1 ratios data"),
  make_option(c('--filter3'), type='character', default=NULL,
              help="path to filter1 ratios data"),
  make_option(c('--filter4'), type='character', default=NULL,
              help="path to filter1 ratios data"),
  make_option(c('--filter5'), type='character', default=NULL,
              help="path to filter1 ratios data"),
  make_option(c('-o1', '--output1'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('-o2', '--output2'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=14.0,
              help='width for figure'),
  make_option(c('--height'), type='double', default=10.0,
              help='height for figure'),
  make_option(c('--units'), type='character', default='in',
              help='units for the figure', metavar="['in', 'cm', 'mm']"),
  make_option(c('--maxHeight'), type='double', default=1,
              help='max height for the plot (ylim).')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$units = tolower(opt$units)

if(is.null(opt$filter1) || is.null(opt$filter2) || is.null(opt$filter3) || is.null(opt$filter4) || is.null(opt$filter5) || is.null(opt$output1) || is.null(opt$output2)) {
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

data_filter1 <- read.delim(opt$filter1, header=TRUE)
data_filter2 <- read.delim(opt$filter2, header=TRUE)
data_filter3 <- read.delim(opt$filter3, header=TRUE)
data_filter4 <- read.delim(opt$filter4, header=TRUE)
data_filter5 <- read.delim(opt$filter5, header=TRUE)
  
pop_list <- c("LWK", "GWD", "ACB", "ASW", "YRI", 
              "ESN", "CEU", "IBS", "TSI", "FIN",
              "GBR", "GIH", "BEB", "ITU", "STU",
              "PJL", "JPT", "CDX", "CHB", "CHS", 
              "KHV", "PUR", "CLM", "MXL", "PEL")

AFR_pops <- c("LWK", "GWD", "ACB", "ASW", "YRI", "ESN")
AFR_color <- rep(as.character(cbbPalette[1]), length(AFR_pops))
EUR_pops <- c("CEU", "IBS", "TSI", "FIN", "GBR")
EUR_color <- rep(as.character(cbbPalette[2]), length(EUR_pops))
SAS_pops <- c("GIH", "BEB", "ITU", "STU", "PJL")
SAS_color <- rep(as.character(cbbPalette[3]), length(SAS_pops))
EAS_pops <- c("JPT", "CDX", "CHB", "CHS", "KHV")
EAS_color <- rep(as.character(cbbPalette[4]), length(EAS_pops))
AMR_pops <- c("PUR", "CLM", "MXL", "PEL")
AMR_color <- rep(as.character(cbbPalette[5]), length(AMR_pops))

denom_pop = 'MSL'

data_filter1$X_A <- data_filter1$nonPAR / data_filter1$chr8
# data_filter1$X_A_l <- data_filter1$nonPAR_l / data_filter1$chr8_l
# data_filter1$X_A_h <- data_filter1$nonPAR_h / data_filter1$chr8_h
data_filter1$X_A_corrected <- data_filter1$X_A / subset(data_filter1,POP==denom_pop)$X_A
# data_filter1$X_A_corrected_l <- data_filter1$X_A_l / subset(data_filter1,POP==denom_pop)$X_A_l
# data_filter1$X_A_corrected_h <- data_filter1$X_A_h / subset(data_filter1,POP==denom_pop)$X_A_h
data_filter1$PAR_A <- data_filter1$PAR1 / data_filter1$chr8
data_filter1$PAR_A_corrected <- data_filter1$PAR_A / subset(data_filter1,POP==denom_pop)$PAR_A
data_filter1$XTR_A <- data_filter1$XTR / data_filter1$chr8
data_filter1$XTR_A_corrected <- data_filter1$XTR_A / subset(data_filter1,POP==denom_pop)$XTR_A
data_filter1$distance <- '0 kb'
data_filter1$distance_num <- 0

data_filter2$X_A <- data_filter2$nonPAR / data_filter2$chr8
# data_filter2$X_A_l <- data_filter2$nonPAR_l / data_filter2$chr8_l
# data_filter2$X_A_h <- data_filter2$nonPAR_h / data_filter2$chr8_h
data_filter2$X_A_corrected <- data_filter2$X_A / subset(data_filter2,POP==denom_pop)$X_A
# data_filter2$X_A_corrected_l <- data_filter2$X_A_l / subset(data_filter2,POP==denom_pop)$X_A_l
# data_filter2$X_A_corrected_h <- data_filter2$X_A_h / subset(data_filter2,POP==denom_pop)$X_A_h
data_filter2$PAR_A <- data_filter2$PAR1 / data_filter2$chr8
data_filter2$PAR_A_corrected <- data_filter2$PAR_A / subset(data_filter2,POP==denom_pop)$PAR_A
data_filter2$XTR_A <- data_filter2$XTR / data_filter2$chr8
data_filter2$XTR_A_corrected <- data_filter2$XTR_A / subset(data_filter2,POP==denom_pop)$XTR_A

data_filter2$distance <- '1 kb'
data_filter2$distance_num <- 1

data_filter3$X_A <- data_filter3$nonPAR / data_filter3$chr8
# data_filter3$X_A_l <- data_filter3$nonPAR_l / data_filter3$chr8_l
# data_filter3$X_A_h <- data_filter3$nonPAR_h / data_filter3$chr8_h
data_filter3$X_A_corrected <- data_filter3$X_A / subset(data_filter3,POP==denom_pop)$X_A
# data_filter3$X_A_corrected_l <- data_filter3$X_A_l / subset(data_filter3,POP==denom_pop)$X_A_l
# data_filter3$X_A_corrected_h <- data_filter3$X_A_h / subset(data_filter3,POP==denom_pop)$X_A_h
data_filter3$PAR_A <- data_filter3$PAR1 / data_filter3$chr8
data_filter3$PAR_A_corrected <- data_filter3$PAR_A / subset(data_filter3,POP==denom_pop)$PAR_A
data_filter3$XTR_A <- data_filter3$XTR / data_filter3$chr8
data_filter3$XTR_A_corrected <- data_filter3$XTR_A / subset(data_filter3,POP==denom_pop)$XTR_A

data_filter3$distance <- '5 kb'
data_filter3$distance_num <- 5

data_filter4$X_A <- data_filter4$nonPAR / data_filter4$chr8
# data_filter4$X_A_l <- data_filter4$nonPAR_l / data_filter4$chr8_l
# data_filter4$X_A_h <- data_filter4$nonPAR_h / data_filter4$chr8_h
data_filter4$X_A_corrected <- data_filter4$X_A / subset(data_filter4,POP==denom_pop)$X_A
# data_filter4$X_A_corrected_l <- data_filter4$X_A_l / subset(data_filter4,POP==denom_pop)$X_A_l
# data_filter4$X_A_corrected_h <- data_filter4$X_A_h / subset(data_filter4,POP==denom_pop)$X_A_h
data_filter4$PAR_A <- data_filter4$PAR1 / data_filter4$chr8
data_filter4$PAR_A_corrected <- data_filter4$PAR_A / subset(data_filter4,POP==denom_pop)$PAR_A
data_filter4$XTR_A <- data_filter4$XTR / data_filter4$chr8
data_filter4$XTR_A_corrected <- data_filter4$XTR_A / subset(data_filter4,POP==denom_pop)$XTR_A

data_filter4$distance <- '10 kb'
data_filter4$distance_num <- 10

data_filter5$X_A <- data_filter5$nonPAR / data_filter5$chr8
# data_filter5$X_A_l <- data_filter5$nonPAR_l / data_filter5$chr8_l
# data_filter5$X_A_h <- data_filter5$nonPAR_h / data_filter5$chr8_h
data_filter5$X_A_corrected <- data_filter5$X_A / subset(data_filter5,POP==denom_pop)$X_A
# data_filter5$X_A_corrected_l <- data_filter5$X_A_l / subset(data_filter5,POP==denom_pop)$X_A_l
# data_filter5$X_A_corrected_h <- data_filter5$X_A_h / subset(data_filter5,POP==denom_pop)$X_A_h
data_filter5$PAR_A <- data_filter5$PAR1 / data_filter5$chr8
data_filter5$PAR_A_corrected <- data_filter5$PAR_A / subset(data_filter5,POP==denom_pop)$PAR_A
data_filter5$XTR_A <- data_filter5$XTR / data_filter5$chr8
data_filter5$XTR_A_corrected <- data_filter5$XTR_A / subset(data_filter5,POP==denom_pop)$XTR_A

data_filter5$distance <- '20 kb'
data_filter5$distance_num <- 20

data <- rbind(data_filter1, data_filter2, data_filter3, data_filter4, data_filter5)

data2 <- subset(data, select = c('SUPERPOP', 'POP', 'distance', 'X_A', 'X_A_corrected', 'PAR_A', 'PAR_A_corrected', 'XTR_A', 'XTR_A_corrected'))
newdata <- subset(data2, POP=='')
for(TEMP_POP in pop_list) {
  temp1 <- subset(data2, POP==TEMP_POP)
  SUPER <- as.character(temp1$SUPERPOP)[1]
  temp2 <- subset(temp1, distance=='0 kb')
  temp3 <- subset(temp1, distance!='0 kb')
  temp4 <- data.frame(SUPERPOP=SUPER, POP=TEMP_POP, distance=temp3$distance,
                      X_A=temp3$X_A - temp2$X_A, 
                      X_A_corrected=temp3$X_A_corrected - temp2$X_A_corrected,
                      PAR_A=temp3$PAR_A - temp2$PAR_A,
                      PAR_A_corrected=temp3$PAR_A_corrected - temp2$PAR_A_corrected,
                      XTR_A=temp3$XTR_A - temp2$XTR_A,
                      XTR_A_corrected=temp3$XTR_A_corrected - temp2$XTR_A_corrected)
  newdata <- rbind(newdata, temp4)
}
newdata$distance <- factor(newdata$distance, levels=c('1 kb', '5 kb', '10 kb', '20 kb'))

#---------------------------------------------------------------
p_XA_AFR <- ggplot(subset(newdata, SUPERPOP=='AFR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AFR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_EUR <- ggplot(subset(newdata, SUPERPOP=='EUR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EUR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_SAS <- ggplot(subset(newdata, SUPERPOP=='SAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=SAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_EAS <- ggplot(subset(newdata, SUPERPOP=='EAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_AMR <- ggplot(subset(newdata, SUPERPOP=='AMR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AMR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))

p_XA_corrected_AFR <- ggplot(subset(newdata, SUPERPOP=='AFR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AFR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_corrected_EUR <- ggplot(subset(newdata, SUPERPOP=='EUR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EUR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_corrected_SAS <- ggplot(subset(newdata, SUPERPOP=='SAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=SAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_corrected_EAS <- ggplot(subset(newdata, SUPERPOP=='EAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XA_corrected_AMR <- ggplot(subset(newdata, SUPERPOP=='AMR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=X_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AMR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"X"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))

p_PARA_AFR <- ggplot(subset(newdata, SUPERPOP=='AFR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AFR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_EUR <- ggplot(subset(newdata, SUPERPOP=='EUR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EUR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_SAS <- ggplot(subset(newdata, SUPERPOP=='SAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=SAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_EAS <- ggplot(subset(newdata, SUPERPOP=='EAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_AMR <- ggplot(subset(newdata, SUPERPOP=='AMR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AMR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))

p_PARA_corrected_AFR <- ggplot(subset(newdata, SUPERPOP=='AFR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AFR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_corrected_EUR <- ggplot(subset(newdata, SUPERPOP=='EUR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EUR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_corrected_SAS <- ggplot(subset(newdata, SUPERPOP=='SAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=SAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_corrected_EAS <- ggplot(subset(newdata, SUPERPOP=='EAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_PARA_corrected_AMR <- ggplot(subset(newdata, SUPERPOP=='AMR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=PAR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AMR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"PAR1"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))

p_XTRA_AFR <- ggplot(subset(newdata, SUPERPOP=='AFR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AFR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_EUR <- ggplot(subset(newdata, SUPERPOP=='EUR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EUR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_SAS <- ggplot(subset(newdata, SUPERPOP=='SAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=SAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_EAS <- ggplot(subset(newdata, SUPERPOP=='EAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_AMR <- ggplot(subset(newdata, SUPERPOP=='AMR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AMR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))

p_XTRA_corrected_AFR <- ggplot(subset(newdata, SUPERPOP=='AFR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AFR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_corrected_EUR <- ggplot(subset(newdata, SUPERPOP=='EUR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EUR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_corrected_SAS <- ggplot(subset(newdata, SUPERPOP=='SAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=SAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_corrected_EAS <- ggplot(subset(newdata, SUPERPOP=='EAS'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=EAS_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))
p_XTRA_corrected_AMR <- ggplot(subset(newdata, SUPERPOP=='AMR'), aes(x=distance)) + theme_pubr() +
  geom_col(aes(y=XTR_A_corrected, fill=POP), position='dodge', color='black') +
  coord_cartesian(ylim=c(-0.12,0.02)) + scale_fill_manual("", values=AMR_color) +
  theme(legend.position="none") +
  labs( y=bquote(bold("Difference in Relative"~"XTR"[pi]*"/"*"A"[pi]))) + 
  theme(axis.text=element_text(size=12), axis.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'mm'))

# create grobs for SUPERPOP labels
a <- text_grob("Africa", size=18, face='bold')
b <- text_grob("Europe", size=18, face='bold')
c <- text_grob("S. Asia", size=18, face='bold')
d <- text_grob("E. Asia", size=18, face='bold')
e <- text_grob("Amer.", size=18, face='bold')
gt <- gtable_row('name-row', widths = unit(c(.2,.2,.2,.2,.2), 'npc'), grobs=list(a,b,c,d,e))

#create grobs for axis labels
a1 <- text_grob(bquote(bold("Difference in"~"X"[pi]*"/"*"A"[pi])), size=16, face='bold', rot=90)
b1 <- text_grob(bquote(bold("Difference in"~"PAR"[pi]*"/"*"A"[pi])), size=16, face='bold', rot=90)
c1 <- text_grob(bquote(bold("Difference in"~"XTR"[pi]*"/"*"A"[pi])), size=16, face='bold', rot=90)
gt1 <- gtable_col('axis-col1', heights = unit(c(1/3,1/3,1/3), 'npc'), grobs=list(a1,b1,c1))

a2 <- text_grob(bquote(bold("Difference in Relative"~"X"[pi]*"/"*"A"[pi])), size=14, face='bold', rot=90)
b2 <- text_grob(bquote(bold("Difference in Relative"~"PAR"[pi]*"/"*"A"[pi])), size=14, face='bold', rot=90)
c2 <- text_grob(bquote(bold("Difference in Relative"~"XTR"[pi]*"/"*"A"[pi])), size=14, face='bold', rot=90)
gt2 <- gtable_col('axis-col1', heights = unit(c(1/3,1/3,1/3), 'npc'), grobs=list(a2,b2,c2))

blank <- rectGrob(gp=gpar(col="white"))
xAxis <- text_grob("Distance Filtered from Genes", size=20, face='bold')

p1_temp <- ggarrange(p_XA_AFR, p_XA_EUR, p_XA_SAS, p_XA_EAS, p_XA_AMR,
                    p_PARA_AFR, p_PARA_EUR, p_PARA_SAS, p_PARA_EAS, p_PARA_AMR,
                    p_XTRA_AFR, p_XTRA_EUR, p_XTRA_SAS, p_XTRA_EAS, p_XTRA_AMR,
                    ncol=5, nrow=3, legend='none', align='hv')
p1 <- ggarrange(blank, as_ggplot(gt), 
                as_ggplot(gt1), p1_temp,
                blank, xAxis,
                nrow=3, ncol=2, heights=c(0.05,1,0.05), widths =c(0.05,1))

p2_temp <- ggarrange(p_XA_corrected_AFR, p_XA_corrected_EUR, p_XA_corrected_SAS, p_XA_corrected_EAS, p_XA_corrected_AMR,
                     p_PARA_corrected_AFR, p_PARA_corrected_EUR, p_PARA_corrected_SAS, p_PARA_corrected_EAS, p_PARA_corrected_AMR,
                     p_XTRA_corrected_AFR, p_XTRA_corrected_EUR, p_XTRA_corrected_SAS, p_XTRA_corrected_EAS, p_XTRA_corrected_AMR,
                     ncol=5, nrow=3, legend='none', align='hv')
p2 <- ggarrange(blank, as_ggplot(gt), 
                as_ggplot(gt2), p2_temp,
                blank, xAxis,
                nrow=3, ncol=2, heights=c(0.05,1,0.05), widths =c(0.05,1))

ggsave(p1, file=opt$output1, height=opt$height, width=opt$width, units=opt$units)
ggsave(p2, file=opt$output2, height=opt$height, width=opt$width, units=opt$units)
