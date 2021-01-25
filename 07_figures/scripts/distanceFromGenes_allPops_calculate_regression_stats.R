suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))
library(broom)


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
  make_option(c('--output1'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--denom_pop'), type='character', default='MSL',
              help="Population code used for demography normalization")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$units = tolower(opt$units)

if( is.null(opt$filter1) || is.null(opt$filter2) || is.null(opt$filter3) || is.null(opt$filter4) || is.null(opt$filter5) ) {
  print_help(opt_parser)
  stop("Input files must be specified.", call.=FALSE)
}

if( is.null(opt$output)) {
  print_help(opt_parser)
  stop("output file must be specified.", call.=FALSE)
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

pop_list <- c("LWK", "GWD", "ACB", "ASW", "MSL", "YRI", "ESN", 
              "CEU", "IBS", "TSI", "FIN", "GBR", 
              "GIH", "BEB", "ITU", "STU", "PJL", 
              "JPT", "CDX", "CHB", "CHS", "KHV", 
              "PUR", "CLM", "MXL", "PEL")

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

denom_pop = opt$denom_pop

data_filter1$XA <- data_filter1$nonPAR / data_filter1$chr8
data_filter1$XA_corrected <- data_filter1$XA / subset(data_filter1,POP==denom_pop)$XA
data_filter1$PARA <- data_filter1$PAR1 / data_filter1$chr8
data_filter1$PARA_corrected <- data_filter1$PARA / subset(data_filter1,POP==denom_pop)$PARA
data_filter1$XTRA <- data_filter1$XTR / data_filter1$chr8
data_filter1$XTRA_corrected <- data_filter1$XTRA / subset(data_filter1,POP==denom_pop)$XTRA
data_filter1$distance <- '0'
data_filter1$distance_num <- 0

data_filter2$XA <- data_filter2$nonPAR / data_filter2$chr8
data_filter2$XA_corrected <- data_filter2$XA / subset(data_filter2,POP==denom_pop)$XA
data_filter2$PARA <- data_filter2$PAR1 / data_filter2$chr8
data_filter2$PARA_corrected <- data_filter2$PARA / subset(data_filter2,POP==denom_pop)$PARA
data_filter2$XTRA <- data_filter2$XTR / data_filter2$chr8
data_filter2$XTRA_corrected <- data_filter2$XTRA / subset(data_filter2,POP==denom_pop)$XTRA

data_filter2$distance <- '1'
data_filter2$distance_num <- 1

data_filter3$XA <- data_filter3$nonPAR / data_filter3$chr8
data_filter3$XA_corrected <- data_filter3$XA / subset(data_filter3,POP==denom_pop)$XA
data_filter3$PARA <- data_filter3$PAR1 / data_filter3$chr8
data_filter3$PARA_corrected <- data_filter3$PARA / subset(data_filter3,POP==denom_pop)$PARA
data_filter3$XTRA <- data_filter3$XTR / data_filter3$chr8
data_filter3$XTRA_corrected <- data_filter3$XTRA / subset(data_filter3,POP==denom_pop)$XTRA

data_filter3$distance <- '5'
data_filter3$distance_num <- 5

data_filter4$XA <- data_filter4$nonPAR / data_filter4$chr8
data_filter4$XA_corrected <- data_filter4$XA / subset(data_filter4,POP==denom_pop)$XA
data_filter4$PARA <- data_filter4$PAR1 / data_filter4$chr8
data_filter4$PARA_corrected <- data_filter4$PARA / subset(data_filter4,POP==denom_pop)$PARA
data_filter4$XTRA <- data_filter4$XTR / data_filter4$chr8
data_filter4$XTRA_corrected <- data_filter4$XTRA / subset(data_filter4,POP==denom_pop)$XTRA

data_filter4$distance <- '10'
data_filter4$distance_num <- 10

data_filter5$XA <- data_filter5$nonPAR / data_filter5$chr8
data_filter5$XA_corrected <- data_filter5$XA / subset(data_filter5,POP==denom_pop)$XA
data_filter5$PARA <- data_filter5$PAR1 / data_filter5$chr8
data_filter5$PARA_corrected <- data_filter5$PARA / subset(data_filter5,POP==denom_pop)$PARA
data_filter5$XTRA <- data_filter5$XTR / data_filter5$chr8
data_filter5$XTRA_corrected <- data_filter5$XTRA / subset(data_filter5,POP==denom_pop)$XTRA

data_filter5$distance <- '20'
data_filter5$distance_num <- 20

data <- rbind(data_filter1, data_filter2, data_filter3, data_filter4, data_filter5)

myDf <- NULL
for (sup in unique(data$SUPERPOP)) {
  subset_data <- data %>% subset(SUPERPOP == sup) %>% 
    select(SUPERPOP, POP, 
           XA, XA_corrected,
           PARA, PARA_corrected, 
           XTRA, XTRA_corrected, distance_num) %>%
    gather(ratio, value, -SUPERPOP, -POP, -distance_num) %>%
    separate(ratio, sep="_", into=c('ratio', 'correction'), fill = 'right') %>%
    mutate(correction=ifelse(is.na(correction), 'uncorrected', correction))
  for (reg in unique(subset_data$ratio)) {
    for (corr in c('uncorrected', 'corrected')) {
      myData <- subset_data %>% subset(ratio == reg) %>%
        subset(correction == corr)
      lm <- lm(value ~ distance_num, data = myData)
      r.sq <- summary(lm)$coefficients[2,1]
      p.val <- summary(lm)$coefficients[2,4]
      temp_df <- data.frame(SUPERPOP=sup, Ratio=reg, Normalization=corr, r.sq=r.sq, p.val=p.val)
      myDf <- rbind(myDf, temp_df)
    }
  }
}

write.table(myDf, file = opt$output)
