suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

max_height = 0.0045

option_list = list(
  make_option(c('--chrX_males'), type='character', default=NULL,
              help="path to chrX males input file"),
  make_option(c('--chrX_females'), type='character', default=NULL,
              help="path to chrX females input file"),
  make_option(c('--chrY'), type='character', default=NULL,
              help="path to chrY input file"),
  make_option(c('-o', '--output'), type='character', default=NULL,
              help="path to output file"),
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

if(is.null(opt$chrX_males) || is.null(opt$chrX_females) || is.null(opt$chrY) || is.null(opt$output)) {
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

# declare window size for use plotting values
win_size <- 1000000 # 1 Mb
win_type <- 'Mb' # 1 Mb

# data frame should consist of cols: chr, start, end, value, called_sites, snp_count
df_Xmales <- read.delim(opt$chrX_males, header=FALSE, colClasses = c(NA,NA,NA,NA,NA,NA))
df_Xmales$sex <- 'Males'
df_Xmales$chr <- 'chrX'
df_Xmales$position <- sapply(df_Xmales$V2, function(x){((x+(x+100000))/2)/win_size})
df_Xmales$pi <- ifelse(df_Xmales$V5 == 'NA', NA, df_Xmales$V4)

df_Xfemales <- read.delim(opt$chrX_females, header=FALSE, colClasses = c(NA,NA,NA,NA,NA,NA))
df_Xfemales$sex <- 'Females'
df_Xfemales$chr <- 'chrX'
df_Xfemales$position <- sapply(df_Xfemales$V2, function(x){((x+(x+100000))/2)/win_size})
df_Xfemales$pi <- ifelse(df_Xfemales$V5 == 'NA', NA, df_Xfemales$V4)

df_Y <- read.delim(opt$chrY, header=FALSE, colClasses = c(NA,NA,NA,NA,NA,NA))
df_Y$sex <- 'Males'
df_Y$chr <- 'chrY'
df_Y$position <- sapply(df_Y$V2, function(x){((x+(x+100000))/2)/win_size})
df_Y$pi <- ifelse(df_Y$V5 == 'NA', NA, df_Y$V4)

df <- rbind(df_Xfemales, df_Xmales)
df <- rbind(df, df_Y)
df <- subset(df, position < 2.699520 + 2.50000)

# create the ggplot
p1 <- ggplot(df, aes(x=position, y=pi, colour=sex, shape=chr))

p1 <- p1 + geom_point()

p1 <- p1 + labs(list(x=paste('Position (', win_type, ')'), 
                     y=expression(paste('Diversity (', pi, ')')))) +
  theme(axis.title=element_text(size=14)) +
  scale_colour_manual("Sex", values=c("Males"="red","Females"="black"))

p1 <- p1 + geom_vline(aes(xintercept = 2.699520), linetype="dashed")

p1 <- p1 + ylim(0,max_height) # keep scales standard for all plots

ggsave(file=opt$output, height=opt$height, width=opt$width, units=opt$units)
