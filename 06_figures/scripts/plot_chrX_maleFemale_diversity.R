suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('--males'), type='character', default=NULL,
              help="path to males input file"),
  make_option(c('--females'), type='character', default=NULL,
              help="path to females input file"),
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

if(is.null(opt$males) || is.null(opt$females) || is.null(opt$output)) {
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
df_males <- read.delim(opt$males, header=FALSE)
df_males$sex <- 'males'
df_males$position <- sapply(df_males$V2, function(x){((x+(x+100000))/2)/win_size})
df_males$pi <- ifelse(df_males$V5 == 'NA', NA, df_males$V4)

df_females <- read.delim(opt$females, header=FALSE)
df_females$sex <- 'females'
df_females$position <- sapply(df_females$V2, function(x){((x+(x+100000))/2)/win_size})
df_females$pi <- ifelse(df_females$V5 == 'NA', NA, df_females$V4)

df <- rbind(df_females, df_males)

# create the ggplot
p1 <- ggplot(df, aes(x=position, y=pi, colour=sex))

p1 <- p1 + geom_point()

p1 <- p1 + labs(list(x=paste('Position (', win_type, ')'), 
                     y=expression(paste('Diversity (', pi, ')')))) +
  theme(axis.title=element_text(size=14)) +
  scale_colour_manual("Sex", values=c("males"="red","females"="black"))

p1 <- p1 + geom_vline(aes(xintercept = 2.699520), linetype="dashed") +
  geom_vline(aes(xintercept = 88.193855), linetype="dashed") +
  geom_vline(aes(xintercept = 93.193855), linetype="dashed") +
  geom_vline(aes(xintercept = 154.931044), linetype="dashed")

p1 <- p1 + ylim(0,0.0045) # keep scales standard for all plots

ggsave(file=opt$output, height=opt$height, width=opt$width, units=opt$units)
