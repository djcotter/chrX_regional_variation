suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c('--input'), type='character', default=NULL,
              help="path to input file"),
  make_option(c('--output'), type='character', default=NULL,
              help="path to output file"),
  make_option(c('--width'), type='double', default=4.8,
              help='width for figure'),
  make_option(c('--height'), type='double', default=10.0,
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

df <- read.table(opt$input, header=T)

df_summ <- df %>% 
  select(gen, pop, region, pi) %>% 
  group_by(gen, pop, region) %>% 
  summarise(mean_pi = mean(pi), pi_l = quantile(pi, 0.025), pi_h = quantile(pi, 0.975))

df_ratio <- df %>% 
  select(gen, pop, region, replicate, pi) %>% 
  spread(region, pi) %>% 
  mutate(r1_r2=r1/r2) %>% 
  select(-replicate) %>% 
  group_by(gen, pop) %>% 
  summarise(mean_pi = mean(r1_r2), pi_l = quantile(r1_r2, 0.025), pi_h = quantile(r1_r2, 0.975))

p1 <- ggplot(df_summ %>% subset(region=='r1'), aes(x=gen, y=mean_pi, color=pop)) + geom_line() + 
  geom_ribbon(aes(ymin=pi_l, ymax=pi_h), lty='dashed', alpha=0.1, show.legend = F) + theme_pubr() + coord_cartesian(ylim=c(0,0.0015)) +
  xlab('Generation') + ylab('Average genetic diversity\n(Region 1)') + labs(color="Population") +
  scale_color_discrete(labels=c('AFR', 'EUR'), 
                       guide = guide_legend(override.aes = list(size=2)))

p2 <- ggplot(df_summ %>% subset(region=='r2'), aes(x=gen, y=mean_pi, color=pop)) + geom_line() + 
  geom_ribbon(aes(ymin=pi_l, ymax=pi_h), lty='dashed', alpha=0.1, show.legend = F) + theme_pubr() + coord_cartesian(ylim=c(0,0.0015)) +
  xlab('Generation') + ylab('Average genetic diversity\n(Region 2)') + labs(color="Population") +
  scale_color_discrete(labels=c('AFR', 'EUR'), 
                       guide = guide_legend(override.aes = list(size=2)))

p3 <- ggplot(df_ratio, aes(x=gen, y=mean_pi, color=pop)) + geom_line() + 
  geom_ribbon(aes(ymin=pi_l, ymax=pi_h), lty='dashed', alpha=0.1, show.legend=F) + theme_pubr() +
  coord_cartesian(ylim=c(0.75,1.75)) + geom_hline(yintercept=1, lty='dashed') +
  xlab('Generation') + ylab('Average genetic diversity\n(Region 1/Region 2)') + labs(color="Population") +
  scale_color_discrete(labels=c('AFR', 'EUR'), 
                       guide = guide_legend(override.aes = list(size=2)))

p <- ggarrange(p1,p2,p3, ncol=1, common.legend=T, labels = c('a)', 'b)', 'c)'))

ggsave(p, filename=opt$output, height=opt$height, width=opt$width, units=opt$units)
